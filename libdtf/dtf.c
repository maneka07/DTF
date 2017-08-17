#include <mpi.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stddef.h>
#include <unistd.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>

#include "dtf.h"
#include "dtf_init_finalize.h"
#include "dtf_util.h"
#include "dtf_nbuf_io.h"
#include "dtf_req_match.h"

int lib_initialized=0;
int gl_verbose;
int gl_my_rank;
struct dtf_config gl_conf;
struct stats gl_stats;

void *gl_msg_buf = NULL;
extern file_info_req_q_t *gl_finfo_req_q;
extern dtf_msg_t *gl_msg_q;

/*
************************************User API**********************************************
*/
/**
  @brief	Function to initialize the library. Should be called from inside
            the application before any other call to the library. Should be called
            after the MPI is initialized.
  @param	filename        Name of the library configuration file.
  @param    module_name     Name of the module calling the init function.
  @return	int             0 if OK, anything else otherwise

 */
_EXTERN_C_ int dtf_init(const char *filename, char *module_name)
{
  //  char* conf_filepath;
    int err, mpi_initialized;
    char* s;
    int verbose;

    if(lib_initialized)
        return 0;

    MPI_Initialized(&mpi_initialized);

    if(!mpi_initialized){
        fprintf(stderr, "DTF Error: dtf cannot be initialized before MPI is initialized. Aborting...\n");
        fflush(stdout);
        fflush(stderr);
        exit(1);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &gl_my_rank);

    if(strlen(module_name)>MAX_COMP_NAME){
        fprintf(stderr, "DTF Error: module name %s too long\n", module_name);
        fflush(stderr);
        exit(1);
    }
    gl_stats.malloc_size = 0;
    gl_stats.data_msg_sz = 0;
    gl_stats.ndata_msg_sent = 0;
    gl_stats.accum_match_time = 0;
    gl_stats.ndb_match = 0;
    gl_stats.walltime = MPI_Wtime();
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: started at %.3f", gl_stats.walltime);
    gl_stats.accum_comm_time = 0;
    gl_stats.nprogress_call = 0;
    gl_stats.nioreqs = 0;
    gl_stats.nbl = 0;
    gl_stats.ngetputcall = 0;
    gl_stats.timer_accum = 0;
    gl_stats.timer_start = 0;
    gl_stats.accum_comm_data_time = 0;
    gl_stats.accum_extract_data_time = 0;
    gl_stats.accum_dbuff_sz = 0;
    gl_stats.accum_dbuff_time = 0;
    gl_stats.accum_rw_var = 0;
    gl_stats.accum_progr_time = 0;
    gl_stats.accum_do_matching_time = 0;
    gl_stats.nfiles = 0;
    gl_stats.accum_send_ioreq_time = 0;
    gl_stats.progr_work_time = 0;
    gl_stats.do_match_idle_time = 0;

    gl_my_comp_name = (char*)dtf_malloc(MAX_COMP_NAME);
    assert(gl_my_comp_name != NULL);
    strcpy(gl_my_comp_name, module_name);

    s = getenv("DTF_VERBOSE_LEVEL");
    if(s == NULL)
        gl_verbose = VERBOSE_DBG_LEVEL;
    else
        gl_verbose = atoi(s);

    s = getenv("DTF_DETECT_OVERLAP");
    if(s == NULL)
        gl_conf.detect_overlap_flag = 0;
    else
        gl_conf.detect_overlap_flag = atoi(s);

    s = getenv("DTF_DATA_MSG_SIZE_LIMIT");
    if(s == NULL)
        gl_conf.data_msg_size_limit = DTF_DATA_MSG_SIZE_LIMIT;
    else
        gl_conf.data_msg_size_limit = atoi(s) * 1024;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Data message size limit set to %d", gl_conf.data_msg_size_limit);

    assert(gl_conf.data_msg_size_limit > 0);

    gl_msg_buf = NULL;

    //during init only root will print out stuff
    if(gl_my_rank != 0){
        verbose = gl_verbose;
        gl_verbose = VERBOSE_ERROR_LEVEL;
    }
    /*Parse ini file and initialize components*/
    err = load_config(filename, module_name);
    if(err) goto panic_exit;

    /*Establish intercommunicators between components*/
    err = init_comp_comm();
    if(err) goto panic_exit;

    lib_initialized = 1;

    //enable print setting for other ranks again
    if(gl_my_rank != 0)
        gl_verbose = verbose;

    DTF_DBG(VERBOSE_DBG_LEVEL, "DTF: Finished initializing");

    return 0;

panic_exit:

    dtf_finalize();
    exit(1);
    return 1;
}

/**
  @brief	Function to finalize the library. Should be called from inside
            the application before the MPI is finalized.
  @return	int             0 if OK, anything else otherwise

 */
_EXTERN_C_ int dtf_finalize()
{
    int mpi_initialized;
    int err, nranks;
    int intsum;
    char *s;
    double dblsum = 0, walltime;
    unsigned long data_sz, lngsum;
    int sclltkf = 0;
    unsigned unsgn;
    struct{
        double dbl;
        int intg;
    }dblint_in, dblint_out;

    if(!lib_initialized) return 0;

    MPI_Initialized(&mpi_initialized);

    if(!mpi_initialized){
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: dtf cannot be finalized after MPI is finalized. Aborting...");
        fflush(stdout);
        fflush(stderr);
        exit(1);
    }

//int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
//               MPI_Op op, int root, MPI_Comm comm)

    while(gl_msg_q != NULL)
        progress_msg_queue();

    /*Send any unsent file notifications
      and delete buf files*/
    finalize_files();
    assert(gl_finfo_req_q == NULL);

    finalize_comp_comm();

    clean_config();

    DTF_DBG(VERBOSE_DBG_LEVEL,"DTF: finalize");
    walltime = MPI_Wtime() - gl_stats.walltime;
    MPI_Comm_size(gl_comps[gl_my_comp_id].comm, &nranks);
    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: walltime %.3f", walltime);


    if(gl_stats.timer_accum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: library-measured I/O time: %.4f", gl_stats.timer_accum);
    if(gl_stats.accum_comm_time > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: time spent in comm: %.5f: %.4f",
                gl_stats.accum_comm_time, (gl_stats.accum_comm_time/gl_stats.timer_accum)*100);

    if(gl_stats.accum_match_time > 0){
         DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: tot match time: %.5f: %.5f",
                 gl_stats.accum_match_time, (gl_stats.accum_match_time/gl_stats.timer_accum)*100);
         DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: do match time: %.5f: %.5f",
                 gl_stats.accum_do_matching_time, (gl_stats.accum_do_matching_time/gl_stats.timer_accum)*100);
         DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: idle do match time: %.5f",
                 gl_stats.do_match_idle_time);
    }
    if(gl_stats.accum_progr_time > 0){
        DTF_DBG(VERBOSE_DBG_LEVEL, "DTF STAT: progr time: %.4f: %.4f", gl_stats.accum_progr_time, (gl_stats.accum_progr_time/gl_stats.timer_accum)*100);
        DTF_DBG(VERBOSE_DBG_LEVEL, "DTF STAT: progr work time: %.4f", gl_stats.progr_work_time);
    }

    if(gl_stats.accum_send_ioreq_time > 0)
        DTF_DBG(VERBOSE_DBG_LEVEL, "DTF STAT: send ioreqs time: %.4f: %.4f", gl_stats.accum_send_ioreq_time, (gl_stats.accum_send_ioreq_time/gl_stats.timer_accum)*100);

    if(gl_stats.accum_rw_var > 0)
        DTF_DBG(VERBOSE_DBG_LEVEL, "DTF STAT: rw var time: %.4f: %.4f", gl_stats.accum_rw_var, (gl_stats.accum_rw_var/gl_stats.timer_accum)*100);

//    if(gl_stats.accum_match_time > 0){
//          DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: db match time: %.5f: %.5f",
//                 gl_stats.accum_db_match_time, (gl_stats.accum_db_match_time/gl_stats.timer_accum)*100);
//         DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: db manage time: %.5f: %.5f",
//                 gl_stats.accum_db_manage_time,  (gl_stats.accum_db_manage_time/gl_stats.timer_accum)*100);
//
//    }
    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: times matched: %u", gl_stats.ndb_match);
    if(gl_stats.nprogress_call > 0)
        DTF_DBG(VERBOSE_DBG_LEVEL, "DTF STAT: times dbmatched %d, progress call %lu", gl_stats.ndb_match, gl_stats.nprogress_call);

    if(gl_stats.nioreqs > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: nreqs: %u", gl_stats.nioreqs);


//    if(gl_stats.timer2_accum > 0)
//        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: user-program-measured I/O time: %.4f", gl_stats.timer2_accum);
    if(gl_stats.nbl > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: Out of %u blocks handled %u noncontig blocks", gl_stats.nbl, gl_stats.ngetputcall);
    if(gl_stats.accum_dbuff_sz > 0){
        DTF_DBG(VERBOSE_DBG_LEVEL, "DTF STAT: buffering time: %.5f: %.4f", gl_stats.accum_dbuff_time,(gl_stats.accum_dbuff_time/gl_stats.timer_accum)*100);
        DTF_DBG(VERBOSE_DBG_LEVEL, "DTF STAT: buffering size: %lu",  gl_stats.accum_dbuff_sz);
    }
    if(gl_stats.ndata_msg_sent > 0){
        data_sz = (unsigned long)(gl_stats.data_msg_sz/gl_stats.ndata_msg_sent);
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: total sent %lu, avg data msg sz %lu (%d msgs)", gl_stats.data_msg_sz, data_sz, gl_stats.ndata_msg_sent);
    } else
        data_sz = 0;

    /*In scale-letkf the last ranks write mean files not treated by dtf, we don't need
      stats for that*/
    s = getenv("DTF_SCALE");
    if(s != NULL)
        sclltkf = atoi(s);
    else
        sclltkf = 0;

    if(sclltkf){
        nranks = nranks - (int)(nranks % (gl_stats.nfiles/2));
        if(gl_my_rank == 0)
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF: for stats consider %d ranks", nranks);
        if(gl_my_rank >= nranks){//Only for because the last set of processes work with mean files
            DTF_DBG(VERBOSE_DBG_LEVEL, "will zero lib time");
            gl_stats.timer_accum = 0;
            walltime = 0;
            gl_stats.accum_match_time = 0;
            gl_stats.accum_comm_data_time = 0;
            gl_stats.accum_comm_time = 0;
            gl_stats.accum_dbuff_sz = 0;
            gl_stats.accum_dbuff_time = 0;
            gl_stats.accum_progr_time = 0;
            gl_stats.accum_rw_var = 0;
            gl_stats.progr_work_time = 0;
            gl_stats.nioreqs = 0;
            gl_stats.progr_work_time = 0;
            gl_stats.do_match_idle_time = 0;
        }
    }

    /*AVERAGE STATS*/
    dblint_in.dbl = walltime;
    dblint_in.intg = gl_my_rank;
    err = MPI_Reduce(&dblint_in, &dblint_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: max walltime: %.4f: rank: %d", dblint_out.dbl, dblint_out.intg);

    dblint_in.dbl = gl_stats.timer_accum;
    dblint_in.intg = gl_my_rank;
    err = MPI_Reduce(&dblint_in, &dblint_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: max libtime: %.4f: rank: %d", dblint_out.dbl, dblint_out.intg);

    err = MPI_Reduce(&(gl_stats.nioreqs), &unsgn, 1, MPI_UNSIGNED, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg nioreqs: %u", (unsigned)(unsgn/nranks));

    err = MPI_Reduce(&gl_stats.nioreqs, &unsgn, 1, MPI_UNSIGNED, MPI_MAX, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: max nioreqs: %u", unsgn);

    err = MPI_Reduce(&walltime, &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg walltime: %.4f", dblsum/nranks);

    err = MPI_Allreduce(&(gl_stats.timer_accum), &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank==0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg library-measured I/O time: %.5f", dblsum/nranks);

    /*Standard deviation*/
    {
        double mydev2;
        double mean = dblsum/nranks;
        mydev2 = (gl_stats.timer_accum - mean)*(gl_stats.timer_accum - mean);

        err = MPI_Reduce(&mydev2, &dblsum, 1, MPI_DOUBLE, MPI_SUM,0, gl_comps[gl_my_comp_id].comm);
        CHECK_MPI(err);

        if(gl_my_rank == 0)
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: standard deviation: %.7f", sqrt(dblsum/nranks));
    }

    err = MPI_Reduce(&(gl_stats.accum_match_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg tot match time: %.4f", dblsum/nranks);

    err = MPI_Reduce(&(gl_stats.accum_do_matching_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg do match time: %.4f", dblsum/nranks);

    err = MPI_Reduce(&(gl_stats.accum_send_ioreq_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg send ioreqs time: %.4f", dblsum/nranks);

    err = MPI_Reduce(&(gl_stats.do_match_idle_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg do match idle time: %.4f", dblsum/nranks);


    err = MPI_Reduce(&(gl_stats.accum_progr_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg progr time: %.4f", dblsum/nranks);

    err = MPI_Reduce(&(gl_stats.progr_work_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg work progr time: %.4f", dblsum/nranks);


    err = MPI_Reduce(&(gl_stats.accum_rw_var), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: rw var time: %.4f", dblsum/nranks);

//    err = MPI_Reduce(&(gl_stats.accum_db_match_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
//    CHECK_MPI(err);
//    if(gl_my_rank == 0 && dblsum > 0)
//        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg db match time: %.5f", dblsum/nranks);

    err = MPI_Reduce(&(gl_stats.accum_comm_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum>0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg comm time: %.5f", dblsum/nranks);

    err = MPI_Reduce(&(gl_stats.accum_extract_data_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank==0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: Avg extract data time: %.3f", (double)(dblsum/nranks));

   err = MPI_Reduce(&(gl_stats.accum_dbuff_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank==0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg buffering time: %.5f", (double)(dblsum/nranks));

    err = MPI_Reduce(&data_sz, &lngsum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank==0 && lngsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: Avg data msg sz acrosss ps: %d", (int)(lngsum/nranks));

//    data_sz = (unsigned long) gl_stats.data_msg_sz;
//    err = MPI_Reduce(&data_sz, &lngsum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
//    CHECK_MPI(err);
//    if(gl_my_rank==0 && lngsum > 0)
//        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: Avg total data sent acrosss ps: %d", (int)(lngsum/nranks));
    intsum = 0;
    err = MPI_Reduce(&(gl_stats.ndata_msg_sent), &intsum, 1, MPI_INT, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank==0 && intsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: Avg num data msg: %d", (int)(intsum/nranks));

    err = MPI_Reduce(&(gl_stats.accum_dbuff_sz), &lngsum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank==0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg buffering size: %lu", (size_t)(lngsum/nranks));

    if(gl_msg_buf != NULL)
        dtf_free(gl_msg_buf, gl_conf.data_msg_size_limit);

    if(gl_stats.malloc_size != MAX_COMP_NAME )
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: DTF memory leak size: %lu", gl_stats.malloc_size - MAX_COMP_NAME);
    assert(gl_finfo_req_q == NULL);
    assert(gl_msg_q == NULL);

    dtf_free(gl_my_comp_name, MAX_COMP_NAME);
    lib_initialized = 0;
    fflush(stdout);
    fflush(stderr);
    return 0;
}

/*called by user to do explicit matching*/
/*
    User must specify either filename or ncid.
    intracomp_io_flag - if set to 1, matching of intracomponent io requests will be
    performed. This flag is intended for for situation when the writer component
    tries to read something from the file it is writing.
*/
_EXTERN_C_ int dtf_match_io(const char *filename, int ncid, int intracomp_io_flag )//, int match_all)
{
    file_buffer_t *fbuf;
    if(!lib_initialized) return 0;
    DTF_DBG(VERBOSE_DBG_LEVEL, "call match io for %s (ncid %d), intra flag %d", filename, ncid, intracomp_io_flag);
    if(intracomp_io_flag){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: scale-letkf hack: skip intracomp matching");
        return 0;
    }

    fbuf = find_file_buffer(gl_filebuf_list, filename, ncid);
    if(fbuf == NULL){

        if( (filename != NULL) && (strlen(filename) == 0) )
            DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: file (%s) with ncid %d is not treated by DTF (not in configuration file). Explicit matching ignored.", filename, ncid);
        else
            DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: file %s (ncid %d) is not treated by DTF (not in configuration file). Explicit matching ignored.", filename, ncid);
        return 0;
    }
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;
    if(!fbuf->explicit_match){
        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: calling dtf_match_io but explicit match for file %s not enabled. Ignored.", filename);
        return 0;
    }

    if( intracomp_io_flag && (gl_my_comp_id != fbuf->writer_id) ){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: dtf_match_io: intracomp_io_flag(%d) can only be set for the writer component", intracomp_io_flag);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    if(fbuf->is_matching_flag){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_match_io is called for file %s, but a matching process has already started before.", fbuf->file_path);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    dtf_tstart();
    match_ioreqs(fbuf, intracomp_io_flag);
    dtf_tend();
    return 0;
}

/*Supposed to be called by the writer process.
  Used to match against several dtf_match_io functions on the reader side*/
_EXTERN_C_ void dtf_match_multiple(int ncid)
{

    if(!lib_initialized) return;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: ncid %d is not treated by DTF( \
                not in configuration file). Explicit matching ignored.", ncid);
        return;
    }
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;
    if(!fbuf->explicit_match){
        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: calling dtf_match_multiple but explicit match for file \
                %s not enabled. Ignored.", fbuf->file_path);
        return;
    }

    if( gl_my_comp_id != fbuf->writer_id){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: dtf_match_multiple can only be called by writer component. Ignoring.");
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Start matching multiple");
    dtf_tstart();
    fbuf->done_match_multiple_flag = 0;
    while(!fbuf->done_match_multiple_flag)
        match_ioreqs(fbuf, 0);
    //reset
    fbuf->done_match_multiple_flag = 0;
    dtf_tend();
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finish matching multiple");
    return;
}

/*Used by reader to notify writer that it can complete dtf_match_multiple*/
_EXTERN_C_ void dtf_complete_multiple(const char *filename, int ncid)
{
    double t_start;
    if(!lib_initialized) return;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, NULL, ncid);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: file %s (ncid %d) is not treated by DTF (not in configuration file). Explicit matching ignored.", filename, ncid);
        return;
    }

    t_start = MPI_Wtime();
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;
    if(fbuf->reader_id != gl_my_comp_id){
        DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: dtf_complete_multiple can only be called by reader");
        return;
    }
    dtf_tstart();
    MPI_Barrier(fbuf->comm);
    notify_complete_multiple(fbuf);
    dtf_tend();

    gl_stats.accum_match_time += MPI_Wtime() - t_start;
}

/*
*******************************Interfaces to be used by PnetCDF***************************************************
*/
_EXTERN_C_ void dtf_write_hdr(const char *filename, MPI_Offset hdr_sz, void *header)
{
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;


    if(hdr_sz == 0){
        DTF_DBG(VERBOSE_DBG_LEVEL, "Header size for file %s is zero", filename);
        return;
    }
    write_hdr(fbuf, hdr_sz, header);
    return;
}

_EXTERN_C_ MPI_Offset dtf_read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
    if(!lib_initialized) return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;
    return read_hdr_chunk(fbuf, offset, chunk_sz, chunk);
}

_EXTERN_C_ void dtf_create(const char *filename, MPI_Comm comm, int ncid)
{


    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL, "Creating file %s. File is not treated by DTF", filename);
        return;
    } else {
        DTF_DBG(VERBOSE_DBG_LEVEL, "Created file %s (ncid %d)", filename, ncid);
    }
    fbuf->ncid = ncid;
    fbuf->comm = comm;
    int root_mst = gl_my_rank;
    int err = MPI_Bcast(&root_mst, 1, MPI_INT, 0, comm);
    CHECK_MPI(err);
    fbuf->root_writer = root_mst;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Root master for file %s (ncid %d) is %d", filename, ncid, fbuf->root_writer);
    if(gl_my_rank == fbuf->root_writer && gl_my_rank != 0){
        /*Notify global rank 0 that I am the root master for this file*/
        DTF_DBG(VERBOSE_DBG_LEVEL, "Will notify global rank 0 that I am master");
        dtf_msg_t *msg = new_dtf_msg(NULL, 0, ROOT_MST_TAG);
        err = MPI_Isend(fbuf->file_path, (int)(MAX_FILE_NAME), MPI_CHAR, 0, ROOT_MST_TAG, gl_comps[gl_my_comp_id].comm, &(msg->req));
        CHECK_MPI(err);
        ENQUEUE_ITEM(msg, gl_msg_q);
    }

    DTF_DBG(VERBOSE_DBG_LEVEL, "Init masters");
    init_req_match_masters(comm, fbuf->mst_info);

    if(fbuf->iomode == DTF_IO_MODE_MEMORY){
        if(fbuf->mst_info->is_master_flag){
            int nranks;
            MPI_Comm_size(comm, &nranks);

            assert(fbuf->mst_info->iodb == NULL);
            init_iodb(fbuf);
            fbuf->mst_info->nwranks_opened = (unsigned int)nranks;
        }

    } else if(fbuf->iomode == DTF_IO_MODE_FILE){
        fbuf->fready_notify_flag = RDR_HASNT_OPENED;

        /*Create symlink to this file (needed for SCALE-LETKF since
          there is no way to execute the script to perform all the file
          movement in the middle of the execution)*/
        if(strlen(fbuf->slink_name)>0 && fbuf->root_writer==gl_my_rank){
            int err;
            size_t slen1, slen2;
            char *dir = NULL;
            char wdir[MAX_FILE_NAME]="\0";
            char slink[MAX_FILE_NAME]="\0";
            char origfile[MAX_FILE_NAME]="\0";

            dir = getcwd(wdir, MAX_FILE_NAME);
            assert(dir != NULL);
            slen1 = strlen(wdir);
            slen2 = strlen(fbuf->slink_name);
            if(slen1+slen2+1 > MAX_FILE_NAME){
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: symlink name of length %ld exceeds max filename of %d",slen1+slen2+1, MAX_FILE_NAME);
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            }

            sprintf(slink, "%s/%s", wdir, fbuf->slink_name);
            sprintf(origfile, "%s/%s", wdir, fbuf->file_path);

            DTF_DBG(VERBOSE_DBG_LEVEL, "Creating symlink %s to %s (%s)", slink, origfile, fbuf->file_path);
            err = symlink(origfile, slink);
            if(err != 0){
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: error creating symlink %s to %s : %s", slink, origfile,  strerror(errno));

            }
        }
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Exit create");
}

/**
  @brief	Called when the corresponding file is opened.

  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void dtf_open(const char *filename, MPI_Comm comm)
{
    int nranks;
    if(!lib_initialized) return;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL, "Opening file %s. File is not treated by DTF", filename);
        return;
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Opening file %s", filename);
    /*A hack for scale-letkf: if we are trying to open an alias,
    create an empty file, otherwise letkf will crash because
    it won't  find the file*/
    if( (fbuf->iomode == DTF_IO_MODE_MEMORY) && (strlen(fbuf->alias_name) > 0) && (strstr(filename, fbuf->alias_name) !=NULL)){
        MPI_File fh;
        int err;
        err = MPI_File_open(comm, (char*)filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
        CHECK_MPI(err);
        err = MPI_File_close(&fh);
        CHECK_MPI(err);
        DTF_DBG(VERBOSE_DBG_LEVEL, "Created a dummy alias file");
    }
    if(fbuf->comm == MPI_COMM_NULL)
        fbuf->comm = comm;
    else
        assert(fbuf->comm == comm);
    MPI_Comm_size(comm, &nranks);
    if(comm == gl_comps[gl_my_comp_id].comm)
        DTF_DBG(VERBOSE_DBG_LEVEL, "File opened in component's communicator (%d nprocs)", nranks);
    else
        DTF_DBG(VERBOSE_DBG_LEVEL, "File opened in subcommunicator (%d nprocs)", nranks);

    open_file(fbuf, comm);
}

_EXTERN_C_ void dtf_enddef(const char *filename)
{
//    if(!lib_initialized) return;
//    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
//    if(fbuf == NULL) return;
//    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return;
//
//    if((fbuf->writer_id == gl_my_comp_id) && (gl_conf.distr_mode == DISTR_MODE_REQ_MATCH))
//        send_file_info(fbuf);

    //We don't know which reader rank needs the info, so file info will be sent when
    //corresponding reader requests it!!
}

_EXTERN_C_ void dtf_set_ncid(const char *filename, int ncid)
{
    if(!lib_initialized) return;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL, "File %s is not treated by DTF. Will not set ncid", filename);
        return;
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Set ncid of file %s to %d (previos value %d)", filename, ncid, fbuf->ncid);
    fbuf->ncid = ncid;
}

/**
  @brief	Called when the corresponding file is closed. If the file is opened in write mode
            (output file) we will notify the reader that the file is ready to be transfered. If it's opened in read mode
            (input file), we will free the buffer.
            The mode is defined from the dtf configuration file.
  @param	filename        file name for the memory buffer
  @return	void

 */
_EXTERN_C_ void dtf_close(const char* filename)
{
    if(!lib_initialized) return;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Closing file %s", filename);
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL, "File not treated by dtf");
        return;
    }
//    if(fbuf->reader_id == gl_my_comp_id)

    MPI_Barrier(fbuf->comm);

    close_file(fbuf);
}

/*called inside wait function in pnetcdf*/
_EXTERN_C_ int dtf_match_ioreqs(const char* filename)
{
    int ret;
    if(!lib_initialized) return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;
    /*User will have to explicitly initiate matching*/
    if(fbuf->explicit_match) return 0;

    if(fbuf->is_matching_flag){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_match_ioreqs is called for file %s, but a matching process has already started before.", fbuf->file_path);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    fbuf->is_matching_flag = 1;
    ret = match_ioreqs(fbuf, 0);
    fbuf->is_matching_flag = 0;

    return ret;
}



///* The user has to state whether the process needs to match all read or write requests.
//   Because the process of matching for a reader and writer is not the same. */
//_EXTERN_C_ void dtf_match_io_all(int rw_flag)
//{
//    file_buffer_t *fbuf = gl_filebuf_list;
//
//    if(!lib_initialized) return;
//    if(gl_conf.distr_mode != DISTR_MODE_REQ_MATCH) return;
//
//    if(rw_flag == DTF_READ){
//        DTF_DBG(VERBOSE_WARNING_LEVEL, "dtf_match_io_all() cannot be used in reader component. Ignoring.");
//        return;
//    }
//    /*First check that no matching is already happening*/
//    while(fbuf != NULL){
//        if(fbuf->is_matching_flag){
//            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dtf_match_io_all is called, but a matching process has already started before.");
//            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
//        }
//        fbuf = fbuf->next;
//    }
//
//    match_ioreqs_all(rw_flag);
//
//    return;
//}


_EXTERN_C_ void dtf_print_data(int varid, int dtype, int ndims, MPI_Offset* count, void* data)
{
    return;
//
//    int i, nelems=1, max_print;
//    if(count == NULL)
//        return;
//    for(i = 0; i < ndims; i++)
//        nelems*=count[i];
//
//    if(nelems < 20)
//        max_print = nelems;
//    else
//        max_print = 20;
//    DTF_DBG(VERBOSE_ERROR_LEVEL, "Data for var %d:", varid);
//    for(i = 0; i < max_print; i++)
//        if(dtype == 0)
//            printf("%.3f\t", ((float*)data)[i]);
//        else if(dtype == 1)
//            printf("%.3f\t", ((double*)data)[i]);
//    printf("\n");
}


_EXTERN_C_ void dtf_print(const char *str)
{
    DTF_DBG(VERBOSE_ERROR_LEVEL, "%s", str);
}

_EXTERN_C_ MPI_Offset dtf_read_write_var(const char *filename,
                                          int varid,
                                          const MPI_Offset *start,
                                          const MPI_Offset *count,
                                          const MPI_Offset *stride,
                                          const MPI_Offset *imap,
                                          MPI_Datatype dtype,
                                          void *buf,
                                          int rw_flag,
                                          int *request)
{
    MPI_Offset ret;

    if(!lib_initialized) return 0;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY)
        return 0;

    if(boundary_check(fbuf, varid, start, count ))
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    if( rw_flag != DTF_READ && rw_flag != DTF_WRITE){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: rw_flag value incorrect (%d)", rw_flag);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    if(rw_flag==DTF_WRITE && fbuf->reader_id == gl_my_comp_id){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: reader component cannot write to the file %s", fbuf->file_path);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }

     if(request == NULL)
        ret = nbuf_read_write_var(fbuf, varid, start, count, stride, imap, dtype, buf, rw_flag, NULL);
     else
        ret = nbuf_read_write_var(fbuf, varid, start, count, stride, imap, dtype, buf, rw_flag, request);

    return ret;
}


/**
    @brief  Returns the I/O mode for this file
    @param  filename    name of the file
    @return the io mode
*/

_EXTERN_C_ int dtf_io_mode(const char* filename)
{
    if(!lib_initialized) return 0;
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL){
        return DTF_IO_MODE_UNDEFINED;
    }
    return fbuf->iomode;
}

_EXTERN_C_ int dtf_def_var(const char* filename, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int i, ret;
    if(!lib_initialized) return 0;
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename, -1);
    if(fbuf == NULL) return 0;
    if(fbuf->iomode != DTF_IO_MODE_MEMORY) return 0;

    DTF_DBG(VERBOSE_ERROR_LEVEL, "def var %d for ncid %d", varid, fbuf->ncid);

    /*For now, can only support unlimited dimension if it's the fasted changing dimension array*/
    if( (ndims > 0) && (shape[ndims - 1] == DTF_UNLIMITED))
        DTF_DBG(VERBOSE_DBG_LEVEL, "var has unlimited dimension");

    for(i = 1; i < ndims-1; i++){
        //we can support unlimited dimension if it's the first dimension
        //else abort
        if(shape[i] == DTF_UNLIMITED){
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: currently cannot support when unlimited dimension is not in the slowest changing dimension (dim 0). Aborting.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
    }
    if(dtype == MPI_DATATYPE_NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: datatype for var %d (file %s) is null. Aborting.", varid, fbuf->file_path);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    ret = def_var(fbuf, varid, ndims, dtype, shape);

    return ret;
}

_EXTERN_C_ void dtf_tstart()
{
    if(gl_stats.timer_start != 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: user timer was started at %.3f and not finished.",
                gl_stats.timer_start - gl_stats.walltime);

    gl_stats.timer_start = MPI_Wtime();
}
_EXTERN_C_ void dtf_tend()
{
    double tt = MPI_Wtime() - gl_stats.timer_start;
    gl_stats.timer_accum += tt;
    gl_stats.timer_start = 0;
    DTF_DBG(VERBOSE_DBG_LEVEL, "time_stat: user time %.4f", tt);
}

/************************************************  Fortran Interfaces  *********************************************************/

_EXTERN_C_ void dtf_tstart_()
{
    dtf_tstart();
}
_EXTERN_C_ void dtf_tend_()
{
    dtf_tend();
}

void dtf_init_(const char *filename, char *module_name, int* ierr)
{
    *ierr = dtf_init(filename, module_name);
}

void dtf_finalize_(int* ierr)
{
    *ierr = dtf_finalize();
}

void dtf_match_io_(const char *filename, int *ncid, int *intracomp_io_flag, int *ierr)
{
    *ierr = dtf_match_io(filename, *ncid, *intracomp_io_flag);
}

//void dtf_match_io_all_(int *rw_flag)
//{
//    dtf_match_io_all(*rw_flag);
//}

void dtf_print_(const char *str)
{
    dtf_print(str);
}

void dtf_match_multiple_(int *ncid)
{
    dtf_match_multiple(*ncid);
}

void dtf_complete_multiple_(const char *filename, int *ncid)
{
    dtf_complete_multiple(filename, *ncid);
}

void dtf_print_data_(int *varid, int *dtype, int *ndims, MPI_Offset* count, void* data)
{
    dtf_print_data(*varid, *dtype, *ndims, count, data);
}
