/*
 * Copyright (C) 2015, Advanced Institute for Computational Science, RIKEN
 * Author: Jianwei Liao(liaotoad@gmail.com)
 */


#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "dtf_util.h"
#include "dtf.h"
#include "dtf_file_buffer.h"
#include "dtf_req_match.h"
#include "rb_red_black_tree.h"

extern file_info_req_q_t *gl_finfo_req_q;
extern dtf_msg_t *gl_msg_q;

static int match_str(char *pattern, const char *filename)
{
	char *next_token, *cur_token;
    int pos, strpos;
    int ret = 0;
    
    next_token = strchr(pattern, '%');
	
	if(next_token == NULL){
		//not a name pattern, just check for inclusion
		if(strstr(pattern, filename)!=NULL || strstr(filename, pattern)!=NULL)
			ret = 1;
	} else {
		strpos = 0;
		cur_token = &pattern[strpos];
		
		while( next_token != NULL){
			pos = next_token-pattern;
			pattern[pos] = '\0';
			
			if(pos - strpos > 0){ //check for inclusion
				//printf("token %s\n", cur_token);
				if(strstr(filename, cur_token) == NULL){
					ret = 0;
					pattern[pos] = '%';
					//printf("->not found\n");
					break;
				} else {
					ret = 1;
					//printf("->found\n");
				}
			}
			pattern[pos] = '%';
			
			strpos = pos+1;
			cur_token = &pattern[strpos];
			next_token = strchr(cur_token, '%');
		}
		if( ret && (strpos != strlen(pattern))){
			//printf("token %s\n", cur_token);
			//check the last token
			if(strstr(filename, cur_token) == NULL){
				ret = 0;
				//printf("->not found\n");
			} else {
				//printf("->found\n");
			}
		}
	}
	return ret;
}

int match_ptrn(char* pattern, const char* filename, char** excl_fnames, int nexcls)
{
	int i;
    if(strlen(pattern)== 0 || strlen(filename)==0)
        return 0;
    
    /*First see if the filename matches against any of the exclude patterns*/
    for(i = 0; i < nexcls; i++)
		if(match_str(excl_fnames[i], filename)) return 0;
	
    return match_str(pattern, filename);
}

MPI_Offset to_1d_index(int ndims, const MPI_Offset *block_start, const MPI_Offset *block_count, const MPI_Offset *coord)
{
      int i, j;
      MPI_Offset idx = 0, mem=0;

      if(ndims == 0) //scalar
        return 0;
      else if(ndims == 1) //1d array
        return *coord - *block_start;

      for(i = 0; i < ndims; i++){
        mem = coord[i] - block_start[i];
        for(j = i+1; j < ndims; j++)
          mem *= block_count[j];
        idx += mem;
      }
    return idx;
}

void notify_file_ready(file_buffer_t *fbuf)
{

    int err;
    char *filename;
    dtf_msg_t *msg;
    if(fbuf->iomode != DTF_IO_MODE_FILE) return;
    if(gl_my_rank != fbuf->root_writer) return;
    
    DTF_DBG(VERBOSE_DBG_LEVEL, "Inside notify_file_ready");
	assert(fbuf->root_reader != -1);
            
	filename = dtf_malloc(MAX_FILE_NAME);
	assert(filename != NULL);
	strcpy(filename, fbuf->file_path);
	msg = new_dtf_msg(filename, MAX_FILE_NAME, FILE_READY_TAG);
	DTF_DBG(VERBOSE_DBG_LEVEL,   "Notify reader root rank %d that file %s is ready", fbuf->root_reader, fbuf->file_path);
	err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader, FILE_READY_TAG, gl_comps[fbuf->reader_id].comm, &(msg->req));
	CHECK_MPI(err);
	ENQUEUE_ITEM(msg, gl_msg_q);
	fbuf->fready_notify_flag = RDR_NOTIF_POSTED;       
    
}

void print_stats()
{
    char *s;
    int err, nranks;
    double dblsum = 0, walltime, avglibt;
    unsigned long data_sz;
    int sclltkf = 0;
    unsigned long unsgn;
    typedef struct{
        double dbl;
        int intg;
    }dblint_t;

    dblint_t dblint_in, dblint_out;
    
    MPI_Barrier(MPI_COMM_WORLD);
	
    walltime = MPI_Wtime() - gl_stats.walltime;
    MPI_Comm_size(gl_comps[gl_my_comp_id].comm, &nranks);

    /*DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: walltime %.3f", walltime);

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
*/
//    if(gl_stats.accum_match_time > 0){
//          DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: db match time: %.5f: %.5f",
//                 gl_stats.accum_db_match_time, (gl_stats.accum_db_match_time/gl_stats.timer_accum)*100);
//         DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: db manage time: %.5f: %.5f",
//                 gl_stats.accum_db_manage_time,  (gl_stats.accum_db_manage_time/gl_stats.timer_accum)*100);
//
//    }

    //DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: times matched: %u", gl_stats.ndb_match);
   // if(gl_stats.ndb_match > 0)
     //   DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: times dbmatched %d, progress call %lu", gl_stats.ndb_match, gl_stats.nprogress_call);

   // if(gl_stats.nioreqs > 0)
     //   DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: nreqs: %lu", gl_stats.nioreqs);

  //  if(gl_stats.nbl > 0)
    //    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: Out of %u blocks handled %u noncontig blocks", gl_stats.nbl, gl_stats.ngetputcall);
    if(gl_stats.accum_dbuff_sz > 0){
    //    DTF_DBG(VERBOSE_DBG_LEVEL, "DTF STAT: buffering time: %.5f: %.4f", gl_stats.accum_dbuff_time,(gl_stats.accum_dbuff_time/gl_stats.timer_accum)*100);
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: buffering size: %lu",  gl_stats.accum_dbuff_sz);
    }
    

    if(gl_stats.timer_accum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: library-measured I/O time: %.4f", gl_stats.timer_accum);
	if(gl_stats.user_timer_accum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: user-measured I/O time: %.4f", gl_stats.user_timer_accum);

    /*In scale-letkf the last ranks write mean files not treated by dtf, we don't need
      stats for that*/
    s = getenv("DTF_SCALE");
    if(s != NULL)
        sclltkf = atoi(s);
    else
        sclltkf = 0;

    if(sclltkf){
		s = getenv("SCALE_ENSEMBLE_SZ");
		if(s == 0) 
			DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: SCALE_ENSEMBLE_SZ is not set. Stats will be skewed");
		else {
			nranks = nranks - atoi(s);
			if(gl_my_rank == 0)
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF: for stats consider %d ranks out of %d", nranks, nranks + atoi(s));
		}
       
        if(gl_my_rank >= nranks){//Only for because the last set of processes work with mean files
            DTF_DBG(VERBOSE_DBG_LEVEL, "will zero stats");
            gl_stats.timer_accum = 0;
            walltime = 0;
            gl_stats.accum_match_time = 0;
            gl_stats.accum_hdr_time = 0;
            gl_stats.accum_comm_time = 0;
            gl_stats.accum_dbuff_sz = 0;
            gl_stats.accum_dbuff_time = 0;
            gl_stats.accum_progr_time = 0;
            gl_stats.accum_rw_var = 0;
            gl_stats.nioreqs = 0;
            gl_stats.user_timer_accum = 0;
        }
    }
	
	if(gl_my_rank == 0)
		rb_print_stats();

    /*AVERAGE STATS*/
    if(gl_stats.iodb_nioreqs > 0 && gl_my_rank == 0){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Stat: nioreqs in iodb %lu", gl_stats.iodb_nioreqs);
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Stat: parse time %.4f (%.3f%%)", gl_stats.parse_ioreq_time, gl_stats.parse_ioreq_time/gl_stats.timer_accum*100);
	}
    
    err = MPI_Reduce(&walltime, &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg walltime: %.4f", dblsum/nranks);

    err = MPI_Allreduce(&(gl_stats.timer_accum), &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    avglibt = dblsum/nranks;
    if(gl_my_rank==0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg library-measured I/O time: %.5f", avglibt);

    /*Standard deviation*/
    {
        double mydev2;
        double mean = avglibt;
        mydev2 = (gl_stats.timer_accum - mean)*(gl_stats.timer_accum - mean);

        err = MPI_Reduce(&mydev2, &dblsum, 1, MPI_DOUBLE, MPI_SUM,0, gl_comps[gl_my_comp_id].comm);
        CHECK_MPI(err);

        if(gl_my_rank == 0)
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: lib I/O standard deviation: %.7f", sqrt(dblsum/nranks));
    }
    
    err = MPI_Allreduce(&(gl_stats.user_timer_accum), &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    avglibt = dblsum/nranks;
    if(gl_my_rank==0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg user-measured I/O time: %.5f", avglibt);

    /*Standard deviation*/
    {
        double mydev2;
        double mean = avglibt;
        mydev2 = (gl_stats.user_timer_accum - mean)*(gl_stats.user_timer_accum - mean);

        err = MPI_Reduce(&mydev2, &dblsum, 1, MPI_DOUBLE, MPI_SUM,0, gl_comps[gl_my_comp_id].comm);
        CHECK_MPI(err);

        if(gl_my_rank == 0)
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: lib I/O standard user deviation: %.7f", sqrt(dblsum/nranks));
    }

    //dblint_in.dbl = walltime;
    //dblint_in.intg = gl_my_rank;
    //err = MPI_Reduce(&dblint_in, &dblint_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, gl_comps[gl_my_comp_id].comm);
    //CHECK_MPI(err);
    //if(gl_my_rank == 0)
        //DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: max walltime: %.4f: rank: %d", dblint_out.dbl, dblint_out.intg);

    dblint_in.dbl = gl_stats.timer_accum;
    dblint_in.intg = gl_my_rank;
    err = MPI_Reduce(&dblint_in, &dblint_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: max libtime: %.4f: rank: %d", dblint_out.dbl, dblint_out.intg);

    err = MPI_Reduce(&(gl_stats.accum_match_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg tot match time: %.4f (%.4f%%)", dblsum/nranks, (dblsum/nranks)/avglibt*100);

    err = MPI_Reduce(&(gl_stats.accum_hdr_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg hdr time: %.4f (%.4f%%)", dblsum/nranks, (dblsum/nranks)/avglibt*100);
    //~ if(gl_my_rank == 0 && gl_stats.accum_hdr_time > 0){
        //~ DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: hdr time: %.4f (%.4f%%), idle do match %.4f",
                //~ gl_stats.accum_do_matching_time, gl_stats.accum_do_matching_time/avglibt*100, gl_stats.idle_do_match_time);
    //~ }

    //DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: my do match time: %.4f (%.4f%%)", gl_stats.accum_do_matching_time, gl_stats.accum_do_matching_time/avglibt*100);

    //err = MPI_Reduce(&(gl_stats.idle_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    //CHECK_MPI(err);
    //if(gl_my_rank == 0 && dblsum > 0)
        //DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg idle time: %.4f (%.4f%%)", dblsum/nranks, (dblsum/nranks)/avglibt*100);


    //err = MPI_Reduce(&(gl_stats.accum_progr_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    //CHECK_MPI(err);
    //if(gl_my_rank == 0 && dblsum > 0)
        //DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg progr time: %.4f (%.4f%%)", dblsum/nranks,  (dblsum/nranks)/avglibt*100);
	if(gl_my_rank == 0)
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: my progr time: %.4f (%.3f%%)", gl_stats.accum_progr_time, gl_stats.accum_progr_time/gl_stats.timer_accum*100);
    err = MPI_Reduce(&(gl_stats.accum_rw_var), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: rw var time: %.4f (%.4f%%)", dblsum/nranks,  (dblsum/nranks)/avglibt*100);

    err = MPI_Reduce(&(gl_stats.accum_comm_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum>0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg comm time: %.5f (%.4f%%)", dblsum/nranks, (dblsum/nranks)/avglibt*100);

   ////err = MPI_Reduce(&(gl_stats.accum_dbuff_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
   //// CHECK_MPI(err);
  ////  if(gl_my_rank==0 && dblsum > 0)
     ////   DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg buffering time: %.5f (%.4f%%)", (double)(dblsum/nranks), (dblsum/nranks)/avglibt*100);

    if(gl_stats.ndata_msg_sent > 0 && gl_my_rank == 0){
        data_sz = (unsigned long)(gl_stats.data_msg_sz/gl_stats.ndata_msg_sent);
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: total sent %lu, avg data msg sz %lu (%d msgs)", gl_stats.data_msg_sz, data_sz, gl_stats.ndata_msg_sent);
    } 

    //~ err = MPI_Reduce(&data_sz, &lngsum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    //~ CHECK_MPI(err);
    //~ if(gl_my_rank==0 && lngsum > 0)
        //~ DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: Avg data msg sz acrosss ps: %lu", (unsigned long)(lngsum/nranks));
    err = MPI_Reduce(&(gl_stats.nioreqs), &unsgn, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg nioreqs: %u", (unsigned)(unsgn/nranks));

//    err = MPI_Reduce(&gl_stats.data_msg_sz, &lngsum, 1, MPI_UNSIGNED_LONG, MPI_MAXLOC, 0, gl_comps[gl_my_comp_id].comm);
//    CHECK_MPI(err);
  //  if(gl_my_rank==0 && gl_stats.data_msg_sz > 0)
    //    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: total data sent acrosss ps by rank 0: %lu", gl_stats.data_msg_sz);
    //~ intsum = 0;
    //~ err = MPI_Reduce(&(gl_stats.ndata_msg_sent), &intsum, 1, MPI_INT, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    //~ CHECK_MPI(err);
    //~ if(gl_my_rank==0 && intsum > 0)
        //~ DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: Avg num data msg: %d", (int)(intsum/nranks));

    //if(gl_stats.master_time > 0){
        //DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: master time %.4f (%.4f%% of libtime)",
        //gl_stats.master_time, gl_stats.master_time/gl_stats.timer_accum*100);
    //}

   // err = MPI_Reduce(&(gl_stats.accum_dbuff_sz), &lngsum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    //CHECK_MPI(err);
    //if(gl_my_rank==0 && dblsum > 0)
      //  DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg buffering size: %lu", (size_t)(lngsum/nranks));
}


//TODO check if didnt miss anything from finalize_files()
void close_file(file_buffer_t *fbuf)
{
    if(fbuf->writer_id == gl_my_comp_id){
		
		fbuf->is_ready = 1;

        if(fbuf->iomode == DTF_IO_MODE_FILE) {
			MPI_Barrier(fbuf->comm);
			//Check for any incoming messages
			progress_io_matching();
            
            if(gl_my_rank == fbuf->root_writer){
				if(fbuf->fready_notify_flag == RDR_NOT_NOTIFIED)
					notify_file_ready(fbuf);
			} 

        } else if(fbuf->iomode == DTF_IO_MODE_MEMORY){
			int finalize = 1;
            DTF_DBG(VERBOSE_DBG_LEVEL, "Cleaning up everything");
            delete_ioreqs(fbuf, finalize);
            if(fbuf->mst_info->iodb != NULL)
                finalize_iodb(fbuf);

            /*File is opened and closed multiple times in SCALE-LETKF
              but it's the same set of processes, hence, don't delete the master related data.
            */
            
        }
    } else if (fbuf->reader_id == gl_my_comp_id){
            assert(fbuf->rreq_cnt == 0);
            assert(fbuf->wreq_cnt == 0);

            if(fbuf->root_reader == gl_my_rank && fbuf->iomode == DTF_IO_MODE_MEMORY){

				while(!fbuf->done_match_confirm_flag){
					progress_io_matching();
				}

                DTF_DBG(VERBOSE_DBG_LEVEL, "Done");
            }
    }
    
    //delete_file_buffer(fbuf);
}


int inquire_root(const char *filename)
{
	int i;
	int root_writer = -1;
	
	FILE *rootf;
	char *glob = getenv("DTF_GLOBAL_PATH");
	assert(glob != NULL);
	char outfname[MAX_FILE_NAME*2];
	char fname[MAX_FILE_NAME];
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Inquire who is the root writer for file %s", filename);
	
	strcpy(fname, filename);
	for(i = 0; i <strlen(filename); i++)
		if(fname[i]=='/' || fname[i]=='\\')
			fname[i]='_';
	strcpy(outfname, glob);
	strcat(outfname, "/");
	strcat(outfname, fname);
	strcat(outfname, ".root");
	
	
	while( access( outfname, F_OK ) == -1 )
		{sleep(1);};
	
	rootf = fopen(outfname, "rb");
	while(1){
		if(!fread(&root_writer, sizeof(int), 1, rootf))
			sleep(1);
		else 
			break;
	}
	fclose(rootf);
	
//	if(remove(outfname))
	//	DTF_DBG(VERBOSE_ERROR_LEVEL,"DTF Warning: couldnt delete the root file %s", outfname);

	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Root writer is %d", root_writer);
	return root_writer;
}

void open_file(file_buffer_t *fbuf, MPI_Comm comm)
{
    DTF_DBG(VERBOSE_DBG_LEVEL,   "Enter dtf_open %s", fbuf->file_path);

    MPI_Status status;
    int rank; //, notif_open=1;

    int err;

    if(fbuf->reader_id == gl_my_comp_id){
        MPI_Comm_rank(comm, &rank);
        if(rank == 0)
            fbuf->root_reader = gl_my_rank;
        err = MPI_Bcast(&fbuf->root_reader, 1, MPI_INT, 0, comm);
        CHECK_MPI(err);

        if(fbuf->iomode == DTF_IO_MODE_FILE){
            if(fbuf->root_writer == -1){
				MPI_Request req;
                if(rank == 0){
					fbuf->root_writer = inquire_root(fbuf->file_path);
                    DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writer root that I am reader root of file %s", fbuf->file_path);
                    err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_writer, FILE_INFO_REQ_TAG, gl_comps[fbuf->writer_id].comm, &req);
                    CHECK_MPI(err);
                }
                DTF_DBG(VERBOSE_DBG_LEVEL, "Broadcast root to other readers");
                err = MPI_Bcast(&fbuf->root_writer, 1, MPI_INT, 0, comm);
                CHECK_MPI(err);
                if(rank == 0){
					err = MPI_Wait(&req, MPI_STATUS_IGNORE);
					CHECK_MPI(err);
				}
            }
            //NOTE: uncomment this for multi-cycle version where we open/close file only once
            //fbuf->is_ready = 0;
            double t_start = MPI_Wtime();
//            if(fbuf->is_ready){
//                notif_open = 0; //already notified before
//                DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: file is already ready");
//            } else

            
			 if(rank == 0){
				DTF_DBG(VERBOSE_DBG_LEVEL, "Waiting for file to become ready");
				/*root reader rank will wait until writer finishes writing the file.
				 then it will broadcast that the file is ready to everyone else*/

				while(!fbuf->is_ready)
					progress_io_matching();
			}
			
            err = MPI_Bcast(&fbuf->is_ready, 1, MPI_INT, 0, comm);
            DTF_DBG(VERBOSE_DBG_LEVEL, "PROFILE: Waiting to open file %.3f", MPI_Wtime()-t_start);
            CHECK_MPI(err);

        } else if(fbuf->iomode == DTF_IO_MODE_MEMORY){
            /*First, find out who is the root master*/
            int bufsz;
            void *buf;
            
            if(fbuf->root_writer == -1){            
				/*Zero rank will inquire the pnetcdf header/dtf vars info/masters info
				from writer's global zero rank and then broadcast this info to other
				readers that opened the file*/
				if(rank == 0){					
					int nranks;
					MPI_Comm_size(comm, &nranks);
                	fbuf->mst_info->nranks_opened = nranks;
					fbuf->root_writer = inquire_root(fbuf->file_path);
					DTF_DBG(VERBOSE_DBG_LEVEL, "Ask writer root for file info for %s", fbuf->file_path);
					err = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_writer, FILE_INFO_REQ_TAG, gl_comps[fbuf->writer_id].comm);
					CHECK_MPI(err);
					
					DTF_DBG(VERBOSE_DBG_LEVEL, "Starting to wait for file info for %s", fbuf->file_path);
					err = MPI_Probe(fbuf->root_writer, FILE_INFO_TAG, gl_comps[fbuf->writer_id].comm, &status);
					CHECK_MPI(err);
					MPI_Get_count(&status, MPI_BYTE, &bufsz);
					
					buf = dtf_malloc(bufsz);
					assert(buf != NULL);
					err = MPI_Recv(buf, bufsz, MPI_BYTE, fbuf->root_writer, FILE_INFO_TAG, gl_comps[fbuf->writer_id].comm, &status);
					CHECK_MPI(err);
				}
			   // MPI_Barrier(comm);
				DTF_DBG(VERBOSE_DBG_LEVEL, "Bcast file info to others");
				err = MPI_Bcast(&bufsz, 1, MPI_INT, 0, comm);
				CHECK_MPI(err);
				assert(bufsz > 0);

				if(rank != 0){
					buf = dtf_malloc(bufsz);
					assert(buf != NULL);
				}
				err = MPI_Bcast(buf, bufsz, MPI_BYTE, 0, comm);
				CHECK_MPI(err);

				unpack_file_info(bufsz, buf);
				dtf_free(buf, bufsz);
			}
        	//TODO this is not the same behavior for mem and file modes. is it ok? is it needed?

			fbuf->is_ready = 1;
			fbuf->done_match_confirm_flag = 1;
        }

    } else if(fbuf->writer_id == gl_my_comp_id){
		//TODO reinit iodb, allocate stuff for witems and ritems
		//TODO do not destroy fbuf when file closed, only destroy in finalize
        //do we need it?
		
      //  assert(0);
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: Writer component (%s) opened file %s instead of creating it", gl_comps[fbuf->writer_id].name, fbuf->file_path);
        /*reset all flags*/
    }

    DTF_DBG(VERBOSE_DBG_LEVEL,   "Exit dtf_open %s", fbuf->file_path);
}

int def_var(file_buffer_t *fbuf, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int el_sz;
    int i;
    dtf_var_t *var = new_var(varid, ndims, dtype, shape);
    add_var(fbuf, var);
    MPI_Type_size(dtype, &el_sz);

    DTF_DBG(VERBOSE_DBG_LEVEL, "varid %d, dim %d, el_sz %d. shape:", varid, ndims, el_sz);
    for(i = 0; i < ndims; i++)
        DTF_DBG(VERBOSE_DBG_LEVEL, "\t%lld", shape[i]);

    return 0;
}

/*Write pnetcdf header*/
void write_hdr(file_buffer_t *fbuf, MPI_Offset hdr_sz, void *header)
{
    DTF_DBG(VERBOSE_DBG_LEVEL, "Writing header (sz %d)", (int)hdr_sz);
    fbuf->hdr_sz = hdr_sz;
    fbuf->header = dtf_malloc(hdr_sz);
    assert(fbuf->header != NULL);
    memcpy(fbuf->header, header, (size_t)hdr_sz);
    return;
}

/*Read pnetcdf header*/
MPI_Offset read_hdr_chunk(file_buffer_t *fbuf, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
    if(offset+chunk_sz > fbuf->hdr_sz){
        DTF_DBG(VERBOSE_ALL_LEVEL, "Warning: trying to read %llu at offt %llu but hdr sz is %llu", chunk_sz, offset, fbuf->hdr_sz);
        chunk_sz = fbuf->hdr_sz - offset;
    }

    memcpy(chunk, (unsigned char*)fbuf->header+offset, (size_t)chunk_sz);
    return chunk_sz;
}

/*Find the biggest subblock of data that can fit into given buffer size sbufsz
  [INOUT] cur_count - the subblock that fits
  [INOUT] cur_nelems - how many elements managed to fit. Zero means didn't fit anything.
  [IN] the rest
*/
void find_fit_block(int ndims,
		    int cur_dim,
		    const MPI_Offset *start,
		    const MPI_Offset *count,
		    MPI_Offset *cur_start,
		    MPI_Offset *cur_count,
		    const size_t sbufsz,
		    const size_t el_sz,
		    MPI_Offset *cur_nelems,
		    MPI_Offset tot_nelems)
{
  int i;

  if(*cur_nelems == tot_nelems)
    return;
  if(sbufsz == 0)
    return;

  /*Compute if current subblock can fit*/
  MPI_Offset full_subbl = 1;
  for(i = 0; i < ndims; i++)
    if(i != cur_dim)
      full_subbl *= cur_count[i];

  DTF_DBG(VERBOSE_DBG_LEVEL, "full subblock %lld, cur dim %d, ndims %d", full_subbl, cur_dim, ndims);

  for(i = count[cur_dim] - (cur_start[cur_dim] - start[cur_dim]); i > 0; i--){
    //DTF_DBG(VERBOSE_DBG_LEVEL, "left %lld, right %lld", full_subbl * i * el_sz, (MPI_Offset)sbufsz);
    if(full_subbl * i * el_sz <= (MPI_Offset)sbufsz)
      break;
  }

  cur_count[cur_dim] = i;
  assert(cur_count[cur_dim] > 0);

  *cur_nelems = 1;
  for(i = 0; i < ndims; i++){
    *cur_nelems *= cur_count[i];
  }

  if(cur_count[cur_dim] == count[cur_dim]){
    if(cur_dim == 0){
      assert(*cur_nelems == tot_nelems);
      return;
    } else {
      DTF_DBG(VERBOSE_DBG_LEVEL, "Go higher. Cur nelems %lld. Cur count:", *cur_nelems);
      for(i=0; i < ndims; i++)
        DTF_DBG(VERBOSE_DBG_LEVEL, "  %lld", cur_count[i]);
      /*Go higher one dimension*/
      find_fit_block(ndims, cur_dim - 1, start, count, cur_start, cur_count, sbufsz, el_sz, cur_nelems, tot_nelems);
    }
  }
}

void shift_coord(int ndims, const MPI_Offset *bl_start,
                 const MPI_Offset *bl_count, MPI_Offset *subbl_start,
                 MPI_Offset *subbl_count, MPI_Offset fit_nelems)
{
    int i;

    /*Shift the start position*/
    if(fit_nelems == 1){ //special case
      subbl_start[ndims-1]++;
    } else {
        for(i = 0; i < ndims; i++)
            if(subbl_count[i] > 1)
                  subbl_start[i] += subbl_count[i];
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "New start before adjustment:");
    for(i = 0; i < ndims; i++)
        DTF_DBG(VERBOSE_DBG_LEVEL, "\t %lld", subbl_start[i]);

    for(i = ndims - 1; i > 0; i--)
        if(subbl_start[i] == bl_start[i] + bl_count[i]){
            subbl_start[i] = bl_start[i];
            if( (subbl_start[i-1] != bl_start[i-1] + bl_count[i-1]) && (subbl_count[i-1] == 1)){
                subbl_start[i-1]++;
            }
        } else
            break;

    DTF_DBG(VERBOSE_DBG_LEVEL, "New start after adjustment:");
    for(i = 0; i < ndims; i++)
        DTF_DBG(VERBOSE_DBG_LEVEL, "\t %lld", subbl_start[i]);


//    DTF_DBG(VERBOSE_DBG_LEVEL, "Copied subblock. Shift start:");
//    for(i = 0; i < var->ndims; i++)
//        DTF_DBG(VERBOSE_DBG_LEVEL, "   %lld\t -->\t %lld", bl_start[i], subbl_start[i]);
}

double compute_checksum(void *arr, int ndims, const MPI_Offset *shape, MPI_Datatype dtype)
{
    double sum = 0;
    unsigned nelems;
    int i;

    if(dtype != MPI_DOUBLE && dtype != MPI_FLOAT){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: checksum supported only for double or float data");
        return 0;
    }

    if(ndims == 0){
        if(dtype == MPI_DOUBLE)
            sum = *(double*)arr;
        else
            sum = *(float*)arr;
        return sum;
    }

    nelems = shape[0];
    for(i = 1; i < ndims; i++)
        nelems *= shape[i];

    for(i = 0; i < nelems; i++)
        if(dtype == MPI_DOUBLE)
            sum += ((double*)arr)[i];
        else
            sum += ((float*)arr)[i];
    return sum;
}

/*only support conversion double<->float*/
void get_put_data(dtf_var_t *var,
                  MPI_Datatype dtype,
                  unsigned char *block_data,
                  const MPI_Offset *block_start,
                  const MPI_Offset *block_count,
                  const MPI_Offset subbl_start[],
                  const MPI_Offset subbl_count[],
                  unsigned char *subbl_data,
                  int get_put_flag,
                  int convert_flag)
{
    int i;
    MPI_Offset *cur_coord;
    int nelems;

    if(var->ndims == 0){
        nelems = 1;
        int el_sz;
        MPI_Type_size(var->dtype, &el_sz);

        if(get_put_flag == DTF_READ){
            if(convert_flag)
                convertcpy(dtype, var->dtype, (void*)block_data, (void*)subbl_data, nelems);
            else
                memcpy(subbl_data, block_data, el_sz);

        } else { /*DTF_WRITE*/
           /*copy data subblock -> block*/
            if(convert_flag)
                convertcpy(var->dtype, dtype, (void*)subbl_data,(void*)block_data, nelems);
            else
                memcpy(block_data, subbl_data, el_sz);
        }
        gl_stats.ngetputcall++;
        return;
    }

    cur_coord = dtf_malloc(var->ndims*sizeof(MPI_Offset));
    for(i = 0; i < var->ndims; i++){
        cur_coord[i] = subbl_start[i];
    }
    DTF_DBG(VERBOSE_DBG_LEVEL, "Call getput data");
    /*read from user buffer to send buffer*/
    recur_get_put_data(var, dtype, block_data, block_start,
                       block_count, subbl_start, subbl_count, 0,
                       cur_coord, subbl_data, get_put_flag, convert_flag);
    dtf_free(cur_coord, var->ndims*sizeof(MPI_Offset));
    gl_stats.ngetputcall++;
}

/*only support conversion double<->float*/
void recur_get_put_data(dtf_var_t *var,
                          MPI_Datatype dtype,
                          unsigned char *block_data,
                          const MPI_Offset *block_start,
                          const MPI_Offset *block_count,
                          const MPI_Offset subbl_start[],
                          const MPI_Offset subbl_count[],
                          int dim,
                          MPI_Offset coord[],
                          unsigned char *subbl_data,
                          int get_put_flag,
                          int convert_flag)
{
    int i;
    if(dim == var->ndims - 1){
        int bl_type_sz, subbl_type_sz;
        MPI_Type_size(dtype, &bl_type_sz);
        MPI_Type_size(var->dtype, &subbl_type_sz);
        MPI_Offset block_offt = to_1d_index(var->ndims, block_start, block_count, coord)*bl_type_sz;
        MPI_Offset subbl_offt = to_1d_index(var->ndims, subbl_start, subbl_count, coord)*subbl_type_sz;
        MPI_Offset nelems = subbl_count[var->ndims-1];
        //MPI_Offset data_sz = subbl_count[var->ndims-1]*subbl_type_sz;    //data_sz

        if(get_put_flag == DTF_READ){
            if(convert_flag){
                convertcpy(dtype, var->dtype, (void*)(block_data+block_offt), (void*)(subbl_data+subbl_offt), (int)nelems);
//                for(i = 0; i < subbl_count[var->ndims-1]; i++)
//                    ((float*)(subbl_data+subbl_offt))[i] = (float)((double*)(block_data+block_offt))[i];
            } else
                /*copy data block -> subblock*/
                memcpy(subbl_data+subbl_offt, block_data+block_offt, nelems*subbl_type_sz);
        } else { /*DTF_WRITE*/
           /*copy data subblock -> block*/
            if(convert_flag)
                convertcpy(var->dtype, dtype, (void*)(subbl_data+subbl_offt),(void*)(block_data+block_offt), nelems);
            else
                memcpy(block_data+block_offt, subbl_data+subbl_offt, nelems*subbl_type_sz);
        }
        return;
    }

    for(i = 0; i < subbl_count[dim]; i++){
        coord[dim] = subbl_start[dim] + i;
        recur_get_put_data(var, dtype, block_data, block_start, block_count,
                           subbl_start, subbl_count, dim+1, coord, subbl_data,
                           get_put_flag, convert_flag);
    }
}

int mpitype2int(MPI_Datatype dtype)
{

    if(dtype == MPI_SIGNED_CHAR)     return 1;
    if(dtype == MPI_CHAR)            return 2;
    if(dtype == MPI_SHORT)           return 3;
    if(dtype == MPI_INT)             return 4;
    if(dtype == MPI_FLOAT)           return 5;
    if(dtype == MPI_DOUBLE)          return 6;
    if(dtype == MPI_UNSIGNED_CHAR)   return 7;
    if(dtype == MPI_UNSIGNED_SHORT)  return 8;
    if(dtype == MPI_UNSIGNED)        return 9;
    if(dtype == MPI_LONG_LONG_INT)   return 10;
    if(dtype == MPI_UNSIGNED_LONG_LONG) return 11;

    DTF_DBG(VERBOSE_ERROR_LEVEL, "Unknown mpi type");
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    return 0;
}


MPI_Datatype int2mpitype(int num)
{
    switch(num){
        case 1 :      return MPI_SIGNED_CHAR;
        case 2 :      return MPI_CHAR;
        case 3 :      return MPI_SHORT;
        case 4 :      return MPI_INT;
        case 5 :      return MPI_FLOAT;
        case 6 :      return MPI_DOUBLE;
        case 7 :      return MPI_UNSIGNED_CHAR;
        case 8 :      return MPI_UNSIGNED_SHORT;
        case 9 :      return MPI_UNSIGNED;
        case 10 :     return MPI_LONG_LONG_INT;
        case 11 :     return MPI_UNSIGNED_LONG_LONG;
        default:
            DTF_DBG(VERBOSE_ERROR_LEVEL, "Unknown mpi type");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }
    return MPI_DATATYPE_NULL;
}

void* dtf_malloc(size_t size)
{
    gl_stats.malloc_size += size;
    return malloc(size);
}

void dtf_free(void *ptr, size_t size)
{
    if(size > gl_stats.malloc_size)
        DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: mem stat negative (left %lu), will free %lu", gl_stats.malloc_size, size);
    gl_stats.malloc_size -= size;
    free(ptr);
    return;
}

void convertcpy(MPI_Datatype type1, MPI_Datatype type2, void* srcbuf, void* dstbuf, int nelems)
{
    int i;
    if(type1 == MPI_FLOAT){
        assert(type2 == MPI_DOUBLE);
        for(i = 0; i < nelems; i++)
            ((double*)dstbuf)[i] = (double)(((float*)srcbuf)[i]);
    } else if(type1 == MPI_DOUBLE){
        assert(type2 == MPI_FLOAT);
        for(i = 0; i < nelems; i++)
            ((float*)dstbuf)[i] = (float)(((double*)srcbuf)[i]);
    }
}

dtf_msg_t *new_dtf_msg(void *buf, size_t bufsz, int tag)
{
    dtf_msg_t *msg = dtf_malloc(sizeof(struct dtf_msg));
    assert(msg != NULL);
    msg->req = MPI_REQUEST_NULL;
    if(bufsz > 0){
        msg->buf = buf;
        //dtf_malloc(bufsz);
        //assert(msg->buf != NULL);
        //memcpy(msg->buf, buf, bufsz);
    } else
        msg->buf = NULL;
    msg->bufsz = bufsz;
    msg->tag = tag;
    msg->next = NULL;
    msg->prev = NULL;

    return msg;
}

void delete_dtf_msg(dtf_msg_t *msg)
{
    if(msg->bufsz > 0)
        dtf_free(msg->buf, msg->bufsz);
    dtf_free(msg, sizeof(msg));
}
