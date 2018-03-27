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
extern dtf_msg_t *gl_out_msg_q;

/*only support conversion double<->float*/
static void recur_get_put_data(dtf_var_t *var,
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
    //~ DTF_DBG(VERBOSE_DBG_LEVEL, "recur getput block:");
    //~ for(i = 0; i < var->ndims; i++)
		//~ DTF_DBG(VERBOSE_DBG_LEVEL, "bl: %lld -> %lld, subbl: %lld -> %lld", block_start[i], block_count[i], subbl_start[i], subbl_count[i]);

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


static void delete_dtf_msg(dtf_msg_t *msg)
{
    if(msg->bufsz > 0)
        dtf_free(msg->buf, msg->bufsz);
    dtf_free(msg, sizeof(msg));
}

void notify_file_ready(file_buffer_t *fbuf)
{
    int err;
    char *filename;
    dtf_msg_t *msg;
    if(fbuf->iomode != DTF_IO_MODE_FILE) return;
    assert(gl_my_rank == fbuf->root_writer);

    DTF_DBG(VERBOSE_DBG_LEVEL, "Inside notify_file_ready");
	assert(fbuf->root_reader != -1);

	filename = dtf_malloc(MAX_FILE_NAME);
	assert(filename != NULL);
	strcpy(filename, fbuf->file_path);
	msg = new_dtf_msg(filename, MAX_FILE_NAME, DTF_UNDEFINED, FILE_READY_TAG);
	DTF_DBG(VERBOSE_DBG_LEVEL,   "Notify reader root rank %d that file %s is ready", fbuf->root_reader, fbuf->file_path);
	err = MPI_Isend(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, fbuf->root_reader, FILE_READY_TAG, gl_comps[fbuf->reader_id].comm, &(msg->req));
	CHECK_MPI(err);
	ENQUEUE_ITEM(msg, gl_out_msg_q);
	fbuf->fready_notify_flag = RDR_NOTIF_POSTED;
	progress_send_queue();
}

/*Standard deviation*/
static double stand_devi(double myval, double sum, int nranks)
{
	int err;
	double tmpsum=0;
	double mean = sum/(double)nranks;
	double tmp = (myval - mean)*(myval - mean);

	if(gl_my_rank >=nranks) tmp = 0;
	err = MPI_Reduce(&tmp, &tmpsum, 1, MPI_DOUBLE, MPI_SUM,0, gl_comps[gl_my_comp_id].comm);
	CHECK_MPI(err);

    return sqrt(tmpsum/(double)nranks);

}

void print_stats()
{
    char *s;
    int err, nranks, enssz;
    double dblsum = 0, walltime, avglibt, dev;
    unsigned long data_sz;
    unsigned long unsgn;
    typedef struct{
        double dbl;
        int intg;
    }dblint_t;

    dblint_t dblint_in, dblint_out;

    walltime = MPI_Wtime() - gl_stats.walltime;
    MPI_Comm_size(gl_comps[gl_my_comp_id].comm, &nranks);

    if(gl_stats.accum_dbuff_sz > 0){
    //    DTF_DBG(VERBOSE_DBG_LEVEL, "DTF STAT: buffering time: %.5f: %.4f", gl_stats.accum_dbuff_time,(gl_stats.accum_dbuff_time/gl_stats.timer_accum)*100);
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: buffering size: %lu",  gl_stats.accum_dbuff_sz);
    }

    if(gl_stats.timer_accum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: pnetcdf I/O time: %.4f", gl_stats.timer_accum);
	if(gl_stats.user_timer_accum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: user-measured I/O time: %.4f", gl_stats.user_timer_accum);

	if(gl_stats.transfer_time > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: transfer time: %.4f", gl_stats.transfer_time);
        
    if(gl_stats.dtf_time > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: total dtf time: %.4f", gl_stats.dtf_time);

	err = MPI_Reduce(&walltime, &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg walltime: %.4f", dblsum/nranks);

	if(gl_my_rank == 0)
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT: my progr time: %.4f", gl_stats.t_progr_comm);

    /*In scale-letkf the last ranks write mean files not treated by dtf, we don't need
      stats for that*/
    if(gl_scale){
		
		err = MPI_Allreduce(&(gl_stats.st_fin), &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
		CHECK_MPI(err);
		dev = stand_devi(gl_stats.st_fin, dblsum, nranks);
		if(gl_my_rank == 0 && dblsum > 0)
			DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stamp avg st_fin, devi: %.4f: %.4f", dblsum/nranks, dev);


		s = getenv("SCALE_ENSEMBLE_SZ");
		if(s == 0)
			DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: SCALE_ENSEMBLE_SZ is not set. Stats will be skewed");
		else {
			enssz = atoi(s);
			nranks = nranks - atoi(s);
			if(gl_my_rank == 0)
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF: for stats consider %d ranks out of %d", nranks, nranks + atoi(s));
		}

        if(gl_my_rank >= nranks){//Only for because the last set of processes work with mean files
            DTF_DBG(VERBOSE_DBG_LEVEL, "will zero stats");
            gl_stats.timer_accum = 0;
            walltime = 0;
            gl_stats.transfer_time = 0;
            gl_stats.t_hdr = 0;
            gl_stats.t_comm = 0;
            gl_stats.accum_dbuff_sz = 0;
            gl_stats.accum_dbuff_time = 0;
            gl_stats.t_progr_comm = 0;
            gl_stats.t_rw_var = 0;
            gl_stats.nioreqs = 0;
            gl_stats.user_timer_accum = 0;
            gl_stats.st_fin=0;
            gl_stats.st_mtch_hist=0;
            gl_stats.st_mtch_rest=0;
            gl_stats.end_mtch_hist=0;
            gl_stats.end_mtch_rest=0;
            gl_stats.t_open_hist = 0;
            gl_stats.t_open_rest = 0;
            gl_stats.dtf_time = 0;
            gl_stats.transfer_time = 0;
        }

		dev = 0;

    	err = MPI_Allreduce(&gl_stats.st_mtch_hist, &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
		CHECK_MPI(err);
		dev = stand_devi(gl_stats.st_mtch_hist, dblsum, nranks);
		if(gl_my_rank == 0 && dblsum > 0)
			DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stamp avg st_mtch_hist, dev: %.4f: %.4f", dblsum/nranks, dev);
			
		err = MPI_Allreduce(&gl_stats.t_open_hist, &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
		CHECK_MPI(err);
		dev = stand_devi(gl_stats.t_open_hist, dblsum, nranks);
		if(gl_my_rank == 0 && dblsum > 0)
			DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stamp avg t_open_hist, dev: %.4f: %.4f", dblsum/nranks, dev);

		err = MPI_Allreduce(&gl_stats.end_mtch_hist, &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
		CHECK_MPI(err);
		dev = stand_devi(gl_stats.end_mtch_hist, dblsum, nranks);
		if(gl_my_rank == 0 && dblsum > 0)
			DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stamp avg end_mtch_hist, dev: %.4f: %.4f", dblsum/nranks, dev);

		err = MPI_Allreduce(&gl_stats.t_open_rest, &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
		CHECK_MPI(err);
		dev = stand_devi(gl_stats.t_open_rest, dblsum, nranks);
		if(gl_my_rank == 0 && dblsum > 0)
			DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stamp avg t_open_rest, dev: %.4f: %.4f", dblsum/nranks, dev);

		err = MPI_Allreduce(&gl_stats.st_mtch_rest, &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
		CHECK_MPI(err);
		dev = stand_devi(gl_stats.st_mtch_rest, dblsum, nranks);
		if(gl_my_rank == 0 && dblsum > 0)
			DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stamp avg st_mtch_rest, dev: %.4f: %.4f", dblsum/nranks, dev);

		err = MPI_Allreduce(&gl_stats.end_mtch_rest, &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
		CHECK_MPI(err);
		dev = stand_devi(gl_stats.end_mtch_rest, dblsum, nranks);
		if(gl_my_rank == 0 && dblsum > 0)
			DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stamp avg end_mtch_rest, dev: %.4f: %.4f", dblsum/nranks, dev);

		//~ err = MPI_Allreduce(&(gl_stats.st_fin), &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
		//~ CHECK_MPI(err);
		//~ dev = stand_devi(gl_stats.st_fin, dblsum, nranks);
		//~ if(gl_my_rank == 0 && dblsum > 0)
			//~ DTF_DBG(VERBOSE_ERROR_LEVEL, "time_stamp avg st_fin except last ens: %.4f: %.4f", dblsum/nranks, dev);
		{
			int i;
			MPI_Request req = MPI_REQUEST_NULL;
			MPI_Status status;
			double tmp;
			
			if(gl_my_rank % enssz == 0){
				err = MPI_Isend(&gl_stats.t_mtch_hist, 1, MPI_DOUBLE, 0, 0, gl_comps[gl_my_comp_id].comm, &req);
				CHECK_MPI(err);
			}
			if(gl_my_rank == 0){
				for(i = 0; i < (int)(nranks/enssz)+1; i++){
					err = MPI_Recv(&tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, gl_comps[gl_my_comp_id].comm, &status);
					CHECK_MPI(err);
					DTF_DBG(VERBOSE_ERROR_LEVEL, "Src %d: match hist %.4f", status.MPI_SOURCE, tmp);
				}
			}
			
			MPI_Wait(&req, MPI_STATUS_IGNORE);
			DTF_DBG(VERBOSE_ERROR_LEVEL, "done");
			if(gl_my_rank % enssz == 0){
				err = MPI_Isend(&gl_stats.t_mtch_rest, 1, MPI_DOUBLE, 0, 1, gl_comps[gl_my_comp_id].comm, &req);
				CHECK_MPI(err);
			}
			if(gl_my_rank == 0){
				for(i = 0; i < (int)(nranks/enssz)+1; i++){
					err = MPI_Recv(&tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, gl_comps[gl_my_comp_id].comm, &status);
					CHECK_MPI(err);
					DTF_DBG(VERBOSE_ERROR_LEVEL, "Src %d: match rest %.4f", status.MPI_SOURCE, tmp);
				}
			}
			MPI_Wait(&req, MPI_STATUS_IGNORE);
		}
    
    }

	if(gl_my_rank == 0)
		rb_print_stats();

    /*AVERAGE STATS*/
    if(gl_stats.iodb_nioreqs > 0 && gl_my_rank == 0){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Stat: nioreqs in iodb %lu", gl_stats.iodb_nioreqs);
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Stat: parse time %.4f (%.3f%%)", gl_stats.parse_ioreq_time, gl_stats.parse_ioreq_time/gl_stats.timer_accum*100);
	}

    err = MPI_Allreduce(&(gl_stats.timer_accum), &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    avglibt = dblsum/nranks;
    dev = stand_devi(gl_stats.timer_accum, dblsum, nranks);

    if(gl_my_rank==0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg pnetcdf time: %.5f : %.5f", avglibt, dev);

    err = MPI_Allreduce(&(gl_stats.user_timer_accum), &dblsum, 1, MPI_DOUBLE, MPI_SUM, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    avglibt = dblsum/nranks;
    dev = stand_devi(gl_stats.user_timer_accum, dblsum, nranks);
    if(gl_my_rank==0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg user timer time: %.5f : %.4f", avglibt, dev);


    err = MPI_Reduce(&(gl_stats.dtf_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg dtf time: %.4f: (%.4f%%)", dblsum/nranks, (dblsum/nranks)/avglibt*100);


    err = MPI_Reduce(&(gl_stats.transfer_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg transfer time: %.4f: (%.4f%%)", dblsum/nranks, (dblsum/nranks)/avglibt*100);

	dblint_in.dbl = gl_stats.dtf_time;
    dblint_in.intg = gl_my_rank;
    err = MPI_Reduce(&dblint_in, &dblint_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: max dtf time: %.4f: rank: %d", dblint_out.dbl, dblint_out.intg);


    //dblint_in.dbl = walltime;
    //dblint_in.intg = gl_my_rank;
    //err = MPI_Reduce(&dblint_in, &dblint_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, gl_comps[gl_my_comp_id].comm);
    //CHECK_MPI(err);
    //if(gl_my_rank == 0)
        //DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: max walltime: %.4f: rank: %d", dblint_out.dbl, dblint_out.intg);

    //~ dblint_in.dbl = gl_stats.timer_accum;
    //~ dblint_in.intg = gl_my_rank;
    //~ err = MPI_Reduce(&dblint_in, &dblint_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, gl_comps[gl_my_comp_id].comm);
    //~ CHECK_MPI(err);
    //~ if(gl_my_rank == 0)
        //~ DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: max libtime: %.4f: rank: %d", dblint_out.dbl, dblint_out.intg);

    {
		char *s=getenv("MAX_WORKGROUP_SIZE");
		assert(s!=NULL);
		int nmst = nranks/atoi(s);
		if (nranks %atoi(s) != 0)
			nmst++;

		err = MPI_Reduce(&(gl_stats.t_do_match), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
		CHECK_MPI(err);
		if(gl_my_rank == 0 && dblsum > 0){
			DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg do match time: %.4f", dblsum/nmst);
		}

		err = MPI_Reduce(&(gl_stats.idle_do_match_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
		CHECK_MPI(err);
		if(gl_my_rank == 0 && dblsum > 0){
			DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg idle do match time: %.4f", dblsum/nmst);
		}
	}




    //DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: my do match time: %.4f (%.4f%%)", gl_stats.t_do_match, gl_stats.t_do_match/avglibt*100);

    //err = MPI_Reduce(&(gl_stats.idle_time), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    //CHECK_MPI(err);
    //if(gl_my_rank == 0 && dblsum > 0)
        //DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg idle time: %.4f (%.4f%%)", dblsum/nranks, (dblsum/nranks)/avglibt*100);


	err = MPI_Reduce(&(gl_stats.t_progr_comm), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg progr time: %.4f", dblsum/nranks);
	
    err = MPI_Reduce(&(gl_stats.t_rw_var), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && dblsum > 0)
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: rw var time: %.4f (%.4f%%)", dblsum/nranks,  (dblsum/nranks)/avglibt*100);

    err = MPI_Reduce(&(gl_stats.t_comm), &dblsum, 1, MPI_DOUBLE, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
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

    err = MPI_Reduce(&(gl_stats.iodb_nioreqs), &unsgn, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, gl_comps[gl_my_comp_id].comm);
    CHECK_MPI(err);
    if(gl_my_rank == 0 && gl_stats.iodb_nioreqs> 0){
		char *s=getenv("MAX_WORKGROUP_SIZE");
		assert(s!=NULL);
		int nmst = nranks/atoi(s);
		if (nranks %atoi(s) != 0)
			nmst++;
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF STAT AVG: avg iodb nioreqs: %u", (unsigned)(unsgn/nmst));
	}

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

void send_mst_info(file_buffer_t *fbuf, int tgt_root, int tgt_comp)
{
	int bufsz, flag;
    void *buf;
    unsigned char *chbuf;
    int err;
	size_t offt = 0;
	//pack my master info and send it to the other component
	bufsz = MAX_FILE_NAME+fbuf->my_mst_info->nmasters*sizeof(int)+sizeof(int);
	buf = dtf_malloc(bufsz);
	assert(buf != NULL);
	chbuf = (unsigned char*)buf;
	
	memcpy(chbuf, fbuf->file_path, MAX_FILE_NAME);
	offt += MAX_FILE_NAME;
	/*number of masters*/
	memcpy(chbuf+offt, &(fbuf->my_mst_info->nmasters), sizeof(int));
	offt += sizeof(int);
	DTF_DBG(VERBOSE_DBG_LEVEL, "pack %d masters", fbuf->my_mst_info->nmasters);
	/*list of masters*/
	memcpy(chbuf+offt, fbuf->my_mst_info->masters, fbuf->my_mst_info->nmasters*sizeof(int));
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Notify writer root that I am reader root of file %s", fbuf->file_path);
	dtf_msg_t *msg = new_dtf_msg(buf, bufsz, DTF_UNDEFINED, FILE_INFO_REQ_TAG);
	err = MPI_Isend(buf, bufsz, MPI_BYTE, tgt_root, FILE_INFO_REQ_TAG, gl_comps[tgt_comp].comm, &(msg->req));
	CHECK_MPI(err);
	err = MPI_Test(&(msg->req), &flag, MPI_STATUS_IGNORE);
	CHECK_MPI(err);
	ENQUEUE_ITEM(msg, gl_out_msg_q);
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

dtf_msg_t *new_dtf_msg(void *buf, size_t bufsz, int src, int tag)
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
    msg->src = src; 
    return msg;
}

void progress_send_queue()
{
    dtf_msg_t *msg;
    int flag, err;
    MPI_Status status;
    double t_st;

    if(gl_out_msg_q == NULL)
       return;

    msg = gl_out_msg_q;
    while(msg != NULL){
        t_st = MPI_Wtime();
        err = MPI_Test(&(msg->req), &flag, &status);
        CHECK_MPI(err);
        if(flag){
            gl_stats.t_comm += MPI_Wtime() - t_st;
            dtf_msg_t *tmp = msg->next;
            DEQUEUE_ITEM(msg, gl_out_msg_q);
            if(msg->tag == FILE_READY_TAG){
                file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, (char*)msg->buf, -1);
                if(fbuf == NULL)
                    DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Error: cant find %s", (char*)msg->buf);
                assert(fbuf != NULL);
                assert(fbuf->fready_notify_flag == RDR_NOTIF_POSTED);
                fbuf->fready_notify_flag = DTF_UNDEFINED; //reset flag
                DTF_DBG(VERBOSE_DBG_LEVEL, "Completed sending file ready notif for %s", (char*)msg->buf);
            }
            delete_dtf_msg(msg);
            msg = tmp;
        } else{
            gl_stats.idle_time += MPI_Wtime() - t_st;
            msg = msg->next;
        }
    }
}

void progress_recv_queue()
{
	int comp;
	dtf_msg_t *msg, *tmp;
	
	for(comp = 0; comp < gl_ncomp; comp++){
		if(gl_comps[comp].in_msg_q == NULL)
			continue;
		msg = gl_comps[comp].in_msg_q;
		while(msg != NULL){
			if(parce_msg(comp, msg->src, msg->tag, msg->buf, msg->bufsz, 1)){
				tmp = msg->next;
				msg->bufsz = 0;
				DEQUEUE_ITEM(msg, gl_comps[comp].in_msg_q);
				delete_dtf_msg(msg);
				msg = tmp;
			} else
				msg = msg->next;
		}
	}
}


int boundary_check(file_buffer_t *fbuf, int varid, const MPI_Offset *start, const MPI_Offset *count )
{
    dtf_var_t *var = fbuf->vars[varid];

    if(var->ndims > 0){
        int i;

      /*  if(frt_indexing){
            for(i = 0; i < var->ndims; i++)
                if(var->shape[i] == DTF_UNLIMITED) //no boundaries for unlimited dimension
                    continue;
                else if(start[i] + count[i] > var->shape[var->ndims-i-1]){
                    DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: var %d, index %llu is out of bounds (shape is %llu)", varid, start[i]+count[i], var->shape[var->ndims-i-1]);
                    return 1;
                }
        } else { */
            for(i = 0; i < var->ndims; i++)
                if(var->shape[i] == DTF_UNLIMITED) //no boundaries for unlimited dimension
                    continue;
                else if(start[i] + count[i] > var->shape[i]){
                            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: var %d, index %llu is out of bounds (shape is %llu)", varid, start[i]+count[i], var->shape[i]);
                            return 1;
                }
       /* } */
    }
//    else {
//        if( (start != NULL) || (count != NULL)){
//            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: var %d is a scalar variable but trying to read an array", varid);
//            return 1;
//        }
//    }
    return 0;
}

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
