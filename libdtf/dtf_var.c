#include "dtf_util.h"
#include "dtf_var.h"
#include "dtf_req_match.h"
#include "dtf.h"

int def_var(struct file_buffer *fbuf, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
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


dtf_var_t* new_var(int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int i;

    dtf_var_t *var = (dtf_var_t*)dtf_malloc(sizeof(dtf_var_t));
    assert(var!=NULL);

    /*Initialize whatever we can initialize at this stage*/
    var->id = varid;
    if(ndims > 0){ 
        var->shape = (MPI_Offset*)dtf_malloc(ndims*sizeof(MPI_Offset));
        for(i=0; i<ndims;i++)
            var->shape[i] = shape[i];
    } else {
        var->shape = NULL;
	}
    var->ndims = ndims;
    var->dtype = dtype;
    var->checksum = 0;
    var->ioreqs = NULL;
    return var;
}

void delete_var(struct file_buffer *fbuf, dtf_var_t* var)
{
	assert(var->ioreqs == NULL);
    dtf_free(var->shape, var->ndims*sizeof(MPI_Offset));
    dtf_free(var, sizeof(dtf_var_t));
    fbuf->nvars--;
}

void add_var(struct file_buffer *fbuf, dtf_var_t *var)
{
    assert(var->id == fbuf->nvars);
    dtf_var_t **tmp = (dtf_var_t**)realloc(fbuf->vars, (fbuf->nvars+1)*sizeof(dtf_var_t*));
    assert(tmp != NULL);
    fbuf->vars = tmp;
    fbuf->nvars++;
    DTF_DBG(VERBOSE_DBG_LEVEL, "var id %d, cnt %d", var->id, fbuf->nvars);

    fbuf->vars[fbuf->nvars-1] = var;
    gl_stats.malloc_size += sizeof(dtf_var_t*);
}

MPI_Offset read_write_var(struct file_buffer *fbuf,
                               int varid,
                               const MPI_Offset *start,
                               const MPI_Offset *count,
                               const MPI_Offset *stride,
                               const MPI_Offset *imap,
                               MPI_Datatype dtype,
                               void *buf,
                               int rw_flag)
{
    MPI_Offset ret;
    io_req_t *req;
    int i;
    int def_el_sz, req_el_sz;
    MPI_Offset nelems;
    double t_start = MPI_Wtime();

    if(rw_flag == DTF_READ){
        if(fbuf->reader_id==gl_my_comp_id){
          assert(fbuf->is_ready);
        } else{
            DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: writer process tries to read file %s (var %d)", fbuf->file_path, varid);
            //MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
    }
    if(imap != NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: writing mapped vars is not impelemented yet. Ignore.");
        return 0;
    }

    if(stride != NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: writing vars at a stride is not impelemented yet. Ignore.");
        return 0;
    }

    dtf_var_t *var = fbuf->vars[varid];
    DTF_DBG(VERBOSE_DBG_LEVEL, "rw call %d for %s (ncid %d) var %d", rw_flag,fbuf->file_path, fbuf->ncid, var->id);
    for(i = 0; i < var->ndims; i++)
			DTF_DBG(VERBOSE_DBG_LEVEL, "  %lld --> %lld", start[i], count[i]);
    /*check number of elements to read*/
    nelems = 0;
    if(var->ndims == 0)
        nelems = 1;
    else
        if(count != NULL){
            int i;
            nelems = count[0];
            for(i = 1; i < var->ndims; i++)
                nelems *= count[i];

            if(nelems == 0){
                DTF_DBG(VERBOSE_DBG_LEVEL, "Nothing to read or write");
                return 0;
            }
        }
    assert(nelems != 0);

    MPI_Type_size(var->dtype, &def_el_sz);
    MPI_Type_size(dtype, &req_el_sz);

	
    if(def_el_sz != req_el_sz){
		
        DTF_DBG(VERBOSE_DBG_LEVEL, "Warning: var %d el_sz mismatch (def %d-bit, access %d).", var->id, def_el_sz, req_el_sz);
        
        /*Small hack for derived data types.*/
        if(dtype != MPI_DOUBLE && dtype != MPI_FLOAT){
			int el_sz = (int)(req_el_sz/nelems);
			assert(el_sz > 0);
			if( el_sz == def_el_sz){
				DTF_DBG(VERBOSE_DBG_LEVEL, "Ioreq dtype is a derived datatype consisting of elements of size %d, same as var type", def_el_sz);
				dtype = var->dtype;
			} else	if(el_sz == def_el_sz*2){
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: assume conversion float->double for var %d in file %s", var->id, fbuf->file_path);
				dtype = MPI_DOUBLE;
			} else if (el_sz == (int)(def_el_sz/2)){
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: assume conversion double->float for var %d in file %s", var->id, fbuf->file_path);
				dtype = MPI_FLOAT;
			} else {
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: cannot handle conversion from derived datatype for var %d in file %s", var->id, fbuf->file_path);
				MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
			}
			MPI_Type_size(dtype, &req_el_sz);
		}
		
    }
    //assert(var->dtype == dtype);

	//~ if(rw_flag == DTF_WRITE){
		//~ DTF_DBG(VERBOSE_DBG_LEVEL, "------------WRITE IOREQ--------:");
		
		//~ for(i = 0; i < nelems; i++)
			//~ printf("%.2f\t", ((double*)buf)[i]);
		//~ printf("\n");
	//~ }
        
	/*NOTE: Because dtype may be a derivative MPI type and differ from var->dtype,
	we ignore it. Start and count parameters are supposed to be with respect to
	element size for var->dtype*/
	int buffered = gl_conf.buffer_data;

	if(rw_flag == DTF_READ)
		buffered = 0;

	if( gl_scale && (var->ndims <= 1) && (rw_flag == DTF_WRITE))
		 /*This is specifically for SCALE-LETKF since they overwrite the
		  user buffer in every time frame iteration */
		buffered = 1;

	req = new_ioreq(fbuf->rreq_cnt+fbuf->wreq_cnt, varid, var->ndims, dtype, start, count, buf, rw_flag, buffered);
	
	if( gl_scale && (var->ndims <= 1) && (rw_flag == DTF_WRITE))
		 /*This is specifically for SCALE-LETKF since they overwrite the
		  user buffer in every time frame iteration */
		req->is_permanent = 1; //dont delete this req when cleaning the list of ioreqs

	if(gl_conf.do_checksum && (rw_flag == DTF_WRITE))
		var->checksum += req->checksum;

	if(rw_flag == DTF_READ)
		fbuf->rreq_cnt++;
	else
		fbuf->wreq_cnt++;
	/*Enqueue the request to the head*/
	if(var->ioreqs == NULL)
		var->ioreqs = req;
	else{
		/*Check if some data is overwritten (just to print out a warning message).
		  Becase the new I/O req is pushed to the head of the queue, the
		  writer will access the newest data.*/
		io_req_t *tmpreq = var->ioreqs;
		while(tmpreq != NULL){
			if(req->rw_flag == DTF_WRITE){
				int overlap = 0;
				for(i = 0; i < var->ndims; i++ )
					if( (req->start[i] >= tmpreq->start[i]) && (req->start[i] < tmpreq->start[i] + tmpreq->count[i]))
						overlap++;
					else
						break;

				if(overlap == var->ndims){
					DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: overwriting var %d data: (old (start,count) --> new (start,count)", var->id);
					for(i = 0; i < var->ndims; i++)
						DTF_DBG(VERBOSE_DBG_LEVEL, "(%lld, %lld) --> (%lld, %lld)", tmpreq->start[i], tmpreq->count[i], req->start[i], req->count[i]);
				}
			}
			tmpreq = tmpreq->next;
		}
		var->ioreqs->prev = req;
		req->next = var->ioreqs;
		var->ioreqs = req;
	}
    
    fbuf->has_unsent_ioreqs = 1;
    if(MPI_Wtime() - fbuf->t_last_sent_ioreqs >= gl_conf.t_send_ioreqs_freq){
		//Send request to master immediately
		if(gl_conf.iodb_build_mode == IODB_BUILD_VARID)
			send_ioreqs_by_var(fbuf);
		else //if(gl_conf.iodb_build_mode == IODB_BUILD_BLOCK)
			send_ioreqs_by_block(fbuf);
		fbuf->t_last_sent_ioreqs = MPI_Wtime();
	}

    ret = nelems*req_el_sz;
    	
    gl_stats.t_rw_var += MPI_Wtime() - t_start;
    return ret;
}

