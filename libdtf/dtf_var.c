#include "dtf_util.h"
#include "dtf_var.h"
#include "dtf_req_match.h"
#include "dtf.h"

static dtype_params_t *get_dtype_params(MPI_Datatype dtype, int ndims, /*out*/ MPI_Datatype *eltype)
{
	int i, err;
	int *array_of_ints = NULL; 
	MPI_Aint *array_of_adds = NULL; 
	MPI_Datatype *array_of_dtypes = NULL; 
	int num_ints=0, num_adds=0, num_dtypes=0, combiner; 
	dtype_params_t *params = NULL;
	*eltype = MPI_DATATYPE_NULL;
	
	err = MPI_Type_get_envelope( dtype, &num_ints, &num_adds, &num_dtypes, &combiner); 
	CHECK_MPI(err);
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "nint %d, nadd %d, ndtype %d, combiner %d", num_ints, num_adds, num_dtypes, combiner);
	
	if(combiner == MPI_COMBINER_NAMED)
		return NULL;
	else if(combiner != MPI_COMBINER_SUBARRAY){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: do not support such derived datatype (currently only support MPI_COMBINER_SUBARRAY)"); 
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	} 
	
	if(num_dtypes != 1){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: derived datatype consists of %d datatypes. Can't handle this.", num_dtypes);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	}
	
	if(num_adds > 0)
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Warning: dispacement addresses are provided in derived data type. Data extraction result might be incorrect.");
	assert(num_ints == ndims*3+2);
	array_of_ints   = (int *)dtf_malloc( num_ints * sizeof(int) ); 
	array_of_adds   =  (MPI_Aint *) dtf_malloc( num_adds * sizeof(MPI_Aint) ); 
	array_of_dtypes = (MPI_Datatype *)dtf_malloc( num_dtypes * sizeof(MPI_Datatype) ); 
   
	err = MPI_Type_get_contents( dtype, num_ints, num_adds, num_dtypes, 
					 array_of_ints, array_of_adds, array_of_dtypes ); 
	CHECK_MPI(err);
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "I/O call for a derived datatype - subblock of array of size:");
	for(i = 0; i < num_ints; i++)
		DTF_DBG(VERBOSE_DBG_LEVEL, "%d", array_of_ints[i]);
	
	params = dtf_malloc(sizeof(dtype_params_t));
	params->orig_array_size = dtf_malloc(ndims*sizeof(MPI_Offset));
	params->orig_start = dtf_malloc(ndims*sizeof(MPI_Offset));
		
	assert(array_of_ints[0] == ndims);
	
	if(array_of_ints[num_ints-1] == MPI_ORDER_C){
		for(i = 0; i < ndims; i++){
			params->orig_array_size[i] = (MPI_Offset)array_of_ints[i+1];
			params->orig_start[i] = (MPI_Offset)array_of_ints[i+1+ndims*2];
		}
	} else { /*MPI_ORDER_FORTRAN*/
		for(i = 0; i < ndims; i++){
			params->orig_array_size[i] = (MPI_Offset)array_of_ints[num_ints-2-ndims*2 - i];
			params->orig_start[i] = (MPI_Offset)array_of_ints[num_ints-2-i];
		}
	}
	DTF_DBG(VERBOSE_DBG_LEVEL, "Original size -> subarray start from:");
	for(i = 0; i < ndims; i++)
		DTF_DBG(VERBOSE_DBG_LEVEL, "%lld --> %lld", params->orig_array_size[i], params->orig_start[i]);
	
	//TODO temporary solution. what about other fortran->ctype conversions?
	*eltype = array_of_dtypes[0];
	if(*eltype == MPI_DOUBLE_PRECISION) *eltype = MPI_DOUBLE;
	else if(*eltype == MPI_REAL) *eltype = MPI_FLOAT;
	
	dtf_free(array_of_ints, num_ints*sizeof(int));
	dtf_free(array_of_dtypes, num_dtypes * sizeof(MPI_Datatype));
	dtf_free(array_of_adds, num_adds * sizeof(MPI_Aint));
	
	/*Finally, make sure that the eltype is not a derived 
	 * datatype itself (don't hande this situation)*/
	num_ints=num_dtypes=num_adds=0;
	err = MPI_Type_get_envelope( *eltype, &num_ints, &num_adds, &num_dtypes, &combiner); 
	CHECK_MPI(err);
	if(combiner != MPI_COMBINER_NAMED){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: Do not support derived datatype that consists of elements of a derived datatype.\n");
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	}
	
	return params;
}

int def_var(struct file_buffer *fbuf, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int el_sz;
    int i;
    char typename[1024];
    int len = 0;
    dtf_var_t *var = new_var(varid, ndims, dtype, shape);
    add_var(fbuf, var);
    
    MPI_Type_size(dtype, &el_sz);
	MPI_Type_get_name(dtype, typename, &len);
    
    DTF_DBG(VERBOSE_DBG_LEVEL, "varid %d, dim %d, type %s, el_sz %d. shape:", varid, ndims, typename, el_sz);
    for(i = 0; i < ndims; i++)
        DTF_DBG(VERBOSE_DBG_LEVEL, "\t%lld", shape[i]);

    return 0;
}


dtf_var_t* new_var(int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape)
{
    int i;

    dtf_var_t *var = (dtf_var_t*)dtf_malloc(sizeof(dtf_var_t));
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
    dtype_params_t *derived_params = NULL;

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
	
	MPI_Type_size(dtype, &req_el_sz);
	
    if(var->dtype != dtype){
		MPI_Datatype eltype = MPI_DATATYPE_NULL;
		MPI_Type_size(var->dtype, &def_el_sz);
		DTF_DBG(VERBOSE_DBG_LEVEL, "Warning: var %d el_sz mismatch (def %d-bit, access %d).", var->id, def_el_sz, req_el_sz);
		derived_params = get_dtype_params(dtype, var->ndims, &eltype); 
		if(derived_params != NULL){
			/*Change original derived data type to the data type of its element*/
			dtype = eltype;
			MPI_Type_size(eltype, &req_el_sz);
		}
	}

	int buffered = gl_conf.buffer_data;

	if(rw_flag == DTF_READ)
		buffered = 0;

	req = new_ioreq(fbuf->rreq_cnt+fbuf->wreq_cnt, var->ndims, dtype, start, count, derived_params, buf, rw_flag, buffered);
	
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
		//~ io_req_t *tmpreq = var->ioreqs;
		//~ while(tmpreq != NULL){
			//~ if(req->rw_flag == DTF_WRITE){
				//~ int overlap = 0;
				//~ for(i = 0; i < var->ndims; i++ )
					//~ if( (req->start[i] >= tmpreq->start[i]) && (req->start[i] < tmpreq->start[i] + tmpreq->count[i]))
						//~ overlap++;
					//~ else
						//~ break;

				//~ if(overlap == var->ndims){
					//~ DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Warning: overwriting var %d data: (old (start,count) --> new (start,count)", var->id);
					//~ for(i = 0; i < var->ndims; i++)
						//~ DTF_DBG(VERBOSE_DBG_LEVEL, "(%lld, %lld) --> (%lld, %lld)", tmpreq->start[i], tmpreq->count[i], req->start[i], req->count[i]);
				//~ }
			//~ }
			//~ tmpreq = tmpreq->next;
		//~ }
		var->ioreqs->prev = req;
		req->next = var->ioreqs;
		var->ioreqs = req;
	}
	
    if(fbuf->t_last_sent_ioreqs == 0)
		fbuf->t_last_sent_ioreqs = MPI_Wtime();
		
    if(MPI_Wtime() - fbuf->t_last_sent_ioreqs >= gl_conf.t_send_ioreqs_freq){
		//Send request to master immediately
		if(gl_conf.iodb_build_mode == IODB_BUILD_VARID)
			send_ioreqs_by_var(fbuf);
		else //if(gl_conf.iodb_build_mode == IODB_BUILD_BLOCK)
			send_ioreqs_by_block(fbuf);
		fbuf->t_last_sent_ioreqs = MPI_Wtime();
	}

    ret = nelems*req_el_sz;
    
   // dtf_log_ioreq(fbuf->file_path, varid, var->ndims, start, count, dtype, buf, rw_flag);	
                                        	
    gl_stats.t_rw_var += MPI_Wtime() - t_start;
    return ret;
}

