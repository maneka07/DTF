#include "pfarb_file_buffer.h"
#include "pfarb_nbuf_io.h"
#include "pfarb_util.h"
#include "pfarb.h"
#include "pfarb_mem.h"
#include "pfarb_common.h"

MPI_Offset nbuf_read_write_var(file_buffer_t *fbuf,
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
    int el_sz;
    io_req_t *req;
    int i;

    if(rw_flag == FARB_READ){
        if(fbuf->reader_id==gl_my_comp_id){
            if(!fbuf->is_ready){
                FARB_DBG(VERBOSE_DBG_LEVEL, "FARB Error trying to read file %s that is not ready", fbuf->file_path);
                assert(fbuf->is_ready);
            }
        } else{
            FARB_DBG(VERBOSE_WARNING_LEVEL, "FARB Warning: writer process tries to read file %s (var %d)", fbuf->file_path, varid);
            //MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }
    }
    if(imap != NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: writing mapped vars is not impelemented yet. Ignore.");
        return 0;
    }

    if(stride != NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: writing vars at a stride is not impelemented yet. Ignore.");
        return 0;
    }

    farb_var_t *var = find_var(fbuf->vars, varid);
    if(var == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: could not find var with id %d", varid);
        return 0;
    }
    FARB_DBG(VERBOSE_DBG_LEVEL, "rw call %d for %s (ncid %d) var %d", rw_flag,fbuf->file_path, fbuf->ncid, var->id);
    /*check number of elements to read*/
    if(count != NULL){
        MPI_Offset nelems;
        int i;
        nelems = count[0];
        for(i = 1; i < var->ndims; i++)
            nelems *= count[i];

        if(nelems == 0){
            FARB_DBG(VERBOSE_DBG_LEVEL, "Nothing to read or write");
            return 0;
        }

    }





    if(var->dtype != dtype){
        int el_sz1, el_sz2;
        MPI_Type_size(var->dtype, &el_sz1);
        MPI_Type_size(dtype, &el_sz2);
        FARB_DBG(VERBOSE_ALL_LEVEL, "Warning: el_sz mismatch (defined as %d-bit, accessed as %d-bit var). Using the original size.", el_sz1, el_sz2);
    }
    //assert(var->dtype == dtype);


    for(i = 0; i < var->ndims; i++)
        FARB_DBG(VERBOSE_DBG_LEVEL, "RW REQ: -> start %d, count %d", (int)start[i], (int)count[i]);

    /*Create an io request*/

    if(frt_indexing){
        FARB_DBG(VERBOSE_DBG_LEVEL, "FARB Warning: reversing indeces because fortran!");
        MPI_Offset *new_start = farb_malloc(var->ndims*sizeof(MPI_Offset));
        assert(new_start != NULL);
        MPI_Offset *new_count = farb_malloc(var->ndims*sizeof(MPI_Offset));
        assert(new_count != NULL);

        for(i = 0; i < var->ndims; i++){
            new_count[var->ndims-1-i] = count[i];
            new_start[var->ndims-1-i] = start[i];
        }
        for(i = 0; i < var->ndims; i++)
            FARB_DBG(VERBOSE_DBG_LEVEL, "RW REQ REVERSE: -> start %d, count %d", (int)new_start[i], (int)new_count[i]);
        req = new_ioreq(fbuf->rreq_cnt+fbuf->wreq_cnt, varid, var->ndims, dtype, new_start, new_count, buf, rw_flag, gl_conf.buffered_req_match);
        farb_free(new_start, var->ndims*sizeof(MPI_Offset));
        farb_free(new_count, var->ndims*sizeof(MPI_Offset));
    } else
        req = new_ioreq(fbuf->rreq_cnt+fbuf->wreq_cnt, varid, var->ndims, dtype, start, count, buf, rw_flag, gl_conf.buffered_req_match);

    if(request != NULL)
        *request = req->id;
    if(rw_flag == FARB_READ)
        fbuf->rreq_cnt++;
    else
        fbuf->wreq_cnt++;

    if(gl_conf.io_db_type == FARB_DB_CHUNKS){
        get_contig_mem_list(var, dtype, start, count, &(req->nchunks), &(req->mem_chunks));
        FARB_DBG(VERBOSE_DBG_LEVEL, "RW_REQ: nchunks %d ------------------", req->nchunks);
        assert(req->mem_chunks != NULL);
    } else {
        req->nchunks = 0;
        req->mem_chunks = NULL;
    }


    /*Enqueue the request*/
    if(fbuf->ioreqs == NULL)
        fbuf->ioreqs = req;
    else{
        fbuf->ioreqs->prev = req;
        req->next = fbuf->ioreqs;
        fbuf->ioreqs = req;
    }

    MPI_Type_size(dtype, &el_sz);
    ret = 1;
    for(i = 0; i < var->ndims; i++)
        ret *= count[i];
    ret *= el_sz;

    return ret;
}
