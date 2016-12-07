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

    if(rw_flag == FARB_READ)
         assert(fbuf->is_ready);

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
    FARB_DBG(VERBOSE_DBG_LEVEL, "rw %d call for ncid %d var %d", rw_flag, fbuf->ncid, var->id);
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
    MPI_Type_size(dtype, &el_sz);
    assert(el_sz > 0);
    if(var->el_sz == 0){
        var->el_sz = (MPI_Offset)el_sz;
    } else {
        if(var->el_sz != (MPI_Offset)el_sz)
        FARB_DBG(VERBOSE_DBG_LEVEL, "Warning: el_sz mismatch (written as %llu-bit,read as %d-bit var). Using the original size.", var->el_sz, el_sz);
        //assert(var->el_sz == (MPI_Offset)el_sz);
    }

    /*Create an io request*/
    req = new_ioreq(fbuf->ioreq_cnt, varid, var->ndims, var->el_sz, start, count, buf, rw_flag, gl_conf.buffered_req_match);
    if(request != NULL)
        *request = req->id;
    fbuf->ioreq_cnt++;
    get_contig_mem_list(var, start, count, &(req->nchunks), &(req->mem_chunks));
    assert(req->mem_chunks != NULL);

    /*Enqueue the request*/
    if(fbuf->ioreqs == NULL)
        fbuf->ioreqs = req;
    else{
        fbuf->ioreqs->prev = req;
        req->next = fbuf->ioreqs;
        fbuf->ioreqs = req;
    }

    ret = 1;
    for(i = 0; i < var->ndims; i++)
        ret *= count[i];
    ret *= var->el_sz;

    return ret;
}
