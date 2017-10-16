#include "dtf_file_buffer.h"
#include "dtf_nbuf_io.h"
#include "dtf_util.h"
#include "dtf.h"
#include "dtf_req_match.h"

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
    int def_el_sz, req_el_sz;
    int type_mismatch = 0;
    MPI_Offset nelems;
    double t_start = MPI_Wtime();

    if(rw_flag == DTF_READ){
        if(fbuf->reader_id==gl_my_comp_id){
            if(!fbuf->is_ready){
                DTF_DBG(VERBOSE_DBG_LEVEL, "DTF Error trying to read file %s that is not ready", fbuf->file_path);
                assert(fbuf->is_ready);
            }
        } else{
            DTF_DBG(VERBOSE_WARNING_LEVEL, "DTF Warning: writer process tries to read file %s (var %d)", fbuf->file_path, varid);
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
        type_mismatch = 1;
        DTF_DBG(VERBOSE_DBG_LEVEL, "Warning: var %d el_sz mismatch (def %d-bit, access %d).", var->id, def_el_sz, req_el_sz);
    }
    //assert(var->dtype == dtype);

    DTF_DBG(VERBOSE_DBG_LEVEL, "------------IOREQ--------:");
    for(i = 0; i < var->ndims; i++)
        DTF_DBG(VERBOSE_DBG_LEVEL, "  %lld --> %lld", start[i], count[i]);

    /*If the process has the data to match a read request completely, then
    copy the data and do not create an I/O request.*/
    int create_ioreq = 1;

    if((rw_flag == DTF_READ) && (fbuf->writer_id == gl_my_comp_id)){

        if(var->ndims > 0){
            int match;
            MPI_Offset matched_els;
            MPI_Offset *_start, *_count;
            _start = dtf_malloc(sizeof(MPI_Offset)*var->ndims);
            assert(_start != NULL);
            _count = dtf_malloc(sizeof(MPI_Offset)*var->ndims);
            assert(_count != NULL);
            memcpy(_start, start, sizeof(MPI_Offset)*var->ndims);
            memcpy(_count, count, sizeof(MPI_Offset)*var->ndims);

            matched_els = 0;
            while(nelems > 0){
                io_req_t *tmp = fbuf->ioreqs;
                while(tmp != NULL){
                    if( (tmp->var_id == varid) && (tmp->rw_flag == DTF_WRITE) ){
                        match = 0;
                        for(i = 0; i < var->ndims; i++)
                            if( (_start[i] >= tmp->start[i]) && (_start[i] < tmp->start[i] + tmp->count[i]))
                                match++;
                            else
                                break;

                        if(match == var->ndims)
                            break;
                    }
                    tmp = tmp->next;
                }
                if((tmp != NULL) && (match == var->ndims)){
                    MPI_Offset *coord = dtf_malloc(var->ndims*sizeof(MPI_Offset));
                    assert(coord != NULL);
                    matched_els = 1;
                    for(i = 0; i < var->ndims; i++){
                        coord[i] = _start[i];
                        if(_start[i] + _count[i] > tmp->start[i] + tmp->count[i])
                            _count[i] = tmp->start[i] + tmp->count[i] - _start[i];
                        matched_els *= _count[i];
                    }
                    recur_get_put_data(var, tmp->dtype, tmp->user_buf, tmp->start, tmp->count, start, count, 0, coord, buf, DTF_READ, type_mismatch);
                    dtf_free(coord, var->ndims*sizeof(MPI_Offset));

                    nelems -= matched_els;
                    /*adjust start coord*/
                    shift_coord(var->ndims, start, count, _start, _count, matched_els);
                } else
                    break;
            } /*while*/
            dtf_free(_start, sizeof(MPI_Offset)*var->ndims);
            dtf_free(_count, sizeof(MPI_Offset)*var->ndims);

        } else {
            io_req_t *tmp = fbuf->ioreqs;
            while(tmp != NULL){
                if( (tmp->var_id == varid) && (tmp->rw_flag == DTF_WRITE) ){

                        assert(tmp->user_buf_sz == (MPI_Offset)def_el_sz);
                        memcpy(buf, tmp->user_buf, (size_t)def_el_sz);
                        nelems--;
                        break;
                }
            }
        }
        if(nelems == 0){
            DTF_DBG(VERBOSE_DBG_LEVEL, "Have read data for var %d, no rreq created", varid);
            create_ioreq = 0;

            if(var->ndims == 1){
                DTF_DBG(VERBOSE_DBG_LEVEL, "Var data: ");
                for(i = 0; i < count[0]; i++)
                    printf("%.3f ", ((double*)(buf))[i]);
                printf("\n");
            }
        } else
            DTF_DBG(VERBOSE_DBG_LEVEL, "Have not read all data for var %d, rreq created", varid);
    }


    if(create_ioreq){
        int sclltkf;
        /*NOTE: Because dtype may be a derivative MPI type and differ from var->dtype,
        we ignore it. Start and count parameters are supposed to be with respect to
        element size for var->dtype*/
        int buffered = gl_conf.buffered_req_match;

        if(rw_flag == DTF_READ)
            buffered = 0;

        char *s = getenv("DTF_SCALE");
        if(s != NULL)
            sclltkf = atoi(s);
        else
            sclltkf = 0;

        if( sclltkf && (var->ndims <= 1) && (rw_flag == DTF_WRITE))
             /*This is specifically for SCALE-LETKF since they overwrite the
              user buffer in every time frame iteration */
            buffered = 1;

        req = new_ioreq(fbuf->rreq_cnt+fbuf->wreq_cnt, varid, var->ndims, dtype, start, count, buf, rw_flag, buffered);
        
        req->is_permanent = 1; //dont delete this req when cleaning the list of ioreqs

        if(gl_conf.do_checksum && (rw_flag == DTF_WRITE))
            var->checksum += req->checksum;

        if(request != NULL)
            *request = req->id;
        if(rw_flag == DTF_READ)
            fbuf->rreq_cnt++;
        else
            fbuf->wreq_cnt++;

        /*Enqueue the request to the head*/
        if(fbuf->ioreqs == NULL)
            fbuf->ioreqs = req;
        else{
            /*Check if some data is overwritten (just to print out a warning message).
              Becase the new I/O req is pushed to the head of the queue, the
              writer will access the newest data.*/
            io_req_t *tmpreq = fbuf->ioreqs;
            while(tmpreq != NULL){
                if( (req->rw_flag == DTF_WRITE) && (tmpreq->var_id == req->var_id)){
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
            fbuf->ioreqs->prev = req;
            req->next = fbuf->ioreqs;
            fbuf->ioreqs = req;
        }
    } /*if(create_ioreq)*/

    MPI_Type_size(dtype, &el_sz);
    ret = 1;
    for(i = 0; i < var->ndims; i++)
        ret *= count[i];
    ret *= el_sz;

    gl_stats.accum_rw_var += MPI_Wtime() - t_start;
    return ret;
}
