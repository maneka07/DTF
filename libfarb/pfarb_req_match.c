#include "pfarb_req_match.h"
#include "pfarb_common.h"
#include "pfarb_util.h"

io_req_t *new_ioreq(int ndims, MPI_Offset *start, MPI_Offset *count, void *buf){

    io_req_t *ioreq = (io_req_t*)malloc(sizeof(io_req_t));
    assert(ioreq != NULL);
    if(ndims > 0){
        ioreq->start = (MPI_Offset*)malloc(sizeof(MPI_Offset)*ndims);
        assert(ioreq->start != NULL);
        memcpy((void*)ioreq->start, (void*)start, sizeof(MPI_Offset)*ndims);
        ioreq->count = (MPI_Offset*)malloc(sizeof(MPI_Offset)*ndims);
        assert(ioreq->count != NULL);
        memcpy((void*)ioreq->count, (void*)count, sizeof(MPI_Offset)*ndims);
    } else{
        ioreq->start = NULL;
        ioreq->count = NULL;
    }
    ioreq->user_buf = buf;
    ioreq->next = NULL;

    return ioreq;
}

void add_ioreq(io_req_t *list, io_req_t *ioreq)
{
    if(list == NULL)
        list = ioreq;
    else{
        io_req_t *tmp = list;
        while(tmp->next!=NULL)
            tmp=tmp->next;
        tmp->next = ioreq;
    }
}

void build_ioreq_db(file_buffer_t *fbuf)
{

}

void send_ioreq_list(file_buffer_t *fbuf)
{
    int bufsz = 0;
    int offt = 0;
    void *sbuf = malloc(1024);
    assert(sbuf != NULL);
    farb_var_t *var = fbuf->vars;

    while(var != NULL){
        if(var->ioreq_cnt == 0)
            continue;

        //TODO continue
        var = var->next;
    }
}

int match_ioreqs(const char* filename)
{
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);

    if(fbuf->writer_id == gl_my_comp_id){

        if(gl_conf.my_master == gl_my_rank) //i am master
            build_ioreq_db(fbuf);
        else
            send_ioreq_list(fbuf);

        while(!fbuf->match_completed)
            progress_io();
    } else if(fbuf->reader_id == gl_my_comp_id){
        //TODO todo
    }
    return 0;
}
