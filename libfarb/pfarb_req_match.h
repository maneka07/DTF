#ifndef PFARB_REQ_MATCH_H_INCLUDED
#define PFARB_REQ_MATCH_H_INCLUDED

#include <mpi.h>

typedef struct io_req{
    void *user_buf;
    MPI_Offset *start;
    MPI_Offset *count;
    struct io_req *next;
}io_req_t;

io_req_t *new_ioreq(int ndims, MPI_Offset *start, MPI_Offset *count, void *buf);
void add_ioreq(io_req_t *list, io_req_t *ioreq);
int match_ioreqs(const char* filename);

#endif // PFARB_REQ_MATCH_H_INCLUDED
