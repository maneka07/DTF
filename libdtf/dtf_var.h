#ifndef DTF_VAR_H_INCLUDED
#define DTF_VAR_H_INCLUDED

#include <mpi.h>

struct file_buffer;
struct io_req;

typedef struct dtf_var{
    int                     id;         /* varid assigned by pnetcdf*/
    MPI_Datatype            dtype;      /*Datatype of the variable*/
    MPI_Offset              *shape;     /* dim->size of each dim */
    int                     ndims;      /* number of dimensions */
    struct io_req           *ioreqs;    /*Read or write I/O requests*/
    double                  checksum;   /*to check if what was written to the var is what was read*/
 }dtf_var_t;

dtf_var_t* 	new_var(int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
void 		add_var(struct file_buffer *fbuf, dtf_var_t *var);
MPI_Offset 	read_write_var(struct file_buffer *fbuf,
								   int varid,
								   const MPI_Offset *start,
								   const MPI_Offset *count,
								   const MPI_Offset *stride,
								   const MPI_Offset *imap,
								   MPI_Datatype dtype,
								   void *buf,
								   int rw_flag);
int 		def_var(struct file_buffer *fbuf, int varid, int ndims, MPI_Datatype dtype, MPI_Offset *shape);
void 		delete_var(struct file_buffer *fbuf, dtf_var_t* var);

#endif
