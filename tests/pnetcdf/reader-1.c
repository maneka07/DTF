  /* 
 * Example created from the pnetcdf example program (pnetcdf)/examples/C/block_cyclic.c 
 * 
 * */
 
 /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example, generalized from column_wise.c, makes a number of nonblocking
 * API calls, each writes a block of columns into a 2D integer array. In other
 * words, the I/O pattern is a blocked cyclic along X dimension.
 *
 * Each process writes NX columns in total. The block length is controlled by
 * block_len. In this example, block_len is set to 2. Blocks are layout in a
 * cyclic fashion in the file. This example can test if PnetCDF can coalesce
 * file offsets and lengths when constructing a merged filetype.
 * 
 * The compile and run commands are given below, together with an ncmpidump of
 * the output file. In this example, block_len = 2.
 *
 *    % mpicc -O2 -o block_cyclic block_cyclic.c -lpnetcdf
 *    % mpiexec -l -n 4 ./block_cyclic /pvfs2/wkliao/testfile.nc
 *    0:  0: NY=10 myNX=  4 myOff=  0
 *    1:  1: NY=10 myNX=  4 myOff=  4
 *    2:  2: NY=10 myNX=  4 myOff=  8
 *    3:  3: NY=10 myNX=  4 myOff= 12
 *    0: [i=0] iput() start=  0   0 count= 10   1
 *    0: [i=1] iput() start=  0   1 count= 10   1
 *    0: [i=2] iput() start=  0   8 count= 10   1
 *    0: [i=3] iput() start=  0   9 count= 10   1
 *    1: [i=0] iput() start=  0   2 count= 10   1
 *    1: [i=1] iput() start=  0   3 count= 10   1
 *    1: [i=2] iput() start=  0  10 count= 10   1
 *    1: [i=3] iput() start=  0  11 count= 10   1
 *    2: [i=0] iput() start=  0   4 count= 10   1
 *    2: [i=1] iput() start=  0   5 count= 10   1
 *    2: [i=2] iput() start=  0  12 count= 10   1
 *    2: [i=3] iput() start=  0  13 count= 10   1
 *    3: [i=0] iput() start=  0   6 count= 10   1
 *    3: [i=1] iput() start=  0   7 count= 10   1
 *    3: [i=2] iput() start=  0  14 count= 10   1
 *    3: [i=3] iput() start=  0  15 count= 10   1
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *            Y = 10 ;
 *            X = 16 ;
 *    variables:
 *            int var(Y, X) ;
 *    data:
 *
 *     var =
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13,
 *      10, 10, 11, 11, 12, 12, 13, 13, 10, 10, 11, 11, 12, 12, 13, 13 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <unistd.h> /* getopt() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>

#include "pfarb.h"

#define NY 3
#define NX 2

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [filename] output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char** argv) {
    extern int optind;
    char *filename="restart.nc";
    int i, j, verbose=1, rank, nprocs, err, num_reqs;
    int ncid, cmode, varid, dimid[2], *reqs, *sts, **buf;
    MPI_Offset myNX, G_NX, myOff, block_start, block_len;
    MPI_Offset start[2], count[2];
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    farb_init("farb.ini", "ireader");
    
    G_NX  = NX * nprocs;
	sleep(5);
  /* open an existing file created earlier for read -----------------------*/
    cmode = NC_NOWRITE;
    err = ncmpi_open(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    ERR

    /* the global array is NY * (NX * nprocs) */
    MPI_Offset global_ny, global_nx;
    err = ncmpi_inq_dimid(ncid, "Y", &dimid[0]);
    ERR
    err = ncmpi_inq_dimlen(ncid, dimid[0], &global_ny);
    ERR
    assert(global_ny == NY);

    err = ncmpi_inq_dimid(ncid, "X", &dimid[1]);
    ERR
    err = ncmpi_inq_dimlen(ncid, dimid[1], &global_nx);
    ERR
    assert(global_nx == G_NX);

    myNX  = global_nx / nprocs;

    err = ncmpi_inq_varid(ncid, "var", &varid);
    ERR

      //      printf("app %d: varid %d, dim0 len %lld, dim1 len %lld\n", rank, varid, global_ny, global_nx);
      
    /* initialize the buffer with -1, so a read error can be pinpointed */
    buf    = (int**) malloc(myNX * sizeof(int*));
    buf[0] = (int*)  malloc(global_ny * myNX * sizeof(int));
    for (i=0; i<myNX; i++) {
        if (i > 0) buf[i] = buf[i-1] + global_ny;
        for (j=0; j<global_ny; j++) buf[i][j] = -1;
    }

    reqs = (int*) malloc(myNX * sizeof(int));
    sts  = (int*) malloc(myNX * sizeof(int));

    /* each proc reads myNX columns of the 2D array, block_len controls the
       number of contiguous columns at a time */
    block_start = 0;
    block_len   = 2;  /* can be 1, 2, 3, ..., myNX */
    if (block_len > myNX) block_len = myNX;

    start[0] = 0;          start[1] = rank * block_len;
    count[0] = global_ny;  count[1] = 1;
    num_reqs = 0;
    for (i=0; i<myNX; i++) {
        err = ncmpi_iget_vara_int(ncid, varid, start, count, buf[i],
                                  &reqs[num_reqs++]);
        ERR

        if (i % block_len == block_len-1)  {
            int stride = MIN(myNX-1-i, block_len);
            block_start += block_len * nprocs;
            start[1] = block_start + stride * rank;
        }
        else
            start[1]++;
    }
    err = ncmpi_wait_all(ncid, num_reqs, reqs, sts);
    ERR

    /* check status of all requests */
    for (i=0; i<num_reqs; i++)
        if (sts[i] != NC_NOERR)
            printf("Error: nonblocking read fails on request %d (%s)\n",
                   i, ncmpi_strerror(sts[i]));

    err = ncmpi_close(ncid);
    ERR

    /* check the read contents */
    for (i=0; i<myNX; i++) {
        for (j=0; j<global_ny; j++)
            if (buf[i][j] != rank+10) {
	      printf("%d: Read contents mismatch at buf[%d][%d] = %d (should be %d)\n", rank, i,j,buf[i][j],rank+10);
            }
    }

    free(sts);
    free(reqs);
    free(buf[0]);
    free(buf);

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }
	farb_finalize();
    MPI_Finalize();
    return 0;
}
