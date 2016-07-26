/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This is a simple example which reads a small dummy array, which was
   written by simple_xy_wr.c. This is intended to illustrate the use
   of the netCDF C API.

   This program is part of the netCDF tutorial:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

   Full documentation of the netCDF C API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c

   $Id: simple_xy_rd.c,v 1.9 2006/08/17 23:00:55 russ Exp $
*/
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <mpi.h>
#include "farb.h"

/* This is the name of the data file we will read. */
#define FILE_NAME "simple_xy.nc"

/* We are reading 2D data, a 6 x 12 grid. */
#define NX 2
#define NY 3

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int
main()
{
  MPI_Init(0, NULL);
  farb_init("farb.ini", "ireader");
	
   /* This will be the netCDF ID for the file and data variable. */
  int ncid;
  int varid;
  int myrank, nranks;
  //int data_in[NX][NY];
  int *data_in, *data_in2;
  char filename[128];
  /* Loop indexes, and error handling. */
  int x, y, retval;
  
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);
   sprintf(filename, "restart_%d.nc", myrank);
   data_in = malloc(nranks * NY * sizeof(int));
   data_in2 = malloc(nranks * NY * sizeof(int));
   /* Open the file. NC_NOWRITE tells netCDF we want read-only access
    * to the file.*/
   if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
      ERR(retval);
   printf("%d: opened ncid %d\n", myrank, ncid);
   /* Get the varid of the data variable, based on its name. */
   if ((retval = nc_inq_varid(ncid, "data", &varid)))
      ERR(retval);


     /* Read the data. */
     //   if ((retval = nc_get_var_int(ncid, varid, &data_in[0][0])))
   if ((retval = nc_get_var_int(ncid, varid, data_in)))
     ERR(retval);
     /* Check the data. */
     /*   for (x = 0; x < NX; x++)
      for (y = 0; y < NY; y++)
      if (data_in[x][y] != x * NY + y)
      return ERRCODE;*/
   for (x = 0; x < nranks; x++)
     for (y = 0; y < NY; y++)
       printf("r%d:%d ",myrank, data_in[x*NY+y]);
   printf("\n");
   
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Alltoall(data_in, NY, MPI_INT, data_in2, NY, MPI_INT, MPI_COMM_WORLD);
     
   for (x = 0; x < nranks; x++)
     for (y = 0; y < NY; y++)
       printf("rr%d:%d ",myrank, data_in2[x*NY+y]);
   printf("\n");

   /* Close the file, freeing all resources. */
   if ((retval = nc_close(ncid)))
     ERR(retval);

   
   printf("*** SUCCESS reading example file %s!\n", filename);
   farb_finalize();
   MPI_Finalize();
   return 0;
}
