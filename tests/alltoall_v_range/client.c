#include <assert.h>
#include "farb_connect.h"
#include "common.h"

void allscatter_col_diffusion(struct client_connection*connexion,
			      int nb_groups,
			      int*buffer,
			      size_t buffer_size[2],
			      int axis) {
  /* receive data from all the simulation processes */
  /* This is equivalent to receiving from 1 process and calling alltoall */
  size_t block_size = buffer_size[1]/(connexion->server_size/nb_groups);
  pub_client_allscatter_2d_range_data(connexion, nb_groups, buffer, MPI_INTEGER, buffer_size, block_size, axis);
}

int main( int argc, char **argv )
{
    char service_name[] = "ocean";
    MPI_Comm server;
    struct client_connection*connexion = NULL;
    int i;
    MPI_Init(&argc, &argv);

    int buffer_size=4;
    if(argc>1) {
      buffer_size = atoi(argv[1]);
    }

    sleep(2);
    conn_init_("sample.ini", "");
    connexion = conn_client_connect(service_name, &server);

    if(connexion->client_rank == 0) {
      printf("------------------------\n");
      printf("BUFFER SIZE: %d bytes\n", buffer_size);
      printf("Warmup: %d loops\n", WARMUP);
      printf("NITER: %d\n", NITER);
      printf("NB of SIM processes: %d\n", connexion->server_size);
      printf("NB of DA processes: %d\n", connexion->client_size);
      printf("------------------------\n");
    }

    /* receive the data */
    //    char* buffer = malloc(buffer_size);
    int buffer[buffer_size*buffer_size];
    size_t buffer_sizes[2] = {buffer_size, buffer_size};
    struct timespec t1, t2;
    int j;

    int axis;
    for(axis=0; axis<2; axis++) {

      for(i=0; i<buffer_size; i++) {
	for(j=0; j<buffer_size; j++) {
	  buffer[(i*buffer_size) + j] = 0;
	}
      }

      int nb_groups;
      for(nb_groups=1; nb_groups < connexion->client_size; nb_groups*=2) {

	for(i=0; i< WARMUP; i++) {
	  allscatter_col_diffusion(connexion, nb_groups, buffer, buffer_sizes, axis);
	}

	clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
	for(i=0; i< NITER; i++) {
	  allscatter_col_diffusion(connexion, nb_groups, buffer, buffer_sizes, axis);
	}
	clock_gettime(CLOCK_MONOTONIC_RAW, &t2);

	if(connexion->client_rank == 0) {
	  double allscatter_duration=TIME_DIFF(t1, t2) / 1e6 / NITER; /* duration in ms */
	  printf("allscatter diffusion took %lf ms\n", allscatter_duration);
	}

#ifdef TEST_DEBUG
	sleep(connexion->client_rank);
	printf("[Client %d] After alltoall (nbgroups = %d, axis = %d): \n", connexion->client_rank, nb_groups, axis);
	int j;
	for(i=0; i<buffer_size; i++) {
	  for(j=0; j<buffer_size; j++) {
	    printf("%d ", buffer[(i*buffer_size)+j]);
	  }
	  printf("\n");
	}
#endif
      }
    }

    conn_client_disconnect(connexion);
    conn_finalize_();
    MPI_Finalize();

    return 0;
}
