#include <assert.h>
#include "farb_connect.h"
#include "common.h"

void two_phases_diffusion(struct client_connection*connexion,
			  char*buffer,
			  size_t buffer_size) {
    int tag = connexion->client_rank;
    int src = connexion->client_rank;
    /* receive data from the simulation process */
    conn_client_recv_data(connexion, buffer, buffer_size, tag, src);
    int sendcount = buffer_size / connexion->client_size;

    /* scatter data to the other assimilation processes */
    int err =  MPI_Alltoall(MPI_IN_PLACE, sendcount, MPI_BYTE,
			    buffer, sendcount, MPI_BYTE, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
}

void alltoall_diffusion(struct client_connection*connexion,
			char*buffer,
			size_t buffer_size) {
  /* receive data from all the simulation processes */
  /* This is equivalent to receiving from 1 process and calling alltoall */
  conn_client_alltoall_data(connexion, buffer, buffer_size/connexion->server_size);
}

void allscatter_diffusion(struct client_connection*connexion,
		       char*buffer,
		       size_t buffer_size) {
  /* receive data from all the simulation processes */
  /* This is equivalent to receiving from 1 process and calling alltoall */
  conn_client_allscatter_data(connexion, buffer, buffer_size/connexion->server_size);
}

void allscatter_nb_diffusion(struct client_connection*connexion,
			     char*buffer,
			     size_t buffer_size) {
  /* receive data from all the simulation processes */
  /* This is equivalent to receiving from 1 process and calling alltoall */
  conn_client_allscatter_nb_data(connexion, buffer, buffer_size/connexion->server_size);
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
    char* buffer = malloc(buffer_size);

    struct timespec t1, t2;

    for(i=0; i< WARMUP; i++) {
      two_phases_diffusion(connexion, buffer, buffer_size);
      alltoall_diffusion(connexion, buffer, buffer_size);
      allscatter_diffusion(connexion, buffer, buffer_size);
      allscatter_nb_diffusion(connexion, buffer, buffer_size);
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
    for(i=0; i< NITER; i++) {
      two_phases_diffusion(connexion, buffer, buffer_size);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &t2);

    if(connexion->client_rank == 0) {
      double two_phases_duration=TIME_DIFF(t1, t2) / 1e6 / NITER; /* duration in ms */
      printf("two phases diffusion took %lf ms\n", two_phases_duration);
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
    for(i=0; i< NITER; i++) {
      alltoall_diffusion(connexion, buffer, buffer_size);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &t2);

    if(connexion->client_rank == 0) {
      double alltoall_duration=TIME_DIFF(t1, t2) / 1e6 / NITER; /* duration in ms */
      printf("alltoall diffusion took %lf ms\n", alltoall_duration);
    }


    clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
    for(i=0; i< NITER; i++) {
      allscatter_diffusion(connexion, buffer, buffer_size);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &t2);

    if(connexion->client_rank == 0) {
      double allscatter_duration=TIME_DIFF(t1, t2) / 1e6 / NITER; /* duration in ms */
      printf("allscatter diffusion took %lf ms\n", allscatter_duration);
    }


    clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
    for(i=0; i< NITER; i++) {
      allscatter_nb_diffusion(connexion, buffer, buffer_size);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &t2);

    if(connexion->client_rank == 0) {
      double allscatter_nb_duration=TIME_DIFF(t1, t2) / 1e6 / NITER; /* duration in ms */
      printf("allscatter NB diffusion took %lf ms\n", allscatter_nb_duration);
    }

#ifdef TEST_DEBUG
    printf("[Client %d] After alltoall: \t", connexion->client_rank);
    for(i=0; i<buffer_size; i++) {
      printf("%d ", buffer[i]);
    }
    printf("\n");
#endif

    conn_client_disconnect(connexion);
    conn_finalize_();
    MPI_Finalize();

    return 0;
}
