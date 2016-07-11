#include "farb_connect.h"
#include "common.h"

void allscatter_3d_diffusion(struct server_connection* connexion,
			      int nb_groups,
			      int*buffer,
			      size_t buffer_size[3],
			      int axis) {
  /* scatter the produced data to all the clients */
  size_t block_size = buffer_size[1]/(connexion->server_size/nb_groups);
  pub_server_allscatter_3d_range_data(connexion, nb_groups, buffer, MPI_INTEGER, buffer_size, block_size, axis);
}

int main( int argc, char **argv ){
    char service_name[] = "ocean";
    char port_name[MPI_MAX_PORT_NAME];
    MPI_Comm client;
    struct server_connection * connexion = NULL;
    MPI_Init(&argc, &argv);

    size_t buffer_size=4;
    if(argc>1) {
      buffer_size = atoi(argv[1]);
    }

    conn_init_("sample.ini", "");
    connexion = conn_serv_publish_port(service_name, port_name, &client);

    /* do something */
    size_t buffer_sizes[3] = {buffer_size, buffer_size, buffer_size};
    int buffer[buffer_size*buffer_size*buffer_size];
    int i, j, k;

    for(i=0;i<buffer_size; i++){
      for(j=0;j<buffer_size; j++){
	for(k=0; k<buffer_size; k++){
	  buffer[(i*buffer_size*buffer_size)+(j*buffer_size)+k] = ((connexion->server_rank * buffer_size*buffer_size*buffer_size)+(i*buffer_size*buffer_size) + (j*buffer_size)+k);
	}
      }
    }

#ifdef TEST_DEBUG
    sleep(connexion->server_rank);
    printf("[Server %d] \n", connexion->server_rank);
    for(i=0; i<buffer_size; i++) {
      printf("i=%d\n", i);
      for(j=0; j<buffer_size; j++) {
	for(k=0; k<buffer_size; k++) {
	  printf("%d ", buffer[(i*buffer_size*buffer_size)+(j*buffer_size)+k]);
	}
	printf("\n");
      }
    }

    if(connexion->server_rank == 0) {

      int index = 2;
      int v1 = buffer[index];
      printf("buffer[%d] = %d\n", index, v1);
      int index2 = index + buffer_sizes[0]*buffer_sizes[1];
      int v2 = buffer[index2];
      printf("buffer[%d] = %d\n", index2, v2);
    }

#endif

    MPI_Barrier(MPI_COMM_WORLD);
    int axis;
    for(axis=0; axis<3; axis++) {
      printf("Sending using axis %d\n", axis);
      int nb_groups;
      for(nb_groups=1; nb_groups < connexion->client_size; nb_groups*=2) {
	printf("%d groups\n", nb_groups);

	/* warmup */
	for(i=0; i< WARMUP; i++) {
	  allscatter_3d_diffusion(connexion, nb_groups, buffer, buffer_sizes, axis);
	}

	for(i=0; i<NITER; i++) {
	  allscatter_3d_diffusion(connexion, nb_groups, buffer, buffer_sizes, axis);
	}
      }
    }
    conn_serv_unpublish_port(connexion);

    MPI_Barrier(MPI_COMM_WORLD);
    conn_finalize_();
    MPI_Finalize();
    return 0;
}
