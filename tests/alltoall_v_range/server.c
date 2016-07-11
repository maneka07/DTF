#include "farb_connect.h"
#include "common.h"

void allscatter_col_diffusion(struct server_connection* connexion,
			      int nb_groups,
			      int*buffer,
			      size_t buffer_size[2],
			      int axis) {
  /* scatter the produced data to all the clients */
  size_t block_size = buffer_size[1]/(connexion->server_size/nb_groups);
  pub_server_allscatter_2d_range_data(connexion, nb_groups, buffer, MPI_INTEGER, buffer_size, block_size, axis);
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
    size_t buffer_sizes[2] = {buffer_size, buffer_size};
    int buffer[buffer_size*buffer_size];
    int i;
    int j;

    for(i=0;i<buffer_size; i++){
      for(j=0;j<buffer_size; j++){
	buffer[(i*buffer_size)+j] = ((connexion->server_rank * buffer_size*buffer_size)+(i*buffer_size) + j);
      }
    }

#ifdef TEST_DEBUG
    sleep(connexion->server_rank);
    printf("[Server %d] \n", connexion->server_rank);
    for(i=0; i<buffer_size; i++) {
      for(j=0; j<buffer_size; j++) {
	printf("%d ", buffer[(i*buffer_size)+j]);
      }
      printf("\n");
    }
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    int axis;
    for(axis=0; axis<2; axis++) {
      printf("Sending using axis %d\n", axis);
      int nb_groups;
      for(nb_groups=1; nb_groups < connexion->client_size; nb_groups*=2) {
	printf("%d groups\n", nb_groups);

	/* warmup */
	for(i=0; i< WARMUP; i++) {
	  allscatter_col_diffusion(connexion, nb_groups, buffer, buffer_sizes, axis);
	}

	for(i=0; i<NITER; i++) {
	  allscatter_col_diffusion(connexion, nb_groups, buffer, buffer_sizes, axis);
	}
      }
    }
    conn_serv_unpublish_port(connexion);

    MPI_Barrier(MPI_COMM_WORLD);
    conn_finalize_();
    MPI_Finalize();
    return 0;
}
