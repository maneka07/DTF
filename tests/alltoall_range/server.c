#include "farb_connect.h"
#include "common.h"

void allscatter_diffusion(struct server_connection* connexion,
			  int nb_groups,
			  char*buffer,
			  size_t buffer_size) {
  /* scatter the produced data to all the clients */
  size_t block_size = buffer_size/(connexion->server_size/nb_groups);
  pub_server_allscatter_range_data(connexion, nb_groups, buffer, block_size);
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
    char* buffer = malloc(buffer_size);
    int i;
    for(i=0;i<buffer_size; i++){
      buffer[i] = (char)((connexion->server_rank * buffer_size)+i);
    }

#ifdef TEST_DEBUG
    printf("[Server %d] \t", connexion->server_rank);
    int j;
    for(j=0; j<buffer_size; j++) {
      printf("%d ", buffer[j]);
    }
    printf("\n");
#endif

    int nb_groups;
    for(nb_groups=1; nb_groups < connexion->client_size; nb_groups*=2) {
      printf("%d groups\n", nb_groups);

      /* warmup */
      for(i=0; i< WARMUP; i++) {
	allscatter_diffusion(connexion, nb_groups, buffer, buffer_size);
      }

      for(i=0; i<NITER; i++) {
	allscatter_diffusion(connexion, nb_groups, buffer, buffer_size);
      }
    }

    conn_serv_unpublish_port(connexion);

    MPI_Barrier(MPI_COMM_WORLD);
    conn_finalize_();
    MPI_Finalize();
    return 0;
}
