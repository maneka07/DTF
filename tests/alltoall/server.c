#include "farb_connect.h"
#include "common.h"

void two_phases_diffusion(struct server_connection* connexion,
			  char*buffer,
			  size_t buffer_size) {
    /* send the produced data to one client */
    int k = connexion->server_rank;
    int tag = k;
    int dest = k;
    conn_server_send_data(connexion, buffer, buffer_size, tag, dest);
}

void alltoall_diffusion(struct server_connection* connexion,
			  char*buffer,
			size_t buffer_size) {
  /* scatter the produced data to all the clients */
  conn_server_alltoall_data(connexion, buffer, buffer_size/connexion->client_size);
}


void allscatter_diffusion(struct server_connection* connexion,
			  char*buffer,
			  size_t buffer_size) {
  /* scatter the produced data to all the clients */
  conn_server_allscatter_data(connexion, buffer, buffer_size/connexion->client_size);
}

void allscatter_nb_diffusion(struct server_connection* connexion,
			  char*buffer,
			  size_t buffer_size) {
  /* scatter the produced data to all the clients */
  conn_server_allscatter_nb_data(connexion, buffer, buffer_size/connexion->client_size);
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

    /* warmup */
    for(i=0; i< WARMUP; i++) {
      two_phases_diffusion(connexion, buffer, buffer_size);
      alltoall_diffusion(connexion, buffer, buffer_size);
      allscatter_diffusion(connexion, buffer, buffer_size);
      allscatter_nb_diffusion(connexion, buffer, buffer_size);
    }

    for(i=0; i<NITER; i++) {
      two_phases_diffusion(connexion, buffer, buffer_size);
    }

    for(i=0; i<NITER; i++) {
      alltoall_diffusion(connexion, buffer, buffer_size);
    }

    for(i=0; i<NITER; i++) {
      allscatter_diffusion(connexion, buffer, buffer_size);
    }

    for(i=0; i<NITER; i++) {
      allscatter_nb_diffusion(connexion, buffer, buffer_size);
    }

    conn_serv_unpublish_port(connexion);

    MPI_Barrier(MPI_COMM_WORLD);
    conn_finalize_();
    MPI_Finalize();
    return 0;
}
