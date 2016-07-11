#include "farb_connect.h"

int main( int argc, char **argv ){
    char service_name[] = "ocean";
    char port_name[MPI_MAX_PORT_NAME];
    MPI_Comm client;
    struct server_connection * connexion = NULL;
    MPI_Init(&argc, &argv);

    conn_init_("sample.ini", "");

    connexion = conn_serv_publish_port(service_name, port_name, &client);

    printf("port: %s\n", connexion->port_name);

    /* do something */
    int buffer[32];
    int i;
    for(i=0;i<32; i++){
      buffer[i] = connexion->server_rank;
    }
    sleep(1);

    /* send the produced data */
    int k = connexion->server_rank;
    for(i=0; i<connexion->client_size; i++) {
      conn_server_send_data(connexion, &k, sizeof(int), i, i);
    }


    /* send 4096 int to each client (sending INTs once at a time) */
    /* Now use pack send */
    for(i=0; i<connexion->client_size; i++) {
      int j;
      for(j=0;j<4096; j++){
	int data = connexion->server_rank;
	/* pack/send data */
	conn_server_send_data(connexion, &data, sizeof(int), 0, i);
      }
    }


    printf("[Server %d] I'm done\n", connexion->server_rank);
    sleep(1);

    printf("[Server %d] fini!\n", connexion->server_rank);

    conn_serv_unpublish_port(connexion);
    printf("[Server %d] finalize\n", connexion->server_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    conn_finalize_();
    MPI_Finalize();
    return 0;
}
