#include <assert.h>
#include "farb_connect.h"

int main( int argc, char **argv )
{
    char service_name[] = "ocean";
    MPI_Comm server;
    struct client_connection*connexion = NULL;
    MPI_Init(&argc, &argv);
    sleep(2);

    conn_init_("sample.ini", "");
    connexion = conn_client_connect(service_name, &server);

    /* receive the data */
    int i;
    int k =-1;
    for(i=0; i<connexion->server_size; i++) {
      printf("recv(%d/%d)\n", i, connexion->server_size);
      conn_client_recv_data(connexion, &k, sizeof(int), connexion->server_size, i);
      printf("Received: %d\n", k);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    /* process the data */

    /* recv 4096 * int from each server */
    for(i=0; i<connexion->server_size; i++) {
      int j;
      printf("Receiving 4096 messages from %d\n", i);
      for(j=0; j<4096; j++) {
	conn_client_recv_data(connexion, &k, sizeof(int), connexion->server_size, i);
      }
    }

    printf("[Client %d] fini!\n", connexion->client_rank);
    conn_client_disconnect(connexion);

    printf("[Client %d] finalize\n", connexion->client_rank);
    conn_finalize_();
    MPI_Finalize();

    return 0;
}
