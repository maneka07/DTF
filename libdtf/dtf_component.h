
#ifndef DTF_COMPONENT_H_INCLUDED
#define DTF_COMPONENT_H_INCLUDED

#define MAX_COMP_NAME 32

#include <mpi.h>

typedef struct component{
    unsigned int    id;
    char            name[MAX_COMP_NAME];
    int             connect_mode; /*0 - I am server, 1 - I am client, -1 - undefined (no interconnection)*/
    struct dtf_msg  *in_msg_q;   /*Queue of incoming messages (those that cannot be processed right now)*/
    struct dtf_msg  *out_msg_q;
    MPI_Comm        comm;   /*intra or inter component communicator*/
    int             finalized;   /*set to one when component calls finalize function*/
}component_t;

int 	init_comp_comm();
void 	finalize_comp_comm();

#endif // DTF_COMPONENT_H_INCLUDED
