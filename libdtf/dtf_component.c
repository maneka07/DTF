#include <unistd.h>

#include "dtf_component.h"
#include "dtf_util.h"
#include "dtf.h"

static int create_intercomm(int comp_id, char* global_path){

    char portname[MPI_MAX_PORT_NAME];
    char service_name[MAX_COMP_NAME*2+1];
    FILE* portfile;
    int err, mode, myrank, l1, l2;
    char portfile_name[MAX_FILE_NAME];
    int wait_timeout;

    mode = gl_comps[comp_id].connect_mode;
    if(mode == DTF_UNDEFINED) return 0;

    portname[0] = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

     if(mode == CONNECT_MODE_CLIENT){
        strcpy(service_name, gl_comps[comp_id].name);
        strcat(service_name, "_");
        strcat(service_name, gl_comps[gl_my_comp_id].name);

     } else if(mode == CONNECT_MODE_SERVER){
        strcpy(service_name, gl_comps[gl_my_comp_id].name);
        strcat(service_name, "_");
        strcat(service_name, gl_comps[comp_id].name);

     } else {
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: connect mode unknown.");
        return 1;
     }

     l1 = strlen(global_path);
     l2 = strlen(service_name);

     if(l1+l2 > MAX_FILE_NAME){
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: Global file path or component name too long.");
        return 1;
     }

     strcpy(portfile_name, global_path);
     strcat(portfile_name, "/port_");
     strcat(portfile_name, service_name);


    if(mode == CONNECT_MODE_CLIENT){
        DTF_DBG(VERBOSE_DBG_LEVEL,   "Trying to connect to comp %s. Open portfile %s", gl_comps[comp_id].name, portfile_name);

        /*To avoid contention on the PFS only rank 0 will try to get the port to connect to.*/
        if(myrank == 0){
            wait_timeout = 0;

            //try to open the file
            sleep(1);
            while((portfile = fopen(portfile_name, "rt")) == NULL){
                sleep(1);
                wait_timeout++;

                if(wait_timeout == DTF_TIMEOUT){
                    DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: timed out waiting for port file %s.", portfile_name);
                    return 1;
                }
            }

            //try to read the file
            wait_timeout = 0;

            while(1){
                fgets(portname, MPI_MAX_PORT_NAME, portfile);
                l1 = strlen(portname);
                if(portname[l1 - 1] == '\n')  //have read complete port
                    break;

                sleep(1);
                wait_timeout++;

                if(wait_timeout == DTF_TIMEOUT){
                    DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: timed out waiting for port name in %s.", portfile_name);
                    return 1;
                }
            }
            fclose(portfile);

            if(portname[0] == '\n'){
                DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: no port specified in the portfile.");
                return 1;
            }
        }

        err = MPI_Bcast(portname, MPI_MAX_PORT_NAME, MPI_CHAR, 0, MPI_COMM_WORLD );
        CHECK_MPI(err);
        DTF_DBG(VERBOSE_DBG_LEVEL,   "%s will connect to service %s.", gl_comps[gl_my_comp_id].name, service_name);

        err = MPI_Comm_connect(portname, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &(gl_comps[comp_id].comm));
        CHECK_MPI(err);
    } else if(mode == CONNECT_MODE_SERVER){
        DTF_DBG(VERBOSE_DBG_LEVEL,   "Creating connection for %s.", service_name);

        err = MPI_Open_port(MPI_INFO_NULL, portname);
        CHECK_MPI(err);

        if(myrank == 0){
           DTF_DBG(VERBOSE_DBG_LEVEL, "Write port to file %s", portfile_name);
           portfile = fopen(portfile_name, "wt");
           fprintf(portfile, "%s\n", portname);
           fclose(portfile);
        }

        DTF_DBG(VERBOSE_DBG_LEVEL,   "%s starts listening for service %s.", gl_comps[gl_my_comp_id].name, service_name);

        err = MPI_Comm_accept(portname, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &(gl_comps[comp_id].comm) );
        CHECK_MPI(err);
        DTF_DBG(VERBOSE_DBG_LEVEL,   "%s accepted connection on service %s.", gl_comps[gl_my_comp_id].name, service_name);

        err = MPI_Close_port(portname);
        CHECK_MPI(err);
        
    }
    MPI_Comm_set_errhandler(gl_comps[comp_id].comm, MPI_ERRORS_RETURN);

	//The two components are roughly synched now. Reset
	//the start time value
	gl_stats.walltime = MPI_Wtime();

    return 0;
}

static void destroy_intercomm(int comp_id){

    int mode;
    char* global_path;
    char portfile_name[MAX_FILE_NAME];

    mode = gl_comps[comp_id].connect_mode;

    if(mode == DTF_UNDEFINED) return;
    if(gl_comps[comp_id].comm == MPI_COMM_NULL) return;
    MPI_Comm_disconnect(&(gl_comps[comp_id].comm));
    //rank 0 of the server component will remove the file
    if(gl_my_rank == 0 && mode == CONNECT_MODE_SERVER){

        global_path = getenv("DTF_GLOBAL_PATH");
        assert(global_path!=NULL);
        strcpy(portfile_name, global_path);
        strcat(portfile_name, "/port_");


         if(mode == CONNECT_MODE_CLIENT){
            strcat(portfile_name, gl_comps[comp_id].name);
            strcat(portfile_name, "_");
            strcat(portfile_name, gl_comps[gl_my_comp_id].name);

         } else if(mode == CONNECT_MODE_SERVER){
            strcat(portfile_name, gl_comps[gl_my_comp_id].name);
            strcat(portfile_name, "_");
            strcat(portfile_name, gl_comps[comp_id].name);

         }

        DTF_DBG(VERBOSE_DBG_LEVEL,   "Removing port file %s.", portfile_name);
        remove(portfile_name);
    }


}

int init_comp_comm(){

    int i, err;
    char *s;

    if(gl_comps == NULL || gl_ncomp == 1)
        return 0;

    s = getenv("DTF_GLOBAL_PATH");
    if(s == NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: please set DTF_GLOBAL_PATH.");
        goto panic_exit;
    }

    for(i = 0; i<gl_ncomp; i++){
        if(i == gl_my_comp_id){
            err = MPI_Comm_dup(MPI_COMM_WORLD, &(gl_comps[i].comm));
            CHECK_MPI(err);
            continue;
        }
        err = create_intercomm(gl_comps[i].id, s);
        if(err)
            goto panic_exit;
    }
    return 0;

panic_exit:
    return 1;
}

void finalize_comp_comm(){

    int i;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finalizing communicators");
    for(i = 0; i<gl_ncomp; i++){
		if(gl_comps[i].in_msg_q != NULL){
			DTF_DBG(VERBOSE_DBG_LEVEL, "Recv msg queue for comp %s not empty. Will discard", gl_comps[i].name);
			dtf_msg_t *msg = gl_comps[i].in_msg_q, *tmp;
			while(msg != NULL){
				//DTF_DBG(VERBOSE_DBG_LEVEL, "%p", (void*)msg);
				DTF_DBG(VERBOSE_DBG_LEVEL, "tag %d",msg->tag);
				tmp = msg->next;
				DEQUEUE_ITEM(msg, gl_comps[i].in_msg_q);
				delete_dtf_msg(msg);
				msg = tmp;
			}
		}
		
       destroy_intercomm(gl_comps[i].id);
    }

    /*Intra-comp communicator will be destroyed in the dtf_finalize() function*/
}

