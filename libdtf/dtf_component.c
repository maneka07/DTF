/*
 *  Copyright (C) 2018, RIKEN
 *  See COPYRIGHT notice in top-level directory.
 * 
 * Written by 
 */
 
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
	double t_con;

    mode = gl_proc.comps[comp_id].connect_mode;
    if(mode == DTF_UNDEFINED) return 0;

	t_con = MPI_Wtime();
	
	if(gl_proc.conf.single_mpirun_mode){
		DTF_DBG(VERBOSE_DBG_LEVEL, "Duplicate intercomm");
		err = MPI_Pcontrol(0);
		assert(err == MPI_SUCCESS);
		err = MPI_Comm_dup(MPI_COMM_WORLD, &(gl_proc.comps[comp_id].comm));
		CHECK_MPI(err);
		MPI_Pcontrol(1);
		goto fn_exit;
	}  

    portname[0] = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

     if(mode == CONNECT_MODE_CLIENT){
        strcpy(service_name, gl_proc.comps[comp_id].name);
        strcat(service_name, "_");
        strcat(service_name, gl_proc.comps[gl_proc.my_comp].name);

     } else if(mode == CONNECT_MODE_SERVER){
        strcpy(service_name, gl_proc.comps[gl_proc.my_comp].name);
        strcat(service_name, "_");
        strcat(service_name, gl_proc.comps[comp_id].name);

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
        DTF_DBG(VERBOSE_DBG_LEVEL,   "Trying to connect to comp %s. Open portfile %s", gl_proc.comps[comp_id].name, portfile_name);

        /*To avoid contention on the PFS only rank 0 will try to get the port to connect to.*/
        if(myrank == 0){
            double t_st = MPI_Wtime();

            //try to open the file
            sleep(1);
            while((portfile = fopen(portfile_name, "rt")) == NULL){
                sleep(1);

                if(MPI_Wtime() - t_st >= gl_proc.conf.timeout){
                    DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: timed out waiting for port file %s.", portfile_name);
                    return 1;
                }
            }

			t_st = MPI_Wtime();
            
            while(1){
                fgets(portname, MPI_MAX_PORT_NAME, portfile);
                l1 = strlen(portname);
                if(portname[l1 - 1] == '\n')  //have read complete port
                    break;

                sleep(1);
              
                if(MPI_Wtime() - t_st >= gl_proc.conf.timeout){
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
        DTF_DBG(VERBOSE_DBG_LEVEL,   "%s will connect to service %s.", gl_proc.comps[gl_proc.my_comp].name, service_name);

        err = MPI_Comm_connect(portname, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &(gl_proc.comps[comp_id].comm));
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

        DTF_DBG(VERBOSE_DBG_LEVEL,   "%s starts listening for service %s.", gl_proc.comps[gl_proc.my_comp].name, service_name);

        err = MPI_Comm_accept(portname, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &(gl_proc.comps[comp_id].comm) );
        CHECK_MPI(err);
        DTF_DBG(VERBOSE_DBG_LEVEL,   "%s accepted connection on service %s.", gl_proc.comps[gl_proc.my_comp].name, service_name);

        err = MPI_Close_port(portname);
        CHECK_MPI(err);

        
		if(myrank == 0)
			unlink(portfile_name);
        
    }
    MPI_Comm_set_errhandler(gl_proc.comps[comp_id].comm, MPI_ERRORS_RETURN);
	//The two components are roughly synched now. Reset
	//the start time value
fn_exit:
	DTF_DBG(VERBOSE_DBG_LEVEL, "Took %.3f secs to esablish intercomm", MPI_Wtime() - t_con);
	gl_proc.walltime = MPI_Wtime();
    DTF_DBG(VERBOSE_DBG_LEVEL, "Intercomm established");
    
    return 0;
}

static void destroy_intercomm(int comp_id){

    int mode, err;
    char* global_path;
    char portfile_name[MAX_FILE_NAME];

    mode = gl_proc.comps[comp_id].connect_mode;

    if(mode == DTF_UNDEFINED) return;
    if(gl_proc.comps[comp_id].comm == MPI_COMM_NULL) return;
	
	if(gl_proc.conf.single_mpirun_mode){

		err = MPI_Comm_free(&(gl_proc.comps[comp_id].comm));
		CHECK_MPI(err);
		return;
	}  
    
    MPI_Comm_disconnect(&(gl_proc.comps[comp_id].comm));
    //rank 0 of the server component will remove the file
    if(gl_proc.myrank == 0 && mode == CONNECT_MODE_SERVER){

        global_path = getenv("DTF_GLOBAL_PATH");
        strcpy(portfile_name, global_path);
        strcat(portfile_name, "/port_");


         if(mode == CONNECT_MODE_CLIENT){
            strcat(portfile_name, gl_proc.comps[comp_id].name);
            strcat(portfile_name, "_");
            strcat(portfile_name, gl_proc.comps[gl_proc.my_comp].name);

         } else if(mode == CONNECT_MODE_SERVER){
            strcat(portfile_name, gl_proc.comps[gl_proc.my_comp].name);
            strcat(portfile_name, "_");
            strcat(portfile_name, gl_proc.comps[comp_id].name);

         }

        DTF_DBG(VERBOSE_DBG_LEVEL,   "Removing port file %s.", portfile_name);
        remove(portfile_name);
    }


}

int init_comp_comm(){

    int i, err;
    char *s;

    if(gl_proc.comps == NULL || gl_proc.ncomps == 1)
        return 0;

    s = getenv("DTF_GLOBAL_PATH");
    if(s == NULL){
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: please set DTF_GLOBAL_PATH.");
        goto panic_exit;
    }

    for(i = 0; i<gl_proc.ncomps; i++){
        if(i == gl_proc.my_comp){
            err = MPI_Comm_dup(MPI_COMM_WORLD, &(gl_proc.comps[i].comm));
            CHECK_MPI(err);
        } else {
			err = create_intercomm(gl_proc.comps[i].id, s);
			
			if(err) goto panic_exit;
		}
    }
    return 0;

panic_exit:
    return 1;
}

void finalize_comp_comm(){

    int i;
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finalizing communicators");
    for(i = 0; i<gl_proc.ncomps; i++){
		if(gl_proc.comps[i].in_msg_q != NULL){
			DTF_DBG(VERBOSE_DBG_LEVEL, "Recv msg queue for comp %s not empty. Will discard", gl_proc.comps[i].name);
			dtf_msg_t *msg = gl_proc.comps[i].in_msg_q, *tmp;
			while(msg != NULL){
				//DTF_DBG(VERBOSE_DBG_LEVEL, "%p", (void*)msg);
				DTF_DBG(VERBOSE_DBG_LEVEL, "tag %d",msg->tag);
				tmp = msg->next;
				DEQUEUE_ITEM(msg, gl_proc.comps[i].in_msg_q);
				delete_dtf_msg(msg);
				msg = tmp;
			}
		}
	   
       destroy_intercomm(gl_proc.comps[i].id);
    }

    /*Intra-comp communicator will be destroyed in the dtf_finalize() function*/
}

