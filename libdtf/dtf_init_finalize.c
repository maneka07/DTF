#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "dtf_init_finalize.h"
#include "dtf_common.h"
#include "dtf_file_buffer.h"
#include "dtf_util.h"
#include "dtf_req_match.h"
#include "dtf.h"


struct file_buffer * gl_filebuf_list = NULL;
struct component *gl_comps = NULL;
int gl_my_comp_id = -1;
int gl_ncomp = 0;
int gl_my_rank;

/*delete blancs in the beginning and end*/
static char* strip_str(char* str){

    static char s[ASCIILINESZ+1];
    int j;
    int len;

    if(str == NULL) return NULL;

    len = strlen(str);
    memset(s, 0, ASCIILINESZ+1);

    /* Get rid of spaces at the beginning of line */
    for(j = 0; j < len; j++){
        if(!isspace(str[j]))
            break;
    }
    if(len - j - 1 == 0){  //empty string
        s[0] = '\n';
        return (char*)s;
    }
    strcpy(s, str+j);
    len = strlen(s);
   // fprintf(stdout, "len %d\n", len);

    while ((len>0) &&
            ((s[len-1]=='\n') || (isspace(s[len-1])))) {
        s[len-1]=0 ;
        len-- ;
    }

    return (char*)s;
}

static char* str_to_lower(char* str){
    static char s[ASCIILINESZ+1];
    int i;

    if(str == NULL) return NULL;
    memset(s, 0, ASCIILINESZ+1);
    i=0 ;
    while (str[i] && i<ASCIILINESZ) {
        s[i] = (char)tolower((int)str[i]);
        i++ ;
    }

    return (char*)s;
}

static void print_config(){

    int i;
    file_buffer_t *fb = gl_filebuf_list;

    if(gl_comps == NULL) return;

    DTF_DBG(VERBOSE_DBG_LEVEL,   "Number of components: %d", gl_ncomp);
    DTF_DBG(VERBOSE_DBG_LEVEL,   "I am comp %d, name %s", gl_my_comp_id, gl_comps[gl_my_comp_id].name);
    for(i = 0; i<gl_ncomp; i++){
        if(i != gl_my_comp_id){
            if(gl_comps[i].connect_mode == CONNECT_MODE_CLIENT)
                DTF_DBG(VERBOSE_DBG_LEVEL,   "Component %s is my server", gl_comps[i].name);
            else if (gl_comps[i].connect_mode == CONNECT_MODE_SERVER)
                DTF_DBG(VERBOSE_DBG_LEVEL,   "Component %s is my client", gl_comps[i].name);
            else
                DTF_DBG(VERBOSE_DBG_LEVEL,   "No I/O with %s\n", gl_comps[i].name);
        }
    }

    while(fb != NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL,   "File %s, version %d, writer %s, reader %s ",
                fb->file_path, fb->version, gl_comps[fb->writer_id].name, gl_comps[fb->reader_id].name);
        switch(fb->iomode){
            case DTF_IO_MODE_UNDEFINED:
                DTF_DBG(VERBOSE_DBG_LEVEL,   "I/O mode: undefined ");
                break;
            case DTF_IO_MODE_FILE:
                DTF_DBG(VERBOSE_DBG_LEVEL,   "I/O mode: file i/o ");
                break;
            case DTF_IO_MODE_MEMORY:
                DTF_DBG(VERBOSE_DBG_LEVEL,   "I/O mode: direct transfer ");
                break;
            default:
                DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: unknown file mode %d", fb->iomode);
        }

        fprintf(stdout, "\n");
        fb=fb->next;
    }

}

static int check_config()
{

    file_buffer_t *fbuf = gl_filebuf_list;

    while(fbuf!= NULL){

        if(fbuf->iomode == DTF_IO_MODE_UNDEFINED){
            DTF_DBG(VERBOSE_ERROR_LEVEL,   " DTF Error: I/O mode for file %s underfined", fbuf->file_path);
            return 1;
        }

        fbuf = fbuf->next;
    }

    return 0;
}

static int create_intercomm(int comp_id, char* global_path){

    char portname[MPI_MAX_PORT_NAME];
    char service_name[MAX_COMP_NAME*2+1];
    FILE* portfile;
    int errno, mode, myrank, l1, l2;
    char portfile_name[MAX_FILE_NAME];
    int wait_timeout;
    int nranks1, nranks2;

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

        errno = MPI_Bcast(portname, MPI_MAX_PORT_NAME, MPI_CHAR, 0, MPI_COMM_WORLD );
        CHECK_MPI(errno);
        DTF_DBG(VERBOSE_DBG_LEVEL,   "%s will connect to service %s.", gl_comps[gl_my_comp_id].name, service_name);

        errno = MPI_Comm_connect(portname, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &(gl_comps[comp_id].comm));
        CHECK_MPI(errno);
    } else if(mode == CONNECT_MODE_SERVER){
        DTF_DBG(VERBOSE_DBG_LEVEL,   "Creating connection for %s.", service_name);

        errno = MPI_Open_port(MPI_INFO_NULL, portname);
        CHECK_MPI(errno);

        if(myrank == 0){
           DTF_DBG(VERBOSE_DBG_LEVEL, "Write port to file %s", portfile_name);
           portfile = fopen(portfile_name, "wt");
           fprintf(portfile, "%s\n", portname);
           fclose(portfile);
        }

        DTF_DBG(VERBOSE_DBG_LEVEL,   "%s starts listening for service %s.", gl_comps[gl_my_comp_id].name, service_name);

        errno = MPI_Comm_accept(portname, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &(gl_comps[comp_id].comm) );
        CHECK_MPI(errno);
        DTF_DBG(VERBOSE_DBG_LEVEL,   "%s accepted connection on service %s.", gl_comps[gl_my_comp_id].name, service_name);

        errno = MPI_Close_port(portname);
        CHECK_MPI(errno);
    }
    MPI_Comm_set_errhandler(gl_comps[comp_id].comm, MPI_ERRORS_RETURN);
    /*Check that the number of ranks is the same in both communicators*/
    MPI_Comm_size(MPI_COMM_WORLD, &nranks1);
    MPI_Comm_remote_size(gl_comps[comp_id].comm, &nranks2);
    if(nranks1 != nranks2){
        DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF can only work if the number of processes in all components is the same. Aborting.");
        return 1;
    }
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


void clean_config(){

    if(gl_comps != NULL)
        dtf_free(gl_comps, gl_ncomp*sizeof(component_t));
}

int load_config(const char *ini_name, const char *comp_name){

  FILE*      in;
  char       line[ASCIILINESZ+1];
  char       param[ASCIILINESZ], value[ASCIILINESZ];
  //char       comp_name[MAX_COMP_NAME];
  int        i, len, lineno=0;
  struct file_buffer* cur_fb = NULL;
 // ini = iniparser_load(ini_name);

  line[0] = 0;
  param[0] = 0;
  value[0] = 0;

  DTF_DBG(VERBOSE_DBG_LEVEL,   "Load config %s for comp %s", ini_name, comp_name);

  if ((in=fopen(ini_name, "r"))==NULL) {
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "Eror: cannot open %s\n", ini_name);
        return 1 ;
  }

  gl_conf.buffered_req_match = 0;

  while(!feof(in)){

    fgets(line, ASCIILINESZ, in);
    lineno++;

    len = strlen(line);
    if(strlen(line) == 0 || line[0] == '#' || line[0] == ';' || line[0] == '!') //comment or empty line
        continue;

    if(line[len-1] != '\n' && !feof(in)){
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: Line %d in the ini file too long.", lineno);
        goto panic_exit;
    }

    strcpy(line, strip_str(line));
    strcpy(line, str_to_lower(line));

    if(strcmp(line, "[file]")==0){

        //check that we have already parsed the [info] section
        if (gl_ncomp == 0){
           DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: [file] section defined before [info] section.");
            goto panic_exit;
        } else{
            for(i = 0; i < gl_ncomp; i++)
                if(gl_comps[i].name[0] == '\0'){
                    DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: not all component names defined.");
                    goto panic_exit;
                }
        }

        if (cur_fb != NULL){
            /* first check that all required params for the
            curent buffer_file have been defined */
            if(cur_fb->file_path[0] == '\0'){
                DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: filename not set.");
                goto panic_exit;
            }

            if( (cur_fb->iomode != DTF_IO_MODE_FILE) && (cur_fb->writer_id == -1 || cur_fb->reader_id == -1) ){
                DTF_DBG(VERBOSE_ERROR_LEVEL,"DTF Error: file reader or writer not set for file %s.", cur_fb->file_path);
                goto panic_exit;
            }
        }
        cur_fb = new_file_buffer();
        assert(cur_fb != NULL);

        add_file_buffer(&gl_filebuf_list, cur_fb);

    } else if (sscanf (line, "%[^=] = \"%[^\"]\"", param, value) == 2
           ||  sscanf (line, "%[^=] = '%[^\']'",   param, value) == 2
           ||  sscanf (line, "%[^=] = %[^;#]",     param, value) == 2) {

        /* Usual key=value */
        /*
         * sscanf cannot handle '' or "" as empty values
         * this is done here
         */
        if (!strcmp(value, "\"\"") || (!strcmp(value, "''"))) {
           continue;
        }

        //if(gl_verbose >= VERBOSE_DBG_LEVEL) fprintf(stdout, "param=%s, val=%s\n", param, value);

        if(strcmp(param, "ncomp") == 0){
            gl_ncomp = atoi(value);
            if(gl_ncomp <= 0){
                DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: invalid number of components.");
                goto panic_exit;
            }
            gl_comps = (struct component*)dtf_malloc(gl_ncomp*sizeof(struct component));
            assert(gl_comps!=NULL);
            for(i = 0; i<gl_ncomp; i++){
                gl_comps[i].id = i;
                gl_comps[i].connect_mode = DTF_UNDEFINED;
                gl_comps[i].name[0] = 0;
                gl_comps[i].comm = MPI_COMM_NULL;
            }

        } else if(strcmp(param, "comp_name") == 0){
            assert(gl_ncomp!=0);
            assert(strlen(value) <= MAX_COMP_NAME);
            for(i = 0; i<gl_ncomp; i++){
                if(gl_comps[i].name[0] == 0){
                    strcpy(gl_comps[i].name, value);
                    if(strcmp(gl_comps[i].name, comp_name) == 0)
                        gl_my_comp_id = i;
                    break;
                }
            }
//        } else if(strcmp(param, "distr_mode") == 0){
//            if(strcmp(value, "static") == 0)
//                gl_conf.distr_mode = DISTR_MODE_STATIC;
//            else if(strcmp(value, "req_match") == 0)
//                gl_conf.distr_mode = DISTR_MODE_REQ_MATCH;
//            else{
//                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: 3 unknown data distribution mode");
//                goto panic_exit;
//            }
        } else if(strcmp(param, "buffered_req_match") == 0){

            gl_conf.buffered_req_match = atoi(value);
            if(gl_conf.buffered_req_match == 0)
                DTF_DBG(VERBOSE_DBG_LEVEL, "DTF: buffered_req_match disabled");
            else if(gl_conf.buffered_req_match == 1)
                DTF_DBG(VERBOSE_DBG_LEVEL, "DTF: buffered_req_match enabled");
            else {
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: Value for buffered_req_match should be 0 or 1.");
                goto panic_exit;
            }
        } else if(strcmp(param, "explicit_match") == 0){
            assert(cur_fb != NULL);
            cur_fb->explicit_match = atoi(value);

            if(cur_fb->explicit_match == 0)
                DTF_DBG(VERBOSE_ALL_LEVEL, "DTF: explicit match disabled for %s", cur_fb->file_path);
            else if(cur_fb->explicit_match == 1)
                DTF_DBG(VERBOSE_ALL_LEVEL, "DTF: explicit match enabled for %s", cur_fb->file_path);
            else {
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: Value for explicit match should be 0 or 1 for %s.", cur_fb->file_path);
                goto panic_exit;
            }

        } else if(strcmp(param, "filename") == 0){
            assert(cur_fb != NULL);
            assert(strlen(value) <= MAX_FILE_NAME);
            strcpy(cur_fb->file_path, value);

        } else if(strcmp(param, "alias") == 0){
            assert(cur_fb != NULL);
            assert(strlen(value) <= MAX_FILE_NAME);
            strcpy(cur_fb->alias_name, value);

        } else if(strcmp(param, "version") == 0){
            assert(cur_fb != NULL);
            cur_fb->version = atoi(value);

        } else if(strcmp(param, "writer") == 0){
            assert(cur_fb != NULL);
            assert(gl_ncomp != 0);
            for(i = 0; i < gl_ncomp; i++){
                if(strcmp(value, gl_comps[i].name) == 0){
                    cur_fb->writer_id = i;
                    break;
                }
            }

        } else if(strcmp(param, "reader") == 0){
            assert(cur_fb != NULL);
            assert(gl_ncomp != 0);

            if(cur_fb->reader_id != -1){
                DTF_DBG(VERBOSE_ERROR_LEVEL, "File %s cannot have multiple readers", cur_fb->file_path);
                goto panic_exit;
            }
            for(i = 0; i < gl_ncomp; i++){
                if(strcmp(value, gl_comps[i].name) == 0){
                    cur_fb->reader_id = i;
                    break;
                }
            }

        } else if(strcmp(param, "mode") == 0){
            assert(cur_fb != NULL);

            if(strcmp(value, "file") == 0)
                cur_fb->iomode = DTF_IO_MODE_FILE;
            else if(strcmp(value, "memory") == 0){

                if(gl_conf.distr_mode == DTF_UNDEFINED){
                    DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: data distribution mode is not defined.");
                    goto panic_exit;
                }

                cur_fb->iomode = DTF_IO_MODE_MEMORY;
            }

        } else {
            DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: unknown parameter %s.", param);
            goto panic_exit;
        }

    } else if (sscanf(line, "%[^=] = %[;#]", param, value)==2
           ||  sscanf(line, "%[^=] %[=]", param, value) == 2) {
        /*
         * Special cases:
         * key=
         * key=;
         * key=#
         */
        //value[0]=0 ;
        continue;
    }
    }

//    for(i = 0; i < gl_ncomp; i++)
//        if(strcmp(gl_comps[i].name, comp_name) == 0)
//            gl_my_comp_id = i;



/*
  During the parsing of the configuarion file, if a file direct data transfer will
  take place between two components, each component is automatically set as a server
  (CONNECT_MODE_SERVER) or client(CONNECT_MODE_CLIENT), respectively, in order to
  esablish an intercommunicatior between these two components later. If there is no
  file flow then the connection mode is not set(DTF_UNDEFINED).
*/
    cur_fb = gl_filebuf_list;
    while(cur_fb != NULL){
        if(cur_fb->writer_id == gl_my_comp_id){
            if(gl_comps[ cur_fb->reader_id ].connect_mode == DTF_UNDEFINED )
               gl_comps[ cur_fb->reader_id ].connect_mode = CONNECT_MODE_SERVER; // I am server

        } else if( gl_comps[ cur_fb->writer_id ].connect_mode == DTF_UNDEFINED){
            if(cur_fb->reader_id == gl_my_comp_id )
               gl_comps[ cur_fb->writer_id ].connect_mode = CONNECT_MODE_CLIENT;  //I am client
        }

        cur_fb = cur_fb->next;
    }

    if(gl_verbose)print_config();
    fclose(in);

    if (check_config())
        goto panic_exit;

    return 0;

panic_exit:
    clean_config();
    fclose(in);
    return 1;
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
            MPI_Comm_dup(MPI_COMM_WORLD, &(gl_comps[i].comm));
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
       destroy_intercomm(gl_comps[i].id);
    }
}

void finalize_files()
{
    int file_cnt = 0, compl_cnt = 0;
    file_buffer_t *fbuf = gl_filebuf_list;
    /*First, writers compute number of files
      they wrote*/
    while(fbuf != NULL){

        DTF_DBG(VERBOSE_DBG_LEVEL, "File %s, close_flag %d, notif_flag %d", fbuf->file_path, fbuf->rdr_closed_flag, fbuf->fready_notify_flag);
        if(fbuf->fready_notify_flag != DTF_UNDEFINED && fbuf->root_writer == gl_my_rank){
            assert(fbuf->writer_id == gl_my_comp_id);
            file_cnt++;
        }
        fbuf = fbuf->next;
    }

    DTF_DBG(VERBOSE_DBG_LEVEL, "Process has to finalize notifications for %d files", file_cnt);
    while(file_cnt != compl_cnt){
        compl_cnt = 0;
        fbuf = gl_filebuf_list;
        while(fbuf != NULL){
            if(fbuf->fready_notify_flag != DTF_UNDEFINED && fbuf->root_writer == gl_my_rank){
                assert(fbuf->writer_id == gl_my_comp_id);

                if(fbuf->iomode == DTF_IO_MODE_FILE){
                        if(fbuf->fready_notify_flag == RDR_NOTIFIED)
                            compl_cnt++;
                        else{
                            assert(fbuf->fready_notify_flag == RDR_NOT_NOTIFIED);
                            while(fbuf->root_reader == -1)
                                progress_io_matching();
                            notify_file_ready(fbuf);
                        }
                } else if(fbuf->iomode == DTF_IO_MODE_MEMORY && gl_conf.distr_mode == DISTR_MODE_REQ_MATCH){
                        if(fbuf->rdr_closed_flag)
                            compl_cnt++;
                        else
                            progress_io_matching();
                }
            }

            fbuf = fbuf->next;
        }
    }
    MPI_Barrier(gl_comps[gl_my_comp_id].comm);
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished finalizing notifications. Will delete file bufs");
    /*Now, delete all file buffers*/
    fbuf = gl_filebuf_list;
    while(fbuf != NULL){
        gl_filebuf_list = fbuf->next;
        delete_file_buffer(&gl_filebuf_list, fbuf);
        fbuf = gl_filebuf_list;
    }
}

