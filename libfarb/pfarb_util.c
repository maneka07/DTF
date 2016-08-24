/*
 * Copyright (C) 2015, Advanced Institute for Computational Science, RIKEN
 * Author: Jianwei Liao(liaotoad@gmail.com)
 */

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>
#include "pfarb_util.h"
#include "pfarb.h"
#include "pfarb_mem.h"

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

void clean_config(){
    file_buffer_t *tmp;

    if(gl_comps != NULL)
        free(gl_comps);

    while(gl_filebuf_list != NULL){
        tmp = gl_filebuf_list;
        delete_file_buffer(&gl_filebuf_list, tmp);
    }
}

static void print_config(){

    int i;
    file_buffer_t *fb = gl_filebuf_list;

    if(gl_comps == NULL) return;

    FARB_DBG(VERBOSE_DBG_LEVEL,   "Number of components: %d", gl_ncomp);
    FARB_DBG(VERBOSE_DBG_LEVEL,   "I am comp %d, name %s", gl_my_comp_id, gl_comps[gl_my_comp_id].name);
    for(i = 0; i<gl_ncomp; i++){
        if(i != gl_my_comp_id){
            if(gl_comps[i].connect_mode == CONNECT_MODE_CLIENT)
                FARB_DBG(VERBOSE_DBG_LEVEL,   "Component %s is my server", gl_comps[i].name);
            else if (gl_comps[i].connect_mode == CONNECT_MODE_SERVER)
                FARB_DBG(VERBOSE_DBG_LEVEL,   "Component %s is my client", gl_comps[i].name);
            else
                FARB_DBG(VERBOSE_DBG_LEVEL,   "No I/O with %s\n", gl_comps[i].name);
        }
    }

    while(fb != NULL){
        FARB_DBG(VERBOSE_DBG_LEVEL,   "File %s, version %d, writer %s, readers: ", fb->file_path, fb->version, gl_comps[fb->writer_id].name);
        for(i = 0; i < fb->nreaders; i++)
            FARB_DBG(VERBOSE_DBG_LEVEL,   "%s, ", gl_comps[(fb->reader_ids)[i]].name);
        switch(fb->mode){
            case FARB_IO_MODE_UNDEFINED:
                FARB_DBG(VERBOSE_DBG_LEVEL,   "I/O mode: undefined ");
                break;
            case FARB_IO_MODE_FILE:
                FARB_DBG(VERBOSE_DBG_LEVEL,   "I/O mode: file i/o ");
                break;
            case FARB_IO_MODE_MEMORY:
                FARB_DBG(VERBOSE_DBG_LEVEL,   "I/O mode: direct transfer ");
                break;
            default:
                FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: unknown file mode %d", fb->mode);
        }

        fprintf(stdout, "\n");
        fb=fb->next;
    }

}

static int check_config()
{

    file_buffer_t *fbuf = gl_filebuf_list;

    while(fbuf!= NULL){

        if(fbuf->mode == FARB_IO_MODE_UNDEFINED){
            FARB_DBG(VERBOSE_ERROR_LEVEL,   " FARB Error: I/O mode for file %s underfined", fbuf->file_path);
            return 1;
        }

        fbuf = fbuf->next;
    }

    return 0;
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


  FARB_DBG(VERBOSE_DBG_LEVEL,   "Load config %s for comp %s", ini_name, comp_name);


  if ((in=fopen(ini_name, "r"))==NULL) {
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "Eror: cannot open %s\n", ini_name);
        return 1 ;
  }

  while(!feof(in)){

    fgets(line, ASCIILINESZ, in);
    lineno++;

    len = strlen(line);
    if(strlen(line) == 0 || line[0] == '#' || line[0] == ';') //comment or empty line
        continue;

    if(line[len-1] != '\n' && !feof(in)){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: Line %d in the ini file too long.", lineno);
        goto panic_exit;
    }

    strcpy(line, strip_str(line));
    strcpy(line, str_to_lower(line));

    if(strcmp(line, "[file]")==0){

        //check that we have already parsed the [info] section
        if (gl_ncomp == 0){
           FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: [file] section defined before [info] section.");
            goto panic_exit;
        } else{
            for(i = 0; i < gl_ncomp; i++)
                if(gl_comps[i].name[0] == '\0'){
                    FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: not all component names defined.");
                    goto panic_exit;
                }
        }

        if (cur_fb != NULL){
        /* first check that all required params for the
        curent buffer_file have been defined */
        if(cur_fb->file_path[0] == '\0'){
            FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: filename not set.");
            goto panic_exit;
          } else if( cur_fb->writer_id == -1 || cur_fb->reader_ids == NULL){
            FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: file reader or writer not set for file %s.", cur_fb->file_path);
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
                FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: invalid number of components.");
                goto panic_exit;
            }
            gl_comps = (struct component*)malloc(gl_ncomp*sizeof(struct component));
            assert(gl_comps!=NULL);
            for(i = 0; i<gl_ncomp; i++){
                gl_comps[i].id = i;
                gl_comps[i].connect_mode = CONNECT_MODE_UNDEFINED;
                gl_comps[i].name[0] = 0;
                gl_comps[i].intercomm = MPI_COMM_NULL;
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
                    if(gl_my_comp_id == i)cur_fb->write_flag = 1;
                    break;
                }
            }

        } else if(strcmp(param, "reader") == 0){
            assert(cur_fb != NULL);
            assert(gl_ncomp != 0);

            if(cur_fb->reader_ids == NULL){
                cur_fb->reader_ids = malloc (sizeof(int));
                cur_fb->nreaders = 1;
                assert(cur_fb->reader_ids != NULL);
            } else {
                cur_fb->reader_ids = realloc(cur_fb->reader_ids, cur_fb->nreaders+1);
                assert(cur_fb->reader_ids != NULL);
                cur_fb->nreaders++;
            }
            for(i = 0; i < gl_ncomp; i++){
                if(strcmp(value, gl_comps[i].name) == 0){
                    cur_fb->reader_ids[cur_fb->nreaders-1] = i;
                    if(gl_my_comp_id == i) cur_fb->read_flag = 1;
                    break;
                }
            }

        } else if(strcmp(param, "mode") == 0){
            assert(cur_fb != NULL);

            if(strcmp(value, "file") == 0)
                cur_fb->mode = FARB_IO_MODE_FILE;
            else if(strcmp(value, "memory") == 0)
                cur_fb->mode = FARB_IO_MODE_MEMORY;

        } else {
            FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: unknown parameter %s.", param);
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
  file flow then the connection mode is not set(CONNECT_MODE_UNDEFINED).
*/
    cur_fb = gl_filebuf_list;
    while(cur_fb != NULL){

        if(cur_fb->writer_id == gl_my_comp_id){
            for(i = 0; i < cur_fb->nreaders; i++){
                if(gl_comps[ (cur_fb->reader_ids)[i] ].connect_mode == CONNECT_MODE_UNDEFINED )
                   gl_comps[ (cur_fb->reader_ids)[i] ].connect_mode = CONNECT_MODE_SERVER; // I am server
            }
        } else if( gl_comps[ cur_fb->writer_id ].connect_mode == CONNECT_MODE_UNDEFINED){
            for(i = 0; i < cur_fb->nreaders; i++){
                if(cur_fb->reader_ids[i] == gl_my_comp_id ){
                   gl_comps[ cur_fb->writer_id ].connect_mode = CONNECT_MODE_CLIENT;  //I am client
                }
            }
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


static int create_intercomm(int comp_id, char* global_path){

    char portname[MPI_MAX_PORT_NAME];
    char service_name[MAX_COMP_NAME*2+1];
    FILE* portfile;
    int errno, mode, myrank, l1, l2;
    char portfile_name[MAX_FILE_NAME];
    int wait_timeout;

    mode = gl_comps[comp_id].connect_mode;
    if(mode == CONNECT_MODE_UNDEFINED) return 0;

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
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: connect mode unknown.");
        goto panic_exit;
     }

     l1 = strlen(global_path);
     l2 = strlen(service_name);

     if(l1+l2 > MAX_FILE_NAME){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: Global file path or component name too long.");
        goto panic_exit;
     }

     strcpy(portfile_name, global_path);
     strcat(portfile_name, "/port_");
     strcat(portfile_name, service_name);


    if(mode == CONNECT_MODE_CLIENT){
        FARB_DBG(VERBOSE_DBG_LEVEL,   "Trying to connect to comp %s.", gl_comps[comp_id].name);

        /*To avoid contention on the PFS only rank 0 will try to get the port to connect to.*/
        if(myrank == 0){
            wait_timeout = 0;

            //try to open the file
            while((portfile = fopen(portfile_name, "rt")) == NULL){
                sleep(1);
                wait_timeout++;

                if(wait_timeout == FARB_TIMEOUT){
                    FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: timed out waiting for port file %s.", portfile_name);
                    goto panic_exit;
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

                if(wait_timeout == FARB_TIMEOUT){
                    FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: timed out waiting for port name in %s.", portfile_name);
                    goto panic_exit;
                }
            }
            fclose(portfile);

            if(portname[0] == '\n'){
                FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: no port specified in the portfile.");
                goto panic_exit;
            }
        }

        errno = MPI_Bcast(portname, MPI_MAX_PORT_NAME, MPI_CHAR, 0, MPI_COMM_WORLD );
        CHECK_MPI(errno);
        FARB_DBG(VERBOSE_DBG_LEVEL,   "%s will connect to service %s.", gl_comps[gl_my_comp_id].name, service_name);

        errno = MPI_Comm_connect(portname, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &(gl_comps[comp_id].intercomm));
        CHECK_MPI(errno);
    } else if(mode == CONNECT_MODE_SERVER){
        FARB_DBG(VERBOSE_DBG_LEVEL,   "Creating connection for %s.", service_name);

        errno = MPI_Open_port(MPI_INFO_NULL, portname);
        CHECK_MPI(errno);

        if(myrank == 0){

           portfile = fopen(portfile_name, "wt");
           fprintf(portfile, "%s\n", portname);
           fclose(portfile);
        }

        FARB_DBG(VERBOSE_DBG_LEVEL,   "%s starts listening for service %s.", gl_comps[gl_my_comp_id].name, service_name);

        errno = MPI_Comm_accept(portname, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &(gl_comps[comp_id].intercomm) );
        CHECK_MPI(errno);
        FARB_DBG(VERBOSE_DBG_LEVEL,   "%s accepted connection on service %s.", gl_comps[gl_my_comp_id].name, service_name);

        errno = MPI_Close_port(portname);
        CHECK_MPI(errno);
    }
    MPI_Errhandler_set(gl_comps[comp_id].intercomm, MPI_ERRORS_RETURN);
    return 0;
panic_exit:

    return 1;
}

int init_comp_comm(){

    int i, err;
    char *s;

    if(gl_comps == NULL || gl_ncomp == 1)
        return 0;

    s = getenv("FARB_GLOBAL_PATH");
    if(s == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: please set FARB_GLOBAL_PATH.");
        goto panic_exit;
    }

    for(i = 0; i<gl_ncomp; i++){

        err = create_intercomm(gl_comps[i].id, s);
        if(err)
            goto panic_exit;
    }

    return 0;

panic_exit:
    return 1;
}


static void destroy_intercomm(int comp_id){

    int mode;
    char* global_path;
    char portfile_name[MAX_FILE_NAME];

    mode = gl_comps[comp_id].connect_mode;

    if(mode == CONNECT_MODE_UNDEFINED) return;
    if(gl_comps[comp_id].intercomm == MPI_COMM_NULL) return;

    MPI_Comm_disconnect(&(gl_comps[comp_id].intercomm));

    //rank 0 of the server component will remove the file
    if(gl_my_rank == 0 && mode == CONNECT_MODE_SERVER){

        global_path = getenv("FARB_GLOBAL_PATH");
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

        FARB_DBG(VERBOSE_DBG_LEVEL,   "Removing port file %s.", portfile_name);
        remove(portfile_name);
    }


}

void finalize_comp_comm(){

    int i;

    for(i = 0; i<gl_ncomp; i++){
       destroy_intercomm(gl_comps[i].id);
    }
}

int get_write_flag(const char* filename)
{
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: file %s not in the configuration file", filename);
        assert(0);
    }
    return fbuf->write_flag;
}

int get_read_flag(const char* filename)
{
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: file %s not in the configuration file", filename);
        assert(0);
    }

    return fbuf->read_flag;
}

int file_buffer_ready(const char* filename)
{
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    if(fbuf == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: file %s not in the configuration file", filename);
        assert(0);
    }
    return fbuf->is_ready;
}

void progress_io()
{
    MPI_Status status;
    int i, j, flag, src, errno;
    char filename[MAX_FILE_NAME];
    file_buffer_t *fbuf;
    farb_var_t *var;
    buffer_node_t *node;
    void *buf;
    int buf_sz;
    MPI_Request sreq, *reqs=NULL;

    for(i = 0; i < gl_ncomp; i++){
        if( i == gl_my_comp_id || gl_comps[i].intercomm == MPI_COMM_NULL)
            continue;

        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_comps[i].intercomm, &flag, &status);
        if(!flag)
            continue;

        switch(status.MPI_TAG){
            case FILE_READY_TAG:
                src = status.MPI_SOURCE;
                errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, FILE_READY_TAG, gl_comps[i].intercomm, &status);
                CHECK_MPI(errno);
                FARB_DBG(VERBOSE_DBG_LEVEL,   "Receive FILE_READY notif for %s", filename);

                fbuf = find_file_buffer(gl_filebuf_list, filename);
                assert(fbuf != NULL);

                if(fbuf->mode == FARB_IO_MODE_MEMORY){
                    /*Notify I am ready to receive*/
                    errno = MPI_Send(filename, MAX_FILE_NAME, MPI_CHAR, src, RECV_READY_TAG, gl_comps[i].intercomm);
                    CHECK_MPI(errno);

                    /*Receive the header*/
                    MPI_Probe(src, HEADER_TAG, gl_comps[i].intercomm, &status);
                    MPI_Get_count(&status, MPI_BYTE, &buf_sz);
                    fbuf->hdr_sz = (MPI_Offset)buf_sz;
                    FARB_DBG(VERBOSE_DBG_LEVEL, "Hdr size to receive %d", buf_sz);
                    fbuf->header = malloc((size_t)buf_sz);
                    assert(fbuf->header != NULL);
                    errno = MPI_Recv(fbuf->header, buf_sz, MPI_BYTE, src, HEADER_TAG, gl_comps[i].intercomm, &status);
                    CHECK_MPI(errno);

                    /*Receive vars*/
                    MPI_Probe(src, VARS_TAG, gl_comps[i].intercomm, &status);
                    MPI_Get_count(&status, MPI_BYTE, &buf_sz);
                    buf = malloc((size_t)buf_sz);
                    assert(buf != NULL);
                    errno = MPI_Recv(buf, buf_sz, MPI_BYTE, src, VARS_TAG, gl_comps[i].intercomm, &status);
                    CHECK_MPI(errno);

                    unpack_vars(buf_sz, buf, &fbuf->vars, &fbuf->var_cnt);
                    free(buf);
                    /*Receive nodes of all variables*/
                    var = fbuf->vars;
                    while(var != NULL){
                        reqs = realloc(reqs, sizeof(MPI_Request)*(int)var->node_cnt);
                        assert(reqs != NULL);

                        node = var->nodes;
                        FARB_DBG(VERBOSE_DBG_LEVEL, "Will recv %d nodes", (int)var->node_cnt);
                        for(j = 0; j < var->node_cnt; j++){
                            FARB_DBG(VERBOSE_DBG_LEVEL, "node sz %d offt %d", (int)node->data_sz, (int)node->offset);
                            errno = MPI_Irecv(node->data, (int)node->data_sz, MPI_BYTE, src, NODE_TAG+j, gl_comps[i].intercomm, &reqs[j]);
                            CHECK_MPI(errno);
                            node = node->next;
                        }
                        errno = MPI_Waitall(var->node_cnt, reqs, MPI_STATUSES_IGNORE);
                        CHECK_MPI(errno);
                        var = var->next;
                    }
                    free(reqs);
                }

                fbuf->is_ready = 1;
                break;
            case RECV_READY_TAG:

                src = status.MPI_SOURCE;
                errno = MPI_Recv(filename, MAX_FILE_NAME, MPI_CHAR, src, RECV_READY_TAG, gl_comps[i].intercomm, &status);

                CHECK_MPI(errno);
                fbuf = find_file_buffer(gl_filebuf_list, filename);
                assert(fbuf != NULL);
                FARB_DBG(VERBOSE_DBG_LEVEL, "Receive RECV_READY_TAG notif for %s", filename);

                /*Send the hearder*/
                errno = MPI_Isend(fbuf->header, (int)fbuf->hdr_sz, MPI_CHAR, src, HEADER_TAG, gl_comps[i].intercomm, &sreq);
                CHECK_MPI(errno);

                /*Pack the vars info and send it*/
                pack_vars(fbuf->var_cnt, fbuf->vars, &buf_sz, &buf);
                errno = MPI_Wait(&sreq, MPI_STATUS_IGNORE);
                CHECK_MPI(errno);

                errno = MPI_Send(buf, buf_sz, MPI_BYTE, status.MPI_SOURCE, VARS_TAG, gl_comps[i].intercomm);
                CHECK_MPI(errno);
                free(buf);
                /*Send buffer nodes*/
                var = fbuf->vars;
                while(var != NULL){
                    reqs = realloc(reqs, sizeof(MPI_Request)*(int)var->node_cnt);
                    assert(reqs != NULL);

                    node = var->nodes;
                    FARB_DBG(VERBOSE_DBG_LEVEL, "Will send %d nodes", (int)var->node_cnt);
                    for(j = 0; j < var->node_cnt; j++){
                        FARB_DBG(VERBOSE_DBG_LEVEL, "node sz %d offt %d", (int)node->data_sz, (int)node->offset);
                        errno = MPI_Isend(node->data, (int)node->data_sz, MPI_BYTE, status.MPI_SOURCE, NODE_TAG+j, gl_comps[i].intercomm, &reqs[j]);
                        CHECK_MPI(errno);
                        node = node->next;
                    }
                    errno = MPI_Waitall(var->node_cnt, reqs, MPI_STATUSES_IGNORE);
                    CHECK_MPI(errno);
                    var = var->next;
                }
                free(reqs);
                fbuf->transfered++;
                break;
            default:
                FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: unknown tag %d", status.MPI_TAG);
                assert(0);
        }
    }

}

//typedef struct farb_var{
//    int                     id;         /* varid assigned by pnetcdf*/
//    MPI_Offset              *shape;     /* dim->size of each dim */
//    int                     ndims;      /* number of dimensions */
//
//    buffer_node_t           *nodes;      /*head buffer node that stores the data*/
//    MPI_Offset              node_cnt;   /*number of allocated buffer nodes*/
//    struct farb_var *next;
//}farb_var_t;

/*
    vars        [OUT]
    var_cnt     [OUT]
    the rest    [IN]
*/
void unpack_vars(int buf_sz, void *buf, farb_var_t **vars, int *var_cnt)
{
    size_t offt = 0;
    MPI_Offset offset, data_sz;
    assert(buf_sz > 0);
    assert(buf != NULL);

    FARB_DBG(VERBOSE_DBG_LEVEL, "Unpacking vars sz %d", buf_sz);

    int i, j;
    farb_var_t *var;
    buffer_node_t *node;

    memcpy(var_cnt, buf, sizeof(int));
    offt += sizeof(int);

    for(i = 0; i < *var_cnt; i++){
        var = malloc(sizeof(farb_var_t));
        assert(var != NULL);
        var->nodes = NULL;
        var->next = NULL;

        memcpy(&var->id, buf+offt, sizeof(var->id));
        offt += sizeof(var->id);
        memcpy(&var->ndims, buf+offt, sizeof(var->ndims));
        offt += sizeof(var->id);
        if(var->ndims > 0){
            var->shape = malloc(var->ndims*sizeof(MPI_Offset));
            assert(var->shape != NULL);
            for(j = 0; j < var->ndims; j++){
                memcpy(&(var->shape[j]), buf+offt, sizeof(MPI_Offset));
                offt += sizeof(MPI_Offset);
            }
        } else
            var->shape = NULL;
        memcpy(&var->node_cnt, buf+offt, sizeof(var->node_cnt));
        offt += sizeof(var->node_cnt);
        for(j = 0; j < var->node_cnt; j++){
            memcpy(&offset, buf+offt, sizeof(MPI_Offset));
            offt += sizeof(MPI_Offset);
            memcpy(&data_sz, buf+offt, sizeof(MPI_Offset));
            offt += sizeof(MPI_Offset);
            node = new_buffer_node(offset, data_sz);
            insert_buffer_node(&var->nodes, node);
        }

        add_var(vars, var);
    }

    assert(offt == buf_sz);
    FARB_DBG(VERBOSE_DBG_LEVEL, "Finished unpacking vars");
}

/*
    *buf_sz     [OUT]
    **buf       [OUT]
    the rest    [IN]
*/
void pack_vars(int var_cnt, farb_var_t *vars, int *buf_sz, void **buf)
{
    int i;
    size_t sz = 0, offt=0;
    buffer_node_t *node;
    /*Count how much space we need. We'll pack necessary fields from farv_var_t
      + will send the list of offsets of all variable nodes*/
    sz += sizeof(var_cnt);
    farb_var_t *var = vars;
    while(var != NULL){
        sz += sizeof(var->id) + sizeof(MPI_Offset)*var->ndims +
              sizeof(var->ndims) + sizeof(var->node_cnt) + var->node_cnt*sizeof(MPI_Offset)*2;
        var = var->next;
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Packing vars: sz %lu", sz);

    *buf = malloc(sz);
    assert(*buf != NULL);
    //write the number of vars
    memcpy(*buf, &var_cnt, sizeof(int));
    offt += sizeof(int);
    var = vars;
    while(var != NULL){
        memcpy(*buf+offt, &var->id, sizeof(var->id));
        offt += sizeof(var->id);
        memcpy(*buf+offt, &var->ndims, sizeof(var->ndims));
        offt += sizeof(var->ndims);
        for(i = 0; i<var->ndims;i++){
            memcpy(*buf+offt, &var->shape[i], sizeof(MPI_Offset));
            offt += sizeof(MPI_Offset);
        }
        memcpy(*buf+offt, &var->node_cnt, sizeof(var->node_cnt));
        offt+=sizeof(var->node_cnt);
        /*Pack node offsets*/
        node = var->nodes;
        while(node != NULL){
            memcpy(*buf+offt, &node->offset, sizeof(MPI_Offset));
            offt+=sizeof(MPI_Offset);
            memcpy(*buf+offt, &node->data_sz, sizeof(MPI_Offset));
            offt+=sizeof(MPI_Offset);
            node = node->next;
        }

        var = var->next;
    }
    FARB_DBG(VERBOSE_DBG_LEVEL, "packing: offt %d, sz %d", (int)offt, (int)sz);
    assert(offt == sz);

    *buf_sz = (int)sz;
    FARB_DBG(VERBOSE_DBG_LEVEL, "Finish packing vars");
}

void notify_file_ready(const char* filename)
{
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);
    int i, comp_id, dest, errno;

    for(i = 0; i < fbuf->nreaders; i++){
        comp_id = fbuf->reader_ids[i];
        dest = gl_my_rank;  //file per process mode
        FARB_DBG(VERBOSE_DBG_LEVEL,   "Notify reader %s that file %s is ready", gl_comps[comp_id].name, fbuf->file_path);
        errno = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, dest, FILE_READY_TAG, gl_comps[comp_id].intercomm);
        CHECK_MPI(errno);

        if(fbuf->mode == FARB_IO_MODE_FILE)
            fbuf->transfered++;
    }

}

void close_file(const char *filename)
{
    FARB_DBG(VERBOSE_DBG_LEVEL, "Closing file %s", filename);
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);

    if(fbuf->write_flag){
        notify_file_ready(filename);
        while(fbuf->transfered != fbuf->nreaders)
            progress_io();
    }

    delete_file_buffer(&gl_filebuf_list, fbuf);
}

int get_io_mode(const char* filename)
{

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    if(fbuf == NULL)
        return FARB_IO_MODE_UNDEFINED;
    return fbuf->mode;
}

int def_var(const char* filename, int varid, int ndims, MPI_Offset *shape)
{
    file_buffer_t *buf = find_file_buffer(gl_filebuf_list, filename);
    assert(buf!=NULL);

    farb_var_t *var = new_var(varid, ndims, shape);

    add_var(&buf->vars, var);
    buf->var_cnt++;

    return 0;
}

void write_hdr(const char *filename, MPI_Offset hdr_sz, void *header)
{
    file_buffer_t *buf = find_file_buffer(gl_filebuf_list, filename);
    assert(buf!=NULL);

    buf->hdr_sz = hdr_sz;
    buf->header = malloc(hdr_sz);
    assert(buf->header != NULL);
    memcpy(buf->header, header, (size_t)hdr_sz);
    return;
}

MPI_Offset read_hdr_chunk(const char *filename, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
    file_buffer_t *buf = find_file_buffer(gl_filebuf_list, filename);
    assert(buf!=NULL);

    if(offset+chunk_sz > buf->hdr_sz){
        FARB_DBG(VERBOSE_DBG_LEVEL, "Warning: trying to read %llu at offt %llu but hdr sz is %llu", chunk_sz, offset, buf->hdr_sz);
        chunk_sz = buf->hdr_sz - offset;
    }

    memcpy(chunk, buf->header+offset, (size_t)chunk_sz);
    return chunk_sz;
}

static MPI_Offset to_1d_index(int ndim, MPI_Offset *shape, MPI_Offset *coord)
{
    int dimid = 0;
    MPI_Offset idx = coord[dimid];

    while(dimid != ndim - 1){
       idx = idx*shape[dimid+1] + coord[dimid+1];
       dimid++;
    }
    return idx;
}

MPI_Offset read_write_var(const char *filename, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, MPI_Datatype dtype, void *buf, int rw_flag)
{
    MPI_Offset el_offt, byte_offt, lead_el_idx, el_cnt;
    int dimid;
    int el_sz;
    MPI_Offset data_sz;
    MPI_Offset buf_offt = 0, tmp, ret;

    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(buf!=NULL);

    if(rw_flag == FARB_READ)
         assert(fbuf->is_ready);

    if(imap != NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: writing mapped vars is not impelemented yet. Ignore.");
        return 0;
    }

    if(stride != NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: writing vars at a stride is not impelemented yet. Ignore.");
        return 0;
    }

    farb_var_t *var = find_var(fbuf->vars, varid);
    if(var == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: could not find var with id %d", varid);
        return 0;
    }

    MPI_Type_size(dtype, &el_sz);
    assert(el_sz > 0);
    if(var->ndims <= 1){ /*scalar or 1d array*/
        el_offt = *start;
        byte_offt = el_offt*(MPI_Offset)el_sz;
        /*Data size to write*/
        data_sz = (*count)*el_sz;

        if(rw_flag == FARB_READ){
            FARB_DBG(VERBOSE_DBG_LEVEL, "Var %d: will read %llu elems of sz %d from element offt %llu ", var->id, *count, el_sz, el_offt);
            ret = mem_read(var, byte_offt, data_sz, buf);
        } else {
            FARB_DBG(VERBOSE_DBG_LEVEL, "Var %d: will write %llu elems of sz %d to element offt %llu ", var->id, *count, el_sz, el_offt);
            ret = mem_write(var, byte_offt, data_sz, buf);
        }

        if(ret != data_sz)
            FARB_DBG(VERBOSE_DBG_LEVEL, "FARB Warning: Meant to read/write %llu bytes but actual val is %llu", data_sz, ret);

    } else { /*multi-dimentional array*/
        MPI_Offset *start1 = malloc(var->ndims*sizeof(MPI_Offset));
        assert(start1 != NULL);
        memcpy(start1, start, var->ndims*sizeof(MPI_Offset));
        FARB_DBG(VERBOSE_DBG_LEVEL, "start el %lld, until %lld", start1[0], start1[0] + count[0]);
        for(lead_el_idx = start[0]; lead_el_idx < start[0] + count[0]; lead_el_idx++){
            /*Calculate element offset within buffer nodes*/
            /*for an array stored in row-major form, elements with a fixed row coordinate
            (e.g. coord X in XYZ dimensions) will be written contiguously in a 1d array.
            So we just need to find an element offset for a given X and then write Y*Z
            contiguous elements. */
            start1[0] = lead_el_idx;
            el_offt = to_1d_index(var->ndims, var->shape, start1);
            byte_offt = el_offt*el_sz;

            /*how many elements to write contiguously?*/
            dimid = 1;
            el_cnt = count[dimid];
            while(dimid != var->ndims - 1){
                el_cnt = el_cnt*count[dimid+1];
                dimid++;
            }

            data_sz = el_cnt*el_sz;


            if(rw_flag == FARB_READ){
                FARB_DBG(VERBOSE_DBG_LEVEL, "Var %d: will read %llu elems of sz %d from element offt %llu ", var->id, el_cnt, el_sz, el_offt);
                tmp = mem_read(var, byte_offt, data_sz, (void*)(buf+buf_offt));
            } else { /*FARB_WRITE*/
                FARB_DBG(VERBOSE_DBG_LEVEL, "Var %d: will write %llu elems of sz %d to element offt %llu ", var->id, el_cnt, el_sz, el_offt);
                tmp = mem_write(var, byte_offt, data_sz, (void*)(buf+buf_offt));
            }

            if(tmp != data_sz)
                FARB_DBG(VERBOSE_DBG_LEVEL, "FARB Warning: Meant to read/write %llu bytes but actual val is %llu", data_sz, tmp);
            buf_offt+=tmp;
        }
        free(start1);
        ret = tmp;
    }

    return ret;
}
