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
        FARB_DBG(VERBOSE_DBG_LEVEL,   "File %s, version %d, writer %s, reader %s ",
                fb->file_path, fb->version, gl_comps[fb->writer_id].name, gl_comps[fb->reader_id].name);
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
    if(strlen(line) == 0 || line[0] == '#' || line[0] == ';' || line[0] == '!') //comment or empty line
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
            }
            if( cur_fb->writer_id == -1 || cur_fb->reader_id == -1){
                FARB_DBG(VERBOSE_ERROR_LEVEL,"FARB Error: file reader or writer not set for file %s.", cur_fb->file_path);
                goto panic_exit;
            }
            if(cur_fb->distr_rule == DISTR_RULE_RANGE){
                if(cur_fb->distr_range == 0){
                    FARB_DBG(VERBOSE_ERROR_LEVEL, "Farb Error: Please specify range for file %s", cur_fb->file_path);
                    goto panic_exit;
                }
                if(cur_fb->distr_pattern == FARB_UNDEFINED){
                    FARB_DBG(VERBOSE_ERROR_LEVEL, "Farb Error: Please specify distribution pattern for file %s", cur_fb->file_path);
                    goto panic_exit;
                }
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
                gl_comps[i].connect_mode = FARB_UNDEFINED;
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

            if(cur_fb->reader_id != -1){
                FARB_DBG(VERBOSE_ERROR_LEVEL, "File %s cannot have multiple readers", cur_fb->file_path);
                goto panic_exit;
            }
            for(i = 0; i < gl_ncomp; i++){
                if(strcmp(value, gl_comps[i].name) == 0){
                    cur_fb->reader_id = i;
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

        } else if(strcmp(param, "distr_rule") == 0){
            if(strcmp(value, "range") == 0)
                cur_fb->distr_rule = DISTR_RULE_RANGE;
            else if(strcmp(value, "ranks") == 0){
                FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: distribution rule <ranks> is not implemented yet. Use <range> instead.");
                goto panic_exit;
            }
        } else if(strcmp(param, "distr_pattern") == 0){
            if(strcmp(value, "scatter") == 0)
                cur_fb->distr_pattern = DISTR_PATTERN_SCATTER;
            else if(strcmp(value, "all") == 0)
                cur_fb->distr_pattern = DISTR_PATTERN_ALL;

            if(cur_fb->distr_pattern == FARB_UNDEFINED){
                FARB_DBG(VERBOSE_ERROR_LEVEL, "Farb Error: unknown distribution pattern %s.", value);
                goto panic_exit;
            }
        } else if(strcmp(param, "distr_range") == 0){
            cur_fb->distr_range = atoi(value);
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
  file flow then the connection mode is not set(FARB_UNDEFINED).
*/
    cur_fb = gl_filebuf_list;
    while(cur_fb != NULL){
        if(cur_fb->writer_id == gl_my_comp_id){
            if(gl_comps[ cur_fb->reader_id ].connect_mode == FARB_UNDEFINED )
               gl_comps[ cur_fb->reader_id ].connect_mode = CONNECT_MODE_SERVER; // I am server

        } else if( gl_comps[ cur_fb->writer_id ].connect_mode == FARB_UNDEFINED){
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


static int create_intercomm(int comp_id, char* global_path){

    char portname[MPI_MAX_PORT_NAME];
    char service_name[MAX_COMP_NAME*2+1];
    FILE* portfile;
    int errno, mode, myrank, l1, l2;
    char portfile_name[MAX_FILE_NAME];
    int wait_timeout;

    mode = gl_comps[comp_id].connect_mode;
    if(mode == FARB_UNDEFINED) return 0;

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

    if(mode == FARB_UNDEFINED) return;
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

int receive_data(file_buffer_t *fbuf, int rank, MPI_Comm intercomm)
{
    MPI_Status status;
    int buf_sz;
    int errno, i;
    void *buf;
    farb_var_t *var;
    MPI_Request *reqs = NULL;
    buffer_node_t *node;

    /*Receive the header*/
    MPI_Probe(rank, HEADER_TAG, intercomm, &status);
    MPI_Get_count(&status, MPI_BYTE, &buf_sz);
    fbuf->hdr_sz = (MPI_Offset)buf_sz;
    FARB_DBG(VERBOSE_DBG_LEVEL, "Hdr size to receive %d", buf_sz);
    fbuf->header = malloc((size_t)buf_sz);
    assert(fbuf->header != NULL);
    errno = MPI_Recv(fbuf->header, buf_sz, MPI_BYTE, rank, HEADER_TAG, intercomm, &status);
    CHECK_MPI(errno);

    /*Receive vars*/
    MPI_Probe(rank, VARS_TAG, intercomm, &status);
    MPI_Get_count(&status, MPI_BYTE, &buf_sz);
    buf = malloc((size_t)buf_sz);
    assert(buf != NULL);
    errno = MPI_Recv(buf, buf_sz, MPI_BYTE, rank, VARS_TAG, intercomm, &status);
    CHECK_MPI(errno);

    unpack_vars(buf_sz, buf, &fbuf->vars, &fbuf->var_cnt);

    /*Receive nodes of all variables*/
    var = fbuf->vars;
    while(var != NULL){
        reqs = realloc(reqs, sizeof(MPI_Request)*(int)var->node_cnt);
        assert(reqs != NULL);

        node = var->nodes;
        FARB_DBG(VERBOSE_DBG_LEVEL, "Will recv %d nodes", (int)var->node_cnt);
        for(i = 0; i < var->node_cnt; i++){
            FARB_DBG(VERBOSE_DBG_LEVEL, "node sz %d offt %d", (int)node->data_sz, (int)node->offset);
            errno = MPI_Irecv(node->data, (int)node->data_sz, MPI_BYTE, rank, NODE_TAG+i, intercomm, &reqs[i]);
            CHECK_MPI(errno);
            node = node->next;
        }
        errno = MPI_Waitall(var->node_cnt, reqs, MPI_STATUSES_IGNORE);
        CHECK_MPI(errno);
        var = var->next;
    }
    free(reqs);
    free(buf);

    return 0;
}

int send_data(file_buffer_t *fbuf, int rank, MPI_Comm intercomm)
{
    int errno,i;
    int buf_sz;
    void *buf;
    farb_var_t *var;
    MPI_Request sreq, *reqs = NULL;
    buffer_node_t *node;

    /*Send the hearder*/
    errno = MPI_Isend(fbuf->header, (int)fbuf->hdr_sz, MPI_CHAR, rank, HEADER_TAG, intercomm, &sreq);
    CHECK_MPI(errno);

    /*Pack the vars info and send it*/
    pack_vars(fbuf->var_cnt, fbuf->vars, &buf_sz, &buf);
    errno = MPI_Wait(&sreq, MPI_STATUS_IGNORE);
    CHECK_MPI(errno);

    errno = MPI_Send(buf, buf_sz, MPI_BYTE, rank, VARS_TAG, intercomm);
    CHECK_MPI(errno);

    /*Send buffer nodes*/
    var = fbuf->vars;
    while(var != NULL){
        reqs = realloc(reqs, sizeof(MPI_Request)*(int)var->node_cnt);
        assert(reqs != NULL);

        node = var->nodes;
        FARB_DBG(VERBOSE_DBG_LEVEL, "Will send %d nodes", (int)var->node_cnt);
        for(i = 0; i < var->node_cnt; i++){
            FARB_DBG(VERBOSE_DBG_LEVEL, "node sz %d offt %d", (int)node->data_sz, (int)node->offset);
            errno = MPI_Isend(node->data, (int)node->data_sz, MPI_BYTE, rank, NODE_TAG + i, intercomm, &reqs[i]);
            CHECK_MPI(errno);
            node = node->next;
        }
        errno = MPI_Waitall(var->node_cnt, reqs, MPI_STATUSES_IGNORE);
        CHECK_MPI(errno);
        var = var->next;
    }
    free(reqs);
    free(buf);
    return 0;
}

void progress_io()
{
    MPI_Status status;
    int i, flag, src, errno;
    char filename[MAX_FILE_NAME];
    file_buffer_t *fbuf;

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
                    receive_data(fbuf, src, gl_comps[i].intercomm);
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

                send_data(fbuf, src, gl_comps[i].intercomm);

                fbuf->distr_ndone++;
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
        var = new_var(0, 0, NULL);
        //memcpy(&var->id, buf+offt, sizeof(var->id));
        var->id = *((int*)(buf+offt));
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
    int comp_id, dest, errno;

    comp_id = fbuf->reader_id;
    dest = gl_my_rank;  //file per process mode
    FARB_DBG(VERBOSE_DBG_LEVEL,   "Notify reader %s that file %s is ready", gl_comps[comp_id].name, fbuf->file_path);
    errno = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, dest, FILE_READY_TAG, gl_comps[comp_id].intercomm);
    CHECK_MPI(errno);

    if(fbuf->mode == FARB_IO_MODE_FILE)
        fbuf->distr_ndone++;


}

void close_file(const char *filename)
{
    FARB_DBG(VERBOSE_DBG_LEVEL, "Closing file %s", filename);
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);

    if(fbuf->write_flag){
        notify_file_ready(filename);
        //while(fbuf->distr_nranks != fbuf->distr_ndone)
        while(fbuf->distr_ndone !=1)
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
/*Write pnetcdf header*/
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

/*Read pnetcdf header*/
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

/*Set how many elements in each dimension to distribute to ranks.
  File's distrib_pattern must be set to scatter */
int set_distr_count(const char* filename, int varid, int count[])
{
    int i;
    MPI_Offset *cnt;
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);

    if(fbuf->distr_pattern != DISTR_PATTERN_SCATTER){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "Farb Warning: cannot set distribute count. File's distribute pattern must be <scatter>.");
        return 0;
    }

    farb_var_t *var = find_var(fbuf->vars, varid);
    if(var == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "No variable with such id (%d)", varid);
        assert(0);
    }

    cnt = (MPI_Offset*)realloc(var->distr_count, var->ndims*sizeof(MPI_Offset));
    assert(cnt != NULL);
    var->distr_count = cnt;
    for(i = 0; i < var->ndims; i++)
        var->distr_count[i] = (MPI_Offset)count[i];

    return 0;
}

int init_data_distr()
{
//    int nranks;
//
//    file_buffer_t *fbuf = gl_filebuf_list;
//
//    while(fbuf != NULL){
//        if(fbuf->distr_rule == DISTR_RULE_RANGE){
//            //if(fbuf->writer_id == gl_my_comp_id)
//                //MPI_Comm_remote_size(gl_comps)
//        }
//        fbuf = fbuf->next;
//    }
    return 0;
}
