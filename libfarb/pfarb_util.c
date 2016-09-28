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
                    FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: Please specify range for file %s", cur_fb->file_path);
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
            else if (strcmp(value, "p2p") == 0){
                cur_fb->distr_rule = DISTR_RULE_P2P;
            } else {
                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: distribution rule %s is unknown", value);
                goto panic_exit;
            }
        } else if(strcmp(param, "distr_pattern") == 0){
            if(strcmp(value, "scatter") == 0)
                cur_fb->distr_pattern = DISTR_PATTERN_SCATTER;
            else if(strcmp(value, "all") == 0)
                cur_fb->distr_pattern = DISTR_PATTERN_ALL;
            else {
                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: unknown distribution pattern %s.", value);
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
    buffer_node_t *node;
    MPI_Offset node_cnt;
    int ncnt = 0;
    MPI_Offset *first_el_coord;

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

    unpack_vars(fbuf, buf_sz, buf, &node_cnt, &first_el_coord);
    free(buf);
    /*Receive memory nodes*/
    var = fbuf->vars;
    if(fbuf->distr_pattern == DISTR_PATTERN_ALL) {
        MPI_Request *reqs = NULL, *tmp;

        while(var != NULL){
            tmp = realloc(reqs, sizeof(MPI_Request)*node_cnt);
            assert(tmp != NULL);
            reqs = tmp;


            FARB_DBG(VERBOSE_DBG_LEVEL, "Will recv %d nodes", (int)node_cnt);
            i = 0;
            node = var->nodes;
            while(i != node_cnt){
                if(node->data != NULL){
                    node = node->next;
                    continue;
                }
                FARB_DBG(VERBOSE_ALL_LEVEL, "Receive node sz %d offt %d", (int)node->data_sz, (int)node->offset);
                node->data = malloc((size_t)node->data_sz);
                assert(node->data != NULL);
                errno = MPI_Irecv(node->data, (int)node->data_sz, MPI_BYTE, rank, NODE_TAG+ncnt, intercomm, &reqs[i]);
                CHECK_MPI(errno);
                ncnt++;
                node = node->next;
                i++;
            }
            errno = MPI_Waitall(node_cnt, reqs, MPI_STATUSES_IGNORE);
            CHECK_MPI(errno);
            var = var->next;
        }
        free(reqs);
    } else { /*DISTR_PATTERN_SCATTER*/
        void *rbuf = NULL, *tmp;
        int bufsz;
        int first_el_coord_sz = 0;
        MPI_Offset written;

        while(var != NULL){
            if(var->distr_count == NULL){
                var = var->next;
                continue;
            }
            bufsz = 1;
            for(i = 0; i < var->ndims; i++)
                bufsz *= (int)var->distr_count[i];
            bufsz *= (int)var->el_sz;
            tmp = realloc(rbuf, bufsz);
            assert(tmp != NULL);
            rbuf = tmp;
            FARB_DBG(VERBOSE_ALL_LEVEL, "Try to recv node of sz %d", (int)bufsz);
            errno = MPI_Recv(rbuf, bufsz, MPI_BYTE, rank, NODE_TAG+ncnt, intercomm, MPI_STATUS_IGNORE);
            CHECK_MPI(errno);

            written = mem_noncontig_write(var, &first_el_coord[first_el_coord_sz], var->distr_count, rbuf);
            assert((int)written == bufsz);
            first_el_coord_sz += (int)var->ndims;
            ncnt++;
            var = var->next;
        }
        free(rbuf);
    }
    assert(ncnt == (int)node_cnt);
    return 0;
}

int send_data(file_buffer_t *fbuf, int rank, MPI_Comm intercomm)
{
    int errno,i;
    int buf_sz;
    void *buf;
    farb_var_t *var;
    MPI_Request sreq;
    buffer_node_t *node;
    MPI_Offset node_count;
    int ncnt = 0;
    int first_el_coord_sz = 0;
    MPI_Offset *first_el_coord;

    /*Send the hearder*/
    errno = MPI_Isend(fbuf->header, (int)fbuf->hdr_sz, MPI_CHAR, rank, HEADER_TAG, intercomm, &sreq);
    CHECK_MPI(errno);

    /*Pack the vars info and send it*/
    pack_vars(fbuf, rank, &buf_sz, &buf, &node_count, &first_el_coord);
    errno = MPI_Wait(&sreq, MPI_STATUS_IGNORE);
    CHECK_MPI(errno);
    errno = MPI_Send(buf, buf_sz, MPI_BYTE, rank, VARS_TAG, intercomm);
    CHECK_MPI(errno);

    /*Send buffer nodes*/
    var = fbuf->vars;

    if(fbuf->distr_pattern == DISTR_PATTERN_ALL) {
        MPI_Request *reqs = NULL, *tmp;
        while(var != NULL){
            tmp = (MPI_Request*)realloc(reqs, sizeof(MPI_Request)*(int)var->node_cnt);
            assert(tmp != NULL);
            reqs = tmp;

            node = var->nodes;
            FARB_DBG(VERBOSE_DBG_LEVEL, "Will send %d nodes", (int)var->node_cnt);
            for(i = 0; i < var->node_cnt; i++){
                FARB_DBG(VERBOSE_ALL_LEVEL, "node sz %d offt %d", (int)node->data_sz, (int)node->offset);
                errno = MPI_Isend(node->data, (int)node->data_sz, MPI_BYTE, rank, NODE_TAG + ncnt, intercomm, &reqs[i]);
                CHECK_MPI(errno);
                node = node->next;
                ncnt++;
            }
            errno = MPI_Waitall(var->node_cnt, reqs, MPI_STATUSES_IGNORE);
            CHECK_MPI(errno);
            var = var->next;
        }
        free(reqs);
    } else { /*SCATTER*/
        void *sbuf = NULL, *tmp;
        int bufsz;
        MPI_Offset readsz;

        while(var != NULL){
            if(var->distr_count == NULL){
                var = var->next;
                continue;
            }
            //pack multi-dim data into 1d contiguous buffer
            bufsz = 1;
            for(i = 0; i < var->ndims; i++)
                bufsz *= (int)var->distr_count[i];
            bufsz *= (int)var->el_sz;

            tmp = realloc(sbuf, bufsz);
            assert(tmp != NULL);
            sbuf = tmp;
            for(i = 0; i < var->ndims; i++){
                FARB_DBG(VERBOSE_ALL_LEVEL, "send data starts at %d", (int)first_el_coord[first_el_coord_sz+i]);
            }
            readsz = mem_noncontig_read(var, &first_el_coord[first_el_coord_sz], var->distr_count, sbuf);
            if( (int)readsz != bufsz ){
                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Warning: not enough data to send to rank %d", rank);
                FARB_DBG(VERBOSE_DBG_LEVEL, "Should read %d but read %d", (int)bufsz, (int)readsz);
            }
            //if(readsz != 0){
            FARB_DBG(VERBOSE_ALL_LEVEL, "Sending node of sz %d", (int)readsz);
            errno = MPI_Send(sbuf, (int)readsz, MPI_BYTE, rank, NODE_TAG+ncnt, intercomm);
            CHECK_MPI(errno);
            ncnt++;
            first_el_coord_sz += var->ndims;
            //}
            var = var->next;
        }

        free(sbuf);
    }
    FARB_DBG(VERBOSE_DBG_LEVEL, "Total sent mem nodes %d, should have sent %d", ncnt, (int)node_count);
    assert((int)node_count == ncnt);
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
                    fbuf->distr_ndone++;
                } else
                    fbuf->distr_ndone++;
                if(fbuf->distr_ndone == fbuf->distr_nranks)
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

static void print_var(farb_var_t *var)
{
    FARB_DBG(VERBOSE_DBG_LEVEL, "Var id %d, ndims %d", var->id, var->ndims);
}
/*
    vars        [OUT]
    var_cnt     [OUT]
    the rest    [IN]
*/
void unpack_vars(file_buffer_t *fbuf, int buf_sz, void *buf, MPI_Offset *node_cnt, MPI_Offset **first_el_coord)
{
    size_t offt = 0;
    MPI_Offset offset, data_sz, ncnt;
    int var_cnt;
    assert(buf_sz > 0);
    assert(buf != NULL);
    *node_cnt = 0;
    *first_el_coord = NULL;
    int first_el_coord_sz = 0;

    FARB_DBG(VERBOSE_DBG_LEVEL, "Unpacking vars sz %d", buf_sz);

    int i, j, varid;
    farb_var_t *var;
    buffer_node_t *node;

    var_cnt = *((int*)buf);
    offt += sizeof(int);
    FARB_DBG(VERBOSE_DBG_LEVEL, "unpack nvars %d", var_cnt);
    for(i = 0; i < var_cnt; i++){

        varid = *((int*)(buf+offt));

        offt += sizeof(int);
        var = find_var(fbuf->vars, varid);
        if(var == NULL){
            var = new_var(varid, 0, NULL);
            var->el_sz = *((MPI_Offset*)(buf+offt));
            offt+=sizeof(MPI_Offset);
            var->ndims = *((int*)(buf+offt));
            offt += sizeof(int);

            if(var->ndims > 0){
                var->shape = malloc(var->ndims*sizeof(MPI_Offset));
                assert(var->shape != NULL);
                for(j = 0; j < var->ndims; j++){
                    var->shape[j] = *((MPI_Offset*)(buf+offt));
                    offt += sizeof(MPI_Offset);
                }
            } else
                var->shape = NULL;

            add_var(&(fbuf->vars), var);
            fbuf->var_cnt++;
        } else {
            /*Skip the fields that have already been defined*/
            offt += sizeof(MPI_Offset) + sizeof(int) + sizeof(MPI_Offset)*var->ndims;
        }

        print_var(var);

        ncnt = *((MPI_Offset*)(buf+offt));
        offt += sizeof(MPI_Offset);
        /*Unpack data about nodes to receive*/
        if(fbuf->distr_pattern == DISTR_PATTERN_ALL) {
            var->node_cnt += ncnt;
            *node_cnt += ncnt;

            /*Unpack node offsets*/
            for(j = 0; j < (int)ncnt; j++){
                offset = *((MPI_Offset*)(buf+offt));
                offt += sizeof(MPI_Offset);
                data_sz = *((MPI_Offset*)(buf+offt));
                offt += sizeof(MPI_Offset);
                node = new_buffer_node(offset, data_sz, 0, NULL, 0); //reader does not need to know coordinate of the first element
                insert_buffer_node(&var->nodes, node);
            }
        } else { /*DIST_PATTERN_SCATTER*/
            if(ncnt != 0) {
                (*node_cnt)++;
                /*Only unpack the info about distr_count[] and first_coord[]*/
                MPI_Offset *tmp = (MPI_Offset*)realloc(var->distr_count, var->ndims*sizeof(MPI_Offset));
                assert(tmp != NULL);
                var->distr_count = tmp;
                memcpy((void*)var->distr_count, buf+offt, var->ndims*sizeof(MPI_Offset));
                offt += var->ndims*sizeof(MPI_Offset);

                first_el_coord_sz += var->ndims;
                tmp = (MPI_Offset*)realloc(*first_el_coord, first_el_coord_sz*sizeof(MPI_Offset));
                assert(tmp != NULL);
                *first_el_coord = tmp;

                // need to remember the coordinate of the first (corner) element to be able to unpack
                // contiguous 1D buffer into multi-dimensional nodes
                for(j = 0; j < var->ndims; j++){
                    (*first_el_coord)[first_el_coord_sz - var->ndims + j] = *((MPI_Offset*)(buf+offt));
                    offt += sizeof(MPI_Offset);
                }
            } else {
                if(var->distr_count != NULL) //could be allocated during recv from another writer rank
                    free(var->distr_count);
                var->distr_count = NULL;
            }
        }
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Unpacking vars: offt is %lu, bufsz is %lu", offt, (size_t)buf_sz);
    assert(offt == (size_t)buf_sz);
    FARB_DBG(VERBOSE_DBG_LEVEL, "Finished unpacking vars");
}


MPI_Offset to_1d_index(int ndims, const MPI_Offset *shape, MPI_Offset *coord)
{
      int i, j;
      MPI_Offset idx = 0, mem=0;

      if(ndims == 0) //scalar
        return 0;
      else if(ndims == 1) //1d array
        return *coord;

      for(i = 0; i < ndims; i++){
        mem = coord[i];
        for(j = i+1; j < ndims; j++)
          mem *= shape[j];
        idx += mem;
      }
    return idx;
}

void pack_vars(file_buffer_t *fbuf, int dst_rank, int *buf_sz, void **buf, MPI_Offset *node_cnt, MPI_Offset **first_el_coord)
{
    int i;
    size_t sz = 0, offt=0;
    buffer_node_t *node;
    sz += sizeof(fbuf->var_cnt);
    *node_cnt = 0;
    *first_el_coord = NULL;
    int first_el_coord_sz = 0;

    farb_var_t *var = fbuf->vars;
    while(var != NULL){
        sz += sizeof(var->id) + sizeof(var->el_sz) + sizeof(var->ndims) + sizeof(MPI_Offset)*var->ndims +
        sizeof(MPI_Offset) + var->node_cnt*sizeof(MPI_Offset)*2;
        var = var->next;
    }

    FARB_DBG(VERBOSE_DBG_LEVEL, "Packing vars: sz %lu", sz);

    *buf = malloc(sz);
    assert(*buf != NULL);
    *node_cnt = 0;
    //write the number of vars
    *((int*)(*buf)) = fbuf->var_cnt;
    offt += sizeof(int);

    var = fbuf->vars;
    while(var != NULL){
        *((int*)(*buf+offt)) = var->id;
        offt += sizeof(int);
        *((MPI_Offset*)(*buf+offt)) = var->el_sz;
        offt += sizeof(MPI_Offset);
        *((int*)(*buf+offt)) = var->ndims;
        offt += sizeof(int);

        for(i = 0; i<var->ndims;i++){
            *((MPI_Offset*)(*buf+offt)) = var->shape[i];
            offt += sizeof(MPI_Offset);
        }

        if(fbuf->distr_pattern == DISTR_PATTERN_ALL) {
            *((MPI_Offset*)(*buf+offt)) = var->node_cnt;
            offt+=sizeof(MPI_Offset);
            *node_cnt += var->node_cnt;
            /*Pack node offsets*/
            node = var->nodes;
            while(node != NULL){
                *((MPI_Offset*)(*buf+offt)) = node->offset;
                offt+=sizeof(MPI_Offset);
                *((MPI_Offset*)(*buf+offt)) = node->data_sz;
                offt+=sizeof(MPI_Offset);
                node = node->next;
            }
        } else { /*DIST_PATTERN_SCATTER*/
            /* We will distribute blocks of count[] elemenets one by one
               to each process in the range. If we fail to read
               exactly count[] elements, it means some values are
               lacking. This will cause the app to abort.
               Writer process may have
               written the blocks non-contiguously, so we need
               to find out the start element of the block that will be sent
               to given dst_rank.
               //TODO figure out with fill values in pnetcdf
            */

            //TODO what if it's a scalar var and we need to send it?
            if(var->distr_count == NULL){
                /*If distr_count wasn't set we won't be distributed this variable to
                  the reader. But print out warning just in case if user forgot
                  to set distr_count for this variable */
                  FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Warning: farb_set_distr_count was not called for variable with id %d. \
                  This variable will not be distributed.", var->id);
                  *((MPI_Offset*)(*buf+offt)) = 0ull;
                  offt+=sizeof(MPI_Offset);
            } else {
                int j;
                *((MPI_Offset*)(*buf+offt)) = 1ull;
                offt+=sizeof(MPI_Offset);
                (*node_cnt)++;


                MPI_Offset offset, index_1d;
                int rank_idx;
                MPI_Offset last_el_idx = to_1d_index(var->ndims, var->shape, var->shape) - 1;

                first_el_coord_sz += var->ndims;
                MPI_Offset *coord = (MPI_Offset*)realloc(*first_el_coord, first_el_coord_sz*sizeof(MPI_Offset));
                assert(coord != NULL);
                *first_el_coord = coord;

                coord = &( (*first_el_coord)[first_el_coord_sz - var->ndims] );
                /*First, find out the coordinate of the start element to send*/
                memcpy((void*)coord, (void*)var->nodes->first_coord, var->ndims*sizeof(MPI_Offset));

                for(j = 0; j < var->ndims; j++){
                    FARB_DBG(VERBOSE_ALL_LEVEL, "start with coord %d", (int)coord[j]);
                }
                offset = var->nodes->offset;
                for(rank_idx = 0; rank_idx < fbuf->distr_nranks; rank_idx++){
                    while(1){
                        index_1d = to_1d_index(var->ndims, var->shape, coord);
                        FARB_DBG(VERBOSE_ALL_LEVEL, "1d idx %d", (int)index_1d);
                        /*Check that we haven' t gone out of bounds*/
                        if(index_1d > last_el_idx){
                            FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: There is no data to send to rank %d. Incorrect config file?", dst_rank);
                            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
                        }
                        offset = var->el_sz * index_1d;
                        FARB_DBG(VERBOSE_ALL_LEVEL, "offset %d", (int)offset);
                        if(data_present_at_offt(var->nodes, offset))
                            break;

                        /*find next offset at which data is present*/
                        for(i = 0; i < var->ndims; i++)
                            coord[i] += var->distr_count[i];
                    }

                    if(fbuf->distr_ranks[rank_idx] == dst_rank)
                        break;
                }

                for(j = 0; j < var->ndims; j++){
                    FARB_DBG(VERBOSE_ALL_LEVEL, "data starts at %d, count %d, shape %d", (int)coord[j], (int)var->distr_count[j], (int)var->shape[j]);
                }
                FARB_DBG(VERBOSE_DBG_LEVEL, "Will send to rank %d data starting at idx %d", dst_rank, (int)index_1d);
                assert(rank_idx != fbuf->distr_nranks);
                /*Will eventually pack multi-dimensional data to one contiguous
                  buffer to send.
                  Save the distr_count[] and first_coord[] values so that
                  the reader could unpack the data correctly later*/

                for(i=0; i < var->ndims; i++){
                    *((MPI_Offset*)(*buf+offt)) = var->distr_count[i];
                    offt+=sizeof(MPI_Offset);
                }
                for(i=0; i < var->ndims; i++){
                    *((MPI_Offset*)(*buf+offt)) = coord[i];
                    offt+=sizeof(MPI_Offset);
                }
            }
        }
        print_var(var);
        var = var->next;
    }

    //assert(offt == sz);

    *buf_sz = (int)offt;
    FARB_DBG(VERBOSE_DBG_LEVEL, "Actually packed %d", *buf_sz);
    FARB_DBG(VERBOSE_DBG_LEVEL, "Finish packing vars");
}

void notify_file_ready(const char* filename)
{
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);
    int comp_id, i, dest, errno;

    comp_id = fbuf->reader_id;
    for(i = 0; i < fbuf->distr_nranks; i++){
        dest = fbuf->distr_ranks[i];
        FARB_DBG(VERBOSE_DBG_LEVEL,   "Notify reader rank %d that file %s is ready", dest, fbuf->file_path);
        errno = MPI_Send(fbuf->file_path, MAX_FILE_NAME, MPI_CHAR, dest, FILE_READY_TAG, gl_comps[comp_id].intercomm);
        CHECK_MPI(errno);

        if(fbuf->mode == FARB_IO_MODE_FILE)
            fbuf->distr_ndone++;
    }

}

void close_file(const char *filename)
{
    FARB_DBG(VERBOSE_DBG_LEVEL, "Closing file %s", filename);
    file_buffer_t *fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);

    if(fbuf->write_flag){
        notify_file_ready(filename);
        while(fbuf->distr_ndone != fbuf->distr_nranks)
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

    if(fbuf->mode != FARB_IO_MODE_MEMORY){
        //nothing to do
        return 0;
    }

    if(fbuf->distr_pattern != DISTR_PATTERN_SCATTER){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Warning: cannot set distribute count. File's distribute pattern must be <scatter>.");
        return 0;
    }

    farb_var_t *var = find_var(fbuf->vars, varid);
    if(var == NULL){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Warning: No variable with such id (%d)", varid);
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
    int nranks;

    file_buffer_t *fbuf = gl_filebuf_list;

    while(fbuf != NULL){
        if(fbuf->mode != FARB_IO_MODE_MEMORY)
            fbuf->distr_rule = DISTR_RULE_P2P;

        if(fbuf->distr_rule == DISTR_RULE_P2P){
            fbuf->distr_nranks = 1;
            fbuf->distr_ranks = malloc(sizeof(int));
            assert(fbuf->distr_ranks != NULL);
            *(fbuf->distr_ranks) = gl_my_rank;
        } else if(fbuf->distr_rule == DISTR_RULE_RANGE){
            if(fbuf->writer_id == gl_my_comp_id)
                MPI_Comm_remote_size(gl_comps[fbuf->reader_id].intercomm, &nranks);
            else
                MPI_Comm_remote_size(gl_comps[fbuf->writer_id].intercomm, &nranks);

            fbuf->distr_nranks = (int)(nranks/fbuf->distr_range);
            assert(fbuf->distr_nranks > 0);
            fbuf->distr_ranks = (int*)malloc(fbuf->distr_nranks*sizeof(int));
            assert(fbuf->distr_ranks != NULL);
            fbuf->distr_ranks[0] = gl_my_rank % fbuf->distr_range; //
            int i = 1;
            while(fbuf->distr_ranks[i-1] + fbuf->distr_range < nranks){
                fbuf->distr_ranks[i] = fbuf->distr_ranks[i-1] + fbuf->distr_range;
                i++;
            }
        }
        FARB_DBG(VERBOSE_DBG_LEVEL, "Number of ranks to distribute data to/from %d", fbuf->distr_nranks);
        fbuf = fbuf->next;
    }
    return 0;
}
