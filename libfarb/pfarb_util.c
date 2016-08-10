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
    int ret, mode, myrank, l1, l2;
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

        ret = MPI_Bcast(portname, MPI_MAX_PORT_NAME, MPI_CHAR, 0, MPI_COMM_WORLD );
        assert(ret == MPI_SUCCESS);
        FARB_DBG(VERBOSE_DBG_LEVEL,   "%s will connect to service %s.", gl_comps[gl_my_comp_id].name, service_name);

        ret = MPI_Comm_connect(portname, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &(gl_comps[comp_id].intercomm));

    } else if(mode == CONNECT_MODE_SERVER){
        FARB_DBG(VERBOSE_DBG_LEVEL,   "Creating connection for %s.", service_name);

        ret = MPI_Open_port(MPI_INFO_NULL, portname);
        assert(ret == MPI_SUCCESS);

        if(myrank == 0){

           portfile = fopen(portfile_name, "wt");
           fprintf(portfile, "%s\n", portname);
           fclose(portfile);
        }

        FARB_DBG(VERBOSE_DBG_LEVEL,   "%s starts listening for service %s.", gl_comps[gl_my_comp_id].name, service_name);

        ret = MPI_Comm_accept(portname, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &(gl_comps[comp_id].intercomm) );
        assert(ret == MPI_SUCCESS);
        FARB_DBG(VERBOSE_DBG_LEVEL,   "%s accepted connection on service %s.", gl_comps[gl_my_comp_id].name, service_name);

        ret = MPI_Close_port(portname);
        assert(ret == MPI_SUCCESS);
    }

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


MPI_Offset mem_write(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data)
{
    MPI_Offset node_id, last_node_id;

    MPI_Offset node_offt, to_cpy, copied=0;
    int i;

    FARB_DBG(VERBOSE_DBG_LEVEL,   "Will write %llu bytes to var %d at offt %llu", data_sz, var->id, offset);

    last_node_id = (offset + data_sz) / gl_sett.node_sz;
    if( (offset + data_sz) % gl_sett.node_sz > 0)
        last_node_id++;

    if(last_node_id > var->bufnode_cnt){
        FARB_DBG(VERBOSE_DBG_LEVEL, "Will allocate %llu nodes", last_node_id - var->bufnode_cnt);
        var->bufnodes = realloc(var->bufnodes, last_node_id * sizeof(buffer_node_t));
        assert(var->bufnodes != NULL);

        for(i = var->bufnode_cnt; i < last_node_id; i++){
            var->bufnodes[i].data = NULL;
            var->bufnodes[i].data_sz = 0;
        }

        var->bufnode_cnt = last_node_id;
    }

    while(copied != data_sz){
        node_id = ( offset + copied ) / gl_sett.node_sz ;
        node_offt = ( offset + copied ) % gl_sett.node_sz;

        if(data_sz < gl_sett.node_sz - node_offt)
            to_cpy = data_sz;
        else
            to_cpy = gl_sett.node_sz - node_offt;

        if(var->bufnodes[node_id].data == NULL){
            var->bufnodes[node_id].data = (void*)malloc(gl_sett.node_sz);
            assert(var->bufnodes[node_id].data != NULL);
        }
        memcpy( (void*)(var->bufnodes[node_id].data+node_offt), (void*)(data+copied), to_cpy );
        copied += to_cpy;
    }
    if(offset+data_sz > var->data_sz)
        var->data_sz = offset+data_sz;
    FARB_DBG(VERBOSE_DBG_LEVEL, "Leave mem_write");
    return copied;
}

MPI_Offset mem_read(farb_var_t *var, MPI_Offset offset,  MPI_Offset data_sz, void *data)
{

    MPI_Offset copied = 0, node_offt, to_cpy;
    MPI_Offset new_data_sz = data_sz;
    MPI_Offset node_id;

    FARB_DBG(VERBOSE_DBG_LEVEL, "mem_read at offt %ld of size %lu", (long int)offset, (long unsigned)data_sz);

  //  FARB_DBG(VERBOSE_DBG_LEVEL, "Will read %ld from offt %ld (file sz %ld)", data_sz, offset, fbuf->data_sz);
    if(offset > var->data_sz){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Warning: ps tries to read %llu from offt %llu but file sz is %llu", data_sz, offset, var->data_sz);
        return 0;
    }

    if(offset + data_sz > var->data_sz){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Warning: ps tries to read %llu from offt %llu but file sz is %llu", data_sz, offset, var->data_sz);
        new_data_sz = var->data_sz - offset;
    }

    while(copied != new_data_sz){
        node_id =  ( offset + copied ) / gl_sett.node_sz;
        node_offt =( offset + copied ) % gl_sett.node_sz;

        if(new_data_sz - copied < gl_sett.node_sz)
            to_cpy = new_data_sz - copied;
        else
            to_cpy = gl_sett.node_sz - node_offt;

        assert(var->bufnodes[node_id].data != NULL);
        //FARB_DBG(VERBOSE_DBG_LEVEL, "Will copy %ld from node %u oft %ld to offt %ld", to_cpy, node_id, node_offt, copied );
        memcpy( (void*)(data+copied), (void*)(var->bufnodes[node_id].data+node_offt), to_cpy );
        copied += to_cpy;
    }

    return copied;
}

void progress_io()
{
    MPI_Status status;
    int i, flag, src, errno;
    char filename[MAX_FILE_NAME];
    file_buffer_t *fbuf;
 //   unsigned int nexpect_msg, nmsg = 0;
    struct msg_ready_notif msg_ready;
 //   off_t node_offt;

    for(i = 0; i < gl_ncomp; i++){
        if( i == gl_my_comp_id || gl_comps[i].intercomm == MPI_COMM_NULL)
            continue;

        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, gl_comps[i].intercomm, &flag, &status);
        if(!flag)
            continue;

        switch(status.MPI_TAG){
            case FILE_READY_TAG:
                src = status.MPI_SOURCE;
                errno = MPI_Recv(&msg_ready, 1, msg_ready_datatype, src, FILE_READY_TAG, gl_comps[i].intercomm, &status);
                assert(errno==MPI_SUCCESS);
                FARB_DBG(VERBOSE_DBG_LEVEL,   "Receive FILE_READY notif for %s", msg_ready.filename);

                //receive the file size
                //errno = MPI_Recv(&file_sz, 1, MPI_UNSIGNED_INT, src, FILE_SZ_TAG, gl_comps[i].intercomm, &status);

                fbuf = find_file_buffer(gl_filebuf_list, msg_ready.filename);
                assert(fbuf != NULL);

                if(fbuf->mode == FARB_IO_MODE_MEMORY){
//                    if(msg_ready.file_sz == 0){
//                        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Warning: size of file %s is zero", msg_ready.filename );
//                        fbuf->is_ready = 1;
//                        return;
//                    }
//
//                    //allocate memory to receive
//                    nnodes = (unsigned int)( msg_ready.file_sz / gl_sett.node_sz);
//                    if( msg_ready.file_sz % gl_sett.node_sz > 0)
//                        nnodes++;
//
//                    assert(fbuf->nodes_tbl == NULL);
//                    fbuf->nodes_tbl = (buffer_node_t*)malloc(nnodes * sizeof(buffer_node_t));
//                    assert(fbuf->nodes_tbl != NULL);
//
//                    fbuf->node_cnt = nnodes;
//                    for(j = 0; j < fbuf->node_cnt; j++){
//                        fbuf->nodes_tbl[j].data = malloc(gl_sett.node_sz);
//                        assert(fbuf->nodes_tbl[j].data !=NULL);
//                        fbuf->nodes_tbl[j].data_sz = 0;
//                    }
//
//                    FARB_DBG(VERBOSE_DBG_LEVEL, "Notify that I am ready to receive");
//                    //notify the file writer that we are ready to receive the data
//                    errno = MPI_Send(msg_ready.filename, MAX_FILE_NAME, MPI_CHAR, src, RECV_READY_TAG, gl_comps[i].intercomm);
//                    assert(errno==MPI_SUCCESS);
//
//                    nexpect_msg = (unsigned int)(msg_ready.file_sz / gl_sett.msg_sz);
//                    if( msg_ready.file_sz  % gl_sett.msg_sz > 0)
//                        nexpect_msg++;
//                    FARB_DBG(VERBOSE_DBG_LEVEL, "Expect to receive %d messages", nexpect_msg);
//                    //TODO remove data_sz in buffer nodes
//                    //receive the actual file
//                    node_id = 0;
//                    nmsg = 0;
//                    while(nmsg != nexpect_msg){
//                        node_offt = (nmsg % (unsigned int)(gl_sett.node_sz / gl_sett.msg_sz)) * (size_t)gl_sett.msg_sz;
//
//                        errno = MPI_Recv(fbuf->nodes_tbl[node_id].data + node_offt, gl_sett.msg_sz, MPI_CHAR, src, DATA_TAG+nmsg, gl_comps[i].intercomm, &status);
//                        assert(errno == MPI_SUCCESS);
//                        nmsg++;
//                        node_id = (unsigned int)((gl_sett.msg_sz*nmsg) / gl_sett.node_sz);
//                        FARB_DBG(VERBOSE_DBG_LEVEL, "Received %d", nmsg);
//                    }
//
//                    fbuf->data_sz = (size_t)msg_ready.file_sz;
                }
                fbuf->is_ready = 1;
                break;
            case RECV_READY_TAG:

                src = status.MPI_SOURCE;
                errno = MPI_Recv(&filename, MAX_FILE_NAME, MPI_CHAR, src, RECV_READY_TAG, gl_comps[i].intercomm, &status);
                assert(errno==MPI_SUCCESS);
                fbuf = find_file_buffer(gl_filebuf_list, filename);
                assert(fbuf != NULL);
                FARB_DBG(VERBOSE_DBG_LEVEL, "Receive RECV_READY_TAG notif for %s", filename);
//                //Send the actual file
//                nexpect_msg = (unsigned int)(fbuf->data_sz / gl_sett.msg_sz);
//                if( fbuf->data_sz % gl_sett.msg_sz > 0)
//                    nexpect_msg++;
//                FARB_DBG(VERBOSE_DBG_LEVEL, "Will send %d messages", nexpect_msg);
//                //send the actual file
//                node_id = 0;
//                nmsg = 0;
//                while(nmsg != nexpect_msg){
//                    node_offt = (nmsg % (unsigned int)(gl_sett.node_sz / gl_sett.msg_sz)) * (size_t)gl_sett.msg_sz;
//
//                    errno = MPI_Send(fbuf->nodes_tbl[node_id].data + node_offt, gl_sett.msg_sz, MPI_CHAR, src, DATA_TAG+nmsg, gl_comps[i].intercomm);
//                    assert(errno == MPI_SUCCESS);
//                    nmsg++;
//                    node_id = (unsigned int)((gl_sett.msg_sz*nmsg) / gl_sett.node_sz);
//                    FARB_DBG(VERBOSE_DBG_LEVEL, "Sent %d", nmsg);
//                }

                fbuf->transfered++;
                break;
            default:
                FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: unknown tag %d", status.MPI_TAG);
                assert(0);
        }

    }

}

void notify_file_ready(const char* filename)
{
    file_buffer_t* fbuf = find_file_buffer(gl_filebuf_list, filename);
    assert(fbuf != NULL);
    int i, comp_id, dest, errno;
    struct msg_ready_notif msg_ready;
    strcpy(msg_ready.filename, fbuf->file_path);
    msg_ready.file_sz = 0;//(unsigned int)fbuf->data_sz;

    for(i = 0; i < fbuf->nreaders; i++){
        comp_id = fbuf->reader_ids[i];
        dest = gl_my_rank;  //file per process mode
        FARB_DBG(VERBOSE_DBG_LEVEL,   "Notify reader %s that file %s is ready", gl_comps[comp_id].name, fbuf->file_path);
        //TODO replace with a nonblocking send???
        FARB_DBG(VERBOSE_DBG_LEVEL, "name %s ,size %u", msg_ready.filename, msg_ready.file_sz);
        errno = MPI_Send(&msg_ready, 1, msg_ready_datatype, dest, FILE_READY_TAG, gl_comps[comp_id].intercomm);
        assert(errno == MPI_SUCCESS);

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

    farb_var_t *var = fbuf->vars;
    while(var != NULL){
        FARB_DBG(VERBOSE_DBG_LEVEL, "varid %d, ndims %d", var->id, var->ndims);
        var = var->next;
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

    if(buf->vars == NULL)
        buf->vars = var;
    else{
        farb_var_t *tmp = buf->vars;
        while(tmp->next != NULL)
            tmp = tmp->next;
        tmp->next = var;
    }
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
            } else {
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
