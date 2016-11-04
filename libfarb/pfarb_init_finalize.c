#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "pfarb_init_finalize.h"
#include "pfarb_common.h"
#include "pfarb_file_buffer.h"
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

static int create_intercomm(int comp_id, char* global_path){

    char portname[MPI_MAX_PORT_NAME];
    char service_name[MAX_COMP_NAME*2+1];
    FILE* portfile;
    int errno, mode, myrank, l1, l2;
    char portfile_name[MAX_FILE_NAME];
    int wait_timeout;
    int nranks1, nranks2;

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
        return 1;
     }

     l1 = strlen(global_path);
     l2 = strlen(service_name);

     if(l1+l2 > MAX_FILE_NAME){
        FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: Global file path or component name too long.");
        return 1;
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

                if(wait_timeout == FARB_TIMEOUT){
                    FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: timed out waiting for port name in %s.", portfile_name);
                    return 1;
                }
            }
            fclose(portfile);

            if(portname[0] == '\n'){
                FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: no port specified in the portfile.");
                return 1;
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
    /*Check that the number of ranks is the same in both communicators*/
    MPI_Comm_size(MPI_COMM_WORLD, &nranks1);
    MPI_Comm_remote_size(gl_comps[comp_id].intercomm, &nranks2);
    if(nranks1 != nranks2){
        FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB can only work if the number of processes in all components is the same. Aborting.");
        return 1;
    }
    return 0;
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


void clean_config(){
    file_buffer_t *tmp;

    if(gl_comps != NULL)
        free(gl_comps);

    while(gl_filebuf_list != NULL){
        tmp = gl_filebuf_list;
        delete_file_buffer(&gl_filebuf_list, tmp);
    }
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

  gl_conf.distr_mode = FARB_UNDEFINED;

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
            //TODO file does not have to have a reader if the io mode is "file"
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
        } else if(strcmp(param, "distr_mode") == 0){
            if(strcmp(value, "static") == 0)
                gl_conf.distr_mode = DISTR_MODE_STATIC;
//            else if(strcmp(value, "buf_req_match") == 0)
//                gl_conf.distr_mode = DISTR_MODE_BUFFERED_REQ_MATCH;
            else if(strcmp(value, "nbuf_req_match") == 0)
                gl_conf.distr_mode = DISTR_MODE_NONBUFFERED_REQ_MATCH;
            else{
                FARB_DBG(VERBOSE_ERROR_LEVEL, "FARB Error: 3 unknown data distribution mode");
                goto panic_exit;
            }
        }else if(strcmp(param, "filename") == 0){
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
                FARB_DBG(VERBOSE_ERROR_LEVEL, "File %s cannot have multiple readers", cur_fb->file_path);
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
                cur_fb->mode = FARB_IO_MODE_FILE;
            else if(strcmp(value, "memory") == 0){

                if(gl_conf.distr_mode == FARB_UNDEFINED){
                    FARB_DBG(VERBOSE_ERROR_LEVEL,   "FARB Error: data distribution mode is not defined.");
                    goto panic_exit;
                }

                cur_fb->mode = FARB_IO_MODE_MEMORY;
            }

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
        if(i == gl_my_comp_id){
            MPI_Comm_dup(MPI_COMM_WORLD, &(gl_comps[i].intercomm));
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

    for(i = 0; i<gl_ncomp; i++){
       destroy_intercomm(gl_comps[i].id);
    }
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
//            if(fbuf->writer_id == gl_my_comp_id)
//                MPI_Comm_remote_size(gl_comps[fbuf->reader_id].intercomm, &nranks);
//            else
//                MPI_Comm_remote_size(gl_comps[fbuf->writer_id].intercomm, &nranks);
            MPI_Comm_size(MPI_COMM_WORLD, &nranks);
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


int init_req_match_masters()
{
    if(gl_conf.distr_mode == FARB_UNDEFINED){
        gl_conf.nmasters = 0;
        gl_conf.masters = NULL;
        gl_conf.my_workgroup_sz = 0;
    } else {

        int wg, nranks, myrank, i;
        char* s = getenv("MAX_WORKGROUP_SIZE");

        if(s == NULL)
            wg = MAX_WORKGROUP_SIZE;
        else
            wg = atoi(s);
        assert(wg > 0);

        /* First, find out my master and, if I am a master, find out the size of my
           workgroup. Create the list of masters. */
        MPI_Comm_size(MPI_COMM_WORLD, &nranks);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        if(nranks <= wg){
            gl_conf.my_master = 0;
            gl_conf.my_workgroup_sz = nranks;
            gl_conf.nmasters = 1;
        } else {
            gl_conf.my_master = (int)(myrank/wg) * wg;
            gl_conf.my_workgroup_sz = wg;
            gl_conf.nmasters = (int)(nranks/wg);
            if(gl_conf.nmasters == 0)
                gl_conf.nmasters++;
            else if(nranks % wg > (int)(wg/2) ){
                if(myrank >= gl_conf.nmasters * wg){
                    gl_conf.my_master = gl_conf.nmasters * wg;
                    gl_conf.my_workgroup_sz = nranks % wg;
                }
                gl_conf.nmasters++;
            } else if ( (nranks % wg > 0) && (myrank >= (gl_conf.nmasters-1)*wg)){
                /*Merge last smaller group with the previous group*/
                gl_conf.my_master = (gl_conf.nmasters-1) * wg;
                gl_conf.my_workgroup_sz = wg + nranks % wg;
            }
        }
        FARB_DBG(VERBOSE_ALL_LEVEL, "My master %d, my wg size %d", gl_conf.my_master, gl_conf.my_workgroup_sz);
        if(gl_my_rank == 0)
            FARB_DBG(VERBOSE_ALL_LEVEL, "Nmasters %d", gl_conf.nmasters);

        gl_conf.masters = (int*)malloc(gl_conf.nmasters * sizeof(int));
        assert(gl_conf.masters != NULL);
        gl_conf.masters[0] = 0;
        for(i = 1; i < gl_conf.nmasters; i++){
            gl_conf.masters[i] = gl_conf.masters[i-1] + wg;
            FARB_DBG(VERBOSE_ALL_LEVEL, "Rank %d is a master", gl_conf.masters[i]);
        }

//        /*For each component, find out to which master I should send read requests*/
//        for(i=0; i < gl_ncomp; i++){
//            if(i == gl_my_comp_id)
//                continue;
//            MPI_Comm_remote_size(gl_comps[i].intercomm, &nrranks);
//
//            if(myrank < nrranks)
//                gl_comps[i].master = gl_conf.my_master; //use same rank
//            else {
//                int nmasters = (int)(nrranks/wg);
//                if( nrranks % wg > (int)(wg/2))
//                   nmasters++;
//                gl_comps[i].master = (nmasters-1)*wg;
//            }
//        }
    }
    return 0;
}

