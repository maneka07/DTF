#include <string.h>
#include <ctype.h>

#include "dtf_config.h"
#include "dtf_util.h"
#include "dtf.h"
#include "dtf_io_pattern.h"
#include "dtf_component.h"

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
    fname_pattern_t *fb = gl_fname_ptrns;

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
        DTF_DBG(VERBOSE_DBG_LEVEL,   "File pattern %s, component 1 %s, component 2 %s ",
                fb->fname, gl_comps[fb->comp1].name, gl_comps[fb->comp2].name);
        switch(fb->iomode){
            case DTF_UNDEFINED:
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
	int ret = 0;
	int i;
    fname_pattern_t *cur_fpat = gl_fname_ptrns;

	for(i = 0; i < gl_ncomp; i++)
		if(gl_comps[i].name[0] == '\0'){
			DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error parsing config file: not all component names defined.");
			ret = 1;
		}

    while(cur_fpat!= NULL){
		if(cur_fpat->fname[0] == '\0'){
			DTF_DBG(VERBOSE_ERROR_LEVEL,"DTF Error parsing config file: file name not set");
			ret = 1;
		}

		if(cur_fpat->ignore_io){
			cur_fpat = cur_fpat->next;
			continue;
		}

		if(cur_fpat->comp1 == -1 || cur_fpat->comp2 == -1 ){
			DTF_DBG(VERBOSE_ERROR_LEVEL,"DTF Error parsing config file: file components not set for file %s", cur_fpat->fname);
			ret = 1;
		}

        if(cur_fpat->iomode == DTF_UNDEFINED){
            DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error parsing config file: I/O mode for file %s underfined", cur_fpat->fname);
            ret = 1;
        }
        
        if( (cur_fpat->iomode == DTF_IO_MODE_FILE) && cur_fpat->replay_io)
			cur_fpat->replay_io = 0; 

        cur_fpat = cur_fpat->next;
    }

    return ret;
}

void clean_config(){
  int i;
  fname_pattern_t *name_pat = gl_fname_ptrns;
  io_pattern_t *iopat;
  rank_pattern_t *rpat;
  void *tmp;
  
  while(name_pat != NULL){

	 for(i = 0; i < name_pat->nexcls; i++)
		dtf_free(name_pat->excl_fnames[i], sizeof(char)*MAX_FILE_NAME);
	 if(name_pat->nexcls > 0)
		dtf_free(name_pat->excl_fnames, sizeof(char*)*name_pat->nexcls);
	
	iopat = name_pat->io_pats;
	while(iopat != NULL){
		rpat = iopat->rank_pats;
		while(rpat != NULL){
			dtf_free(rpat->data, rpat->datasz);
		    tmp = rpat->next;
		    dtf_free(rpat, sizeof(rank_pattern_t));
		    rpat = (rank_pattern_t*)tmp;
		}
		tmp = iopat->next;
		dtf_free(iopat, sizeof(io_pattern_t));
		iopat = (io_pattern_t*)tmp;
	}

	 gl_fname_ptrns = gl_fname_ptrns->next;
	 dtf_free(name_pat, sizeof(fname_pattern_t));
	 name_pat = gl_fname_ptrns;
  }

  if(gl_comps != NULL)
	dtf_free(gl_comps, gl_ncomp*sizeof(component_t));

}

int load_config(const char *ini_name, const char *comp_name){

  FILE*      in;
  char       line[ASCIILINESZ+1];
  char       param[ASCIILINESZ], value[ASCIILINESZ];
  //char       comp_name[MAX_COMP_NAME];
  int        i, len, lineno=0;
  struct fname_pattern* cur_fpat = NULL;
  char *tmp_ininame;
 // ini = iniparser_load(ini_name);

  line[0] = 0;
  param[0] = 0;
  value[0] = 0;

  DTF_DBG(VERBOSE_DBG_LEVEL,   "Load config %s for comp %s", ini_name, comp_name);

  tmp_ininame = getenv("DTF_INI_FILE");
  if(tmp_ininame == NULL){
	if ((in=fopen(ini_name, "r"))==NULL) {
		DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: cannot open %s\n", ini_name);
		return 1 ;
	}
  } else {
	  	if ((in=fopen(tmp_ininame, "r"))==NULL) {
		DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error: cannot open %s\n", ini_name);
		return 1 ;
	}
  }


  gl_conf.log_ioreqs = 0;
  gl_conf.buffered_req_match = 0;
  gl_conf.do_checksum = 0;
  gl_conf.iodb_build_mode = IODB_BUILD_BLOCK; // default

  while(!feof(in)){

    fgets(line, ASCIILINESZ, in);
    lineno++;

    len = strlen(line);
    if(strlen(line) == 0 || line[0] == '#' || line[0] == ';' || line[0] == '!') //comment or empty line
        continue;

    if(line[len-1] != '\n' && !feof(in)){
        DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error parsing config file: Line %d in the ini file too long.", lineno);
        goto panic_exit;
    }

    strcpy(line, strip_str(line));
    strcpy(line, str_to_lower(line));

    if(strcmp(line, "[file]")==0){

        //check that we have already parsed the [info] section
        if (gl_ncomp == 0){
           DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error parsing config file: [file] section defined before [info] section.");
            goto panic_exit;
        }
        cur_fpat = new_fname_pattern();
        assert(cur_fpat != NULL);

        if(gl_fname_ptrns == NULL)
			gl_fname_ptrns = cur_fpat;
		else {
			cur_fpat->next = gl_fname_ptrns;
			gl_fname_ptrns = cur_fpat;
		}

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
                DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error parsing config file: invalid number of components.");
                goto panic_exit;
            }
            gl_comps = (struct component*)dtf_malloc(gl_ncomp*sizeof(struct component));
            assert(gl_comps!=NULL);
            for(i = 0; i<gl_ncomp; i++){
                gl_comps[i].id = i;
                gl_comps[i].connect_mode = DTF_UNDEFINED;
                gl_comps[i].in_msg_q = NULL;
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
        } else if(strcmp(param, "iodb_build_mode") == 0){
            if(strcmp(value, "varid") == 0)
               gl_conf.iodb_build_mode = IODB_BUILD_VARID;
            else if(strcmp(value, "range") == 0)
                gl_conf.iodb_build_mode = IODB_BUILD_BLOCK;
            else {
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error parsing config file: unknown iodb build mode (%s)", value);
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            }
        } else if(strcmp(param, "iodb_range") == 0){
			gl_conf.iodb_range = (MPI_Offset)atoi(value);
			assert(gl_conf.iodb_range > 0);
		} else if(strcmp(param, "buffered_req_match") == 0){

            gl_conf.buffered_req_match = atoi(value);
            if(gl_conf.buffered_req_match == 0)
                DTF_DBG(VERBOSE_DBG_LEVEL, "buffered_req_match disabled");
            else if(gl_conf.buffered_req_match == 1)
                DTF_DBG(VERBOSE_DBG_LEVEL, "buffered_req_match enabled");
            else {
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error parsing config file: Value for buffered_req_match should be 0 or 1.");
                goto panic_exit;
            }
        } else if(strcmp(param, "do_checksum") == 0){

            gl_conf.do_checksum = atoi(value);
            assert((gl_conf.do_checksum == 0) || (gl_conf.do_checksum == 1));

        } else if(strcmp(param, "log_ioreqs") == 0){

            gl_conf.log_ioreqs = atoi(value);
            assert((gl_conf.log_ioreqs == 0) || (gl_conf.log_ioreqs == 1));

        } else if(strcmp(param, "filename") == 0){
            assert(cur_fpat != NULL);
			assert(strlen(value) <= MAX_FILE_NAME);
			strcpy(cur_fpat->fname, value);

        } else if(strcmp(param, "exclude_name") == 0){
            assert(cur_fpat != NULL);
            assert(strlen(value) <= MAX_FILE_NAME);
            char **tmp = realloc(cur_fpat->excl_fnames, sizeof(char*)*(cur_fpat->nexcls+1));
            assert(tmp != NULL);
            gl_stats.malloc_size += sizeof(char*);
            cur_fpat->excl_fnames = tmp;
            cur_fpat->excl_fnames[cur_fpat->nexcls] = dtf_malloc(MAX_FILE_NAME*sizeof(char));
            assert(cur_fpat->excl_fnames[cur_fpat->nexcls] != NULL);
            strcpy(cur_fpat->excl_fnames[cur_fpat->nexcls], value);
            cur_fpat->nexcls++;

        } else if(strcmp(param, "comp1") == 0){
            assert(cur_fpat != NULL);
            assert(gl_ncomp != 0);
            for(i = 0; i < gl_ncomp; i++){
                if(strcmp(value, gl_comps[i].name) == 0){
                    cur_fpat->comp1 = i;
                    break;
                }
            }

        } else if(strcmp(param, "comp2") == 0){
            assert(cur_fpat != NULL);
            assert(gl_ncomp != 0);

            if(cur_fpat->comp2 != -1){
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error parsing config file: File %s cannot have multiple readers (%d)", cur_fpat->fname, cur_fpat->comp2);
                goto panic_exit;
            }
            for(i = 0; i < gl_ncomp; i++){
                if(strcmp(value, gl_comps[i].name) == 0){
                    cur_fpat->comp2 = i;
                    break;
                }
            }
		} else if(strcmp(param, "mode") == 0){
            assert(cur_fpat != NULL);

            if(strcmp(value, "file") == 0)
                cur_fpat->iomode = DTF_IO_MODE_FILE;
            else if(strcmp(value, "transfer") == 0){

                cur_fpat->iomode = DTF_IO_MODE_MEMORY;

                if(gl_msg_buf == NULL){
                    gl_msg_buf = dtf_malloc(gl_conf.data_msg_size_limit);
                    assert(gl_msg_buf != NULL);
                }
            } else {
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error parsing config file: unknown I/O mode: %s.", value);
				goto panic_exit;
			}
		}  else if(strcmp(param, "replay_io") == 0){
			cur_fpat->replay_io = atoi(value);
			assert(cur_fpat->replay_io==0 ||  cur_fpat->replay_io==1);
			
        }  else if(strcmp(param, "write_only") == 0){
			cur_fpat->write_only = atoi(value);
			assert(cur_fpat->write_only==0 ||  cur_fpat->write_only==1);
			
        }  else if(strcmp(param, "ignore_io") == 0){
			cur_fpat->ignore_io = atoi(value);
			assert(cur_fpat->ignore_io==0 ||  cur_fpat->ignore_io==1);

		} else {
            DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error parsing config file: unknown parameter %s.", param);
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
    cur_fpat = gl_fname_ptrns;
    while(cur_fpat != NULL){
		if(cur_fpat->ignore_io){
			cur_fpat = cur_fpat->next;
			continue;
		}
        if(cur_fpat->comp1 == gl_my_comp_id){
            if(gl_comps[ cur_fpat->comp2 ].connect_mode == DTF_UNDEFINED )
               gl_comps[ cur_fpat->comp2 ].connect_mode = CONNECT_MODE_SERVER; // I am server

        } else if( gl_comps[ cur_fpat->comp1 ].connect_mode == DTF_UNDEFINED){
            if(cur_fpat->comp2 == gl_my_comp_id )
               gl_comps[ cur_fpat->comp1 ].connect_mode = CONNECT_MODE_CLIENT;  //I am client
        }

        cur_fpat = cur_fpat->next;
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

