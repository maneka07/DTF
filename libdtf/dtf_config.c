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
    fname_pattern_t *fb = gl_proc.fname_ptrns;

    if(gl_proc.comps == NULL) return;

    DTF_DBG(VERBOSE_DBG_LEVEL,   "Number of components: %d", gl_proc.ncomps);
    DTF_DBG(VERBOSE_DBG_LEVEL,   "I am comp %d, name %s", gl_proc.my_comp, gl_proc.comps[gl_proc.my_comp].name);
    for(i = 0; i<gl_proc.ncomps; i++){
        if(i != gl_proc.my_comp){
            if(gl_proc.comps[i].connect_mode == CONNECT_MODE_CLIENT)
                DTF_DBG(VERBOSE_DBG_LEVEL,   "Component %s is my server", gl_proc.comps[i].name);
            else if (gl_proc.comps[i].connect_mode == CONNECT_MODE_SERVER)
                DTF_DBG(VERBOSE_DBG_LEVEL,   "Component %s is my client", gl_proc.comps[i].name);
            else
                DTF_DBG(VERBOSE_DBG_LEVEL,   "No I/O with %s\n", gl_proc.comps[i].name);
        }
    }

    while(fb != NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL,   "File pattern %s, component 1 %s, component 2 %s ",
                fb->fname, gl_proc.comps[fb->comp1].name, gl_proc.comps[fb->comp2].name);
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
    fname_pattern_t *cur_fpat = gl_proc.fname_ptrns;

	for(i = 0; i < gl_proc.ncomps; i++)
		if(gl_proc.comps[i].name[0] == '\0'){
			DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error parsing config file: not all component names defined.");
			ret = 1;
		}

    while(cur_fpat!= NULL){
		if(cur_fpat->fname[0] == '\0'){
			DTF_DBG(VERBOSE_ERROR_LEVEL,"DTF Error parsing config file: file name not set");
			ret = 1;
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
  fname_pattern_t *name_pat = gl_proc.fname_ptrns;
  io_pattern_t *iopat;
  rank_pattern_t *rpat;
  void *tmp;
  
  while(name_pat != NULL){

	 for(i = 0; i < name_pat->nexcls; i++)
		dtf_free(name_pat->excl_fnames[i], sizeof(char)*MAX_FILE_NAME);
	 if(name_pat->nexcls > 0)
		dtf_free(name_pat->excl_fnames, sizeof(char*)*name_pat->nexcls);
	if(name_pat->finfo_sz > 0)
		dtf_free(name_pat->finfo, name_pat->finfo_sz);
		
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

	 gl_proc.fname_ptrns = gl_proc.fname_ptrns->next;
	 dtf_free(name_pat, sizeof(fname_pattern_t));
	 name_pat = gl_proc.fname_ptrns;
  }

  if(gl_proc.comps != NULL)
	dtf_free(gl_proc.comps, gl_proc.ncomps*sizeof(component_t));

}

static int parse_config(const char *ini_name, const char *comp_name){

  FILE*      in;
  char       line[ASCIILINESZ+1];
  char       param[ASCIILINESZ], value[ASCIILINESZ];
  //char       comp_name[MAX_COMP_NAME];
  int        i, len, lineno=0;
  fname_pattern_t* cur_fpat = NULL;
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
        if (gl_proc.ncomps == 0){
           DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error parsing config file: [file] section defined before [info] section.");
            goto panic_exit;
        }
        cur_fpat = new_fname_pattern();
        assert(cur_fpat != NULL);

        if(gl_proc.fname_ptrns == NULL)
			gl_proc.fname_ptrns = cur_fpat;
		else {
			cur_fpat->next = gl_proc.fname_ptrns;
			gl_proc.fname_ptrns = cur_fpat;
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
            gl_proc.ncomps = atoi(value);
            if(gl_proc.ncomps <= 0){
                DTF_DBG(VERBOSE_ERROR_LEVEL,   "DTF Error parsing config file: invalid number of components.");
                goto panic_exit;
            }
            gl_proc.comps = (struct component*)dtf_malloc(gl_proc.ncomps*sizeof(struct component));
            for(i = 0; i<gl_proc.ncomps; i++){
                gl_proc.comps[i].id = i;
                gl_proc.comps[i].connect_mode = DTF_UNDEFINED;
                gl_proc.comps[i].in_msg_q = NULL;
                gl_proc.comps[i].name[0] = 0;
                gl_proc.comps[i].comm = MPI_COMM_NULL;
                gl_proc.comps[i].finalized = 0;
                gl_proc.comps[i].out_msg_q = NULL;
            }

        } else if(strcmp(param, "comp_name") == 0){
            assert(gl_proc.ncomps!=0);
            assert(strlen(value) <= MAX_COMP_NAME);
            for(i = 0; i<gl_proc.ncomps; i++){
                if(gl_proc.comps[i].name[0] == 0){
                    strcpy(gl_proc.comps[i].name, value);
                    if(strcmp(gl_proc.comps[i].name, comp_name) == 0)
                        gl_proc.my_comp = i;
                    break;
                }
            }
        } else if(strcmp(param, "iodb_build_mode") == 0){
            if(strcmp(value, "varid") == 0)
               gl_proc.conf.iodb_build_mode = IODB_BUILD_VARID;
            else if(strcmp(value, "range") == 0)
                gl_proc.conf.iodb_build_mode = IODB_BUILD_BLOCK;
            else {
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error parsing config file: unknown iodb build mode (%s)", value);
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            }
        } else if(strcmp(param, "iodb_range") == 0){
			gl_proc.conf.iodb_range = (MPI_Offset)atoi(value);
			assert(gl_proc.conf.iodb_range > 0);
		} else if(strcmp(param, "buffer_data") == 0){

            gl_proc.conf.buffer_data = atoi(value);
            if(gl_proc.conf.buffer_data == 0)
                DTF_DBG(VERBOSE_DBG_LEVEL, "buffer_data disabled");
            else if(gl_proc.conf.buffer_data == 1)
                DTF_DBG(VERBOSE_DBG_LEVEL, "buffer_data enabled");
            else {
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error parsing config file: Value for buffer_data should be 0 or 1.");
                goto panic_exit;
            }
        } else if(strcmp(param, "do_checksum") == 0){

            gl_proc.conf.do_checksum = atoi(value);
            assert((gl_proc.conf.do_checksum == 0) || (gl_proc.conf.do_checksum == 1));

        } else if(strcmp(param, "log_ioreqs") == 0){

            gl_proc.conf.log_ioreqs = atoi(value);
            assert((gl_proc.conf.log_ioreqs == 0) || (gl_proc.conf.log_ioreqs == 1));

        } else if(strcmp(param, "filename") == 0){
            assert(cur_fpat != NULL);
			assert(strlen(value) <= MAX_FILE_NAME);
			strcpy(cur_fpat->fname, value);

        } else if(strcmp(param, "exclude_name") == 0){
            assert(cur_fpat != NULL);
            assert(strlen(value) <= MAX_FILE_NAME);
            char **tmp = realloc(cur_fpat->excl_fnames, sizeof(char*)*(cur_fpat->nexcls+1));
            assert(tmp != NULL);
            gl_proc.stats_info.malloc_size += sizeof(char*);
            cur_fpat->excl_fnames = tmp;
            cur_fpat->excl_fnames[cur_fpat->nexcls] = dtf_malloc(MAX_FILE_NAME*sizeof(char));
            strcpy(cur_fpat->excl_fnames[cur_fpat->nexcls], value);
            cur_fpat->nexcls++;

        } else if(strcmp(param, "comp1") == 0){
            assert(cur_fpat != NULL);
            assert(gl_proc.ncomps != 0);
            for(i = 0; i < gl_proc.ncomps; i++){
                if(strcmp(value, gl_proc.comps[i].name) == 0){
                    cur_fpat->comp1 = i;
                    break;
                }
            }

        } else if(strcmp(param, "comp2") == 0){
            assert(cur_fpat != NULL);
            assert(gl_proc.ncomps != 0);

            if(cur_fpat->comp2 != -1){
                DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error parsing config file: File %s cannot have multiple readers (%d)", cur_fpat->fname, cur_fpat->comp2);
                goto panic_exit;
            }
            for(i = 0; i < gl_proc.ncomps; i++){
                if(strcmp(value, gl_proc.comps[i].name) == 0){
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
                if(gl_proc.conf.use_msg_buffer && gl_proc.msgbuf == NULL)
					gl_proc.msgbuf = dtf_malloc(gl_proc.conf.data_msg_size_limit);
            } else {
				DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error parsing config file: unknown I/O mode: %s.", value);
				goto panic_exit;
			}
		}  else if(strcmp(param, "replay_io") == 0){
			cur_fpat->replay_io = atoi(value);
			assert(cur_fpat->replay_io==0 ||  cur_fpat->replay_io==1);
			
        } else if(strcmp(param, "mirror_io_root") == 0){
			cur_fpat->mirror_io_root = atoi(value);
			assert(cur_fpat->mirror_io_root==0 ||  cur_fpat->mirror_io_root==1);
			
        } else if(strcmp(param, "write_only") == 0){
			cur_fpat->write_only = atoi(value);
			assert(cur_fpat->write_only==0 ||  cur_fpat->write_only==1);
			
        } else if(strcmp(param, "num_sessions") == 0){
			cur_fpat->num_sessions = atoi(value);
			assert(cur_fpat->num_sessions > 0);

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

//    for(i = 0; i < gl_proc.ncomps; i++)
//        if(strcmp(gl_proc.comps[i].name, comp_name) == 0)
//            gl_proc.my_comp = i;



/*
  During the parsing of the configuarion file, if a file direct data transfer will
  take place between two components, each component is automatically set as a server
  (CONNECT_MODE_SERVER) or client(CONNECT_MODE_CLIENT), respectively, in order to
  esablish an intercommunicatior between these two components later. If there is no
  file flow then the connection mode is not set(DTF_UNDEFINED).
*/
    cur_fpat = gl_proc.fname_ptrns;
    while(cur_fpat != NULL){
        if(cur_fpat->comp1 == gl_proc.my_comp){
            if(gl_proc.comps[ cur_fpat->comp2 ].connect_mode == DTF_UNDEFINED )
               gl_proc.comps[ cur_fpat->comp2 ].connect_mode = CONNECT_MODE_SERVER; // I am server

        } else if( gl_proc.comps[ cur_fpat->comp1 ].connect_mode == DTF_UNDEFINED){
            if(cur_fpat->comp2 == gl_proc.my_comp )
               gl_proc.comps[ cur_fpat->comp1 ].connect_mode = CONNECT_MODE_CLIENT;  //I am client
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

void pack_config(void **buf, int *offt1)
{
	int i;
	fname_pattern_t *pat;
	unsigned char *chbuf;
	
	int offt = 0, sz = 0;
	
	sz += 7*sizeof(int)+gl_proc.ncomps*MAX_COMP_NAME;
	pat = gl_proc.fname_ptrns;
	while(pat != NULL){
		sz += 8*sizeof(int)+ MAX_FILE_NAME*(pat->nexcls+1);
		pat = pat->next;
	}
	
	*buf = dtf_malloc(sz);
	chbuf = (unsigned char*)*buf;
	
    *(int*)(chbuf) = gl_proc.ncomps;
    offt+=sizeof(int);
    
    for(i = 0; i < gl_proc.ncomps; i++){
		strcpy((char*)(chbuf+offt), gl_proc.comps[i].name);
		offt += MAX_COMP_NAME; 
	}
	
	*(int*)(chbuf+offt) = gl_proc.conf.buffer_data;
	offt += sizeof(int);
	*(int*)(chbuf+offt) = gl_proc.conf.iodb_build_mode;
	offt += sizeof(int);
	*(int*)(chbuf+offt) = (int)gl_proc.conf.iodb_range;
	offt += sizeof(int);
	*(int*)(chbuf+offt) = gl_proc.conf.do_checksum;
	offt += sizeof(int);
	*(int*)(chbuf+offt) = gl_proc.conf.log_ioreqs;
	offt += sizeof(int);
	*(int*)(chbuf+offt) = gl_proc.stats_info.num_fpats;
	offt += sizeof(int);
	
	pat = gl_proc.fname_ptrns;
	while(pat != NULL){
		
		strcpy((char*)(chbuf+offt), pat->fname);
		offt+=MAX_FILE_NAME;
		
		*(int*)(chbuf+offt) = pat->comp1;
		offt+=sizeof(int);
		*(int*)(chbuf+offt) = pat->comp2;
		offt+=sizeof(int);
		*(int*)(chbuf+offt) = pat->iomode;
		offt+=sizeof(int);
		*(int*)(chbuf+offt) = pat->num_sessions;
		offt+=sizeof(int);
		*(int*)(chbuf+offt) = pat->write_only;
		offt+=sizeof(int);
		*(int*)(chbuf+offt) = pat->mirror_io_root;
		offt+=sizeof(int);
		*(int*)(chbuf+offt) = pat->replay_io;
		offt+=sizeof(int);
		*(int*)(chbuf+offt) = pat->nexcls;
		offt+=sizeof(int);
		
		for(i=0; i<pat->nexcls;i++){
			strcpy((char*)(chbuf+offt), pat->excl_fnames[i]);
			offt += MAX_FILE_NAME;
		}
		
		pat = pat->next;
	}
	assert(offt == sz);
	*offt1 = offt;
}

void unpack_config(void *buf,const char* comp_name)
{
	int i, j;
	fname_pattern_t *pat;
	unsigned char *chbuf;
	
	int npats;
	
	int offt = 0;
	assert(buf != NULL);
	chbuf = (unsigned char*)buf;
	
    gl_proc.ncomps = *(int*)(chbuf);
    offt+=sizeof(int);
    
    gl_proc.comps = (struct component*)dtf_malloc(gl_proc.ncomps*sizeof(struct component));
	for(i = 0; i<gl_proc.ncomps; i++){
		gl_proc.comps[i].id = i;
		gl_proc.comps[i].connect_mode = DTF_UNDEFINED;
		gl_proc.comps[i].in_msg_q = NULL;
		
		strcpy(gl_proc.comps[i].name, (char*)(chbuf+offt));
		offt += MAX_COMP_NAME; 

		if(strcmp(gl_proc.comps[i].name, comp_name) == 0)
			gl_proc.my_comp = i;
                        
		gl_proc.comps[i].comm = MPI_COMM_NULL;
		gl_proc.comps[i].finalized = 0;
		gl_proc.comps[i].out_msg_q = NULL;
	}
	
	gl_proc.conf.buffer_data = *(int*)(chbuf+offt);
	offt += sizeof(int);
	gl_proc.conf.iodb_build_mode = *(int*)(chbuf+offt);
	offt += sizeof(int);
	gl_proc.conf.iodb_range = (MPI_Offset)(*(int*)(chbuf+offt));
	offt += sizeof(int);
	gl_proc.conf.do_checksum = *(int*)(chbuf+offt);
	offt += sizeof(int);
	gl_proc.conf.log_ioreqs = *(int*)(chbuf+offt);
	offt += sizeof(int);
	
	npats = *(int*)(chbuf+offt);
	offt += sizeof(int);
	
	for(i = 0; i < npats; i++){
		pat = new_fname_pattern();
        assert(pat != NULL);

        if(gl_proc.fname_ptrns == NULL)
			gl_proc.fname_ptrns = pat;
		else {
			pat->next = gl_proc.fname_ptrns;
			gl_proc.fname_ptrns = pat;
		}
		
		strcpy(pat->fname, (char*)(chbuf+offt));
		offt+=MAX_FILE_NAME;
		
		pat->comp1 = *(int*)(chbuf+offt);
		offt+=sizeof(int);
		pat->comp2 = *(int*)(chbuf+offt);
		offt+=sizeof(int);
		pat->iomode = *(int*)(chbuf+offt);
		offt+=sizeof(int);
		if(pat->iomode == DTF_IO_MODE_MEMORY && gl_proc.conf.use_msg_buffer && gl_proc.msgbuf == NULL)
			gl_proc.msgbuf = dtf_malloc(gl_proc.conf.data_msg_size_limit);
		pat->num_sessions = *(int*)(chbuf+offt);
		offt+=sizeof(int);
		pat->write_only = *(int*)(chbuf+offt);
		offt+=sizeof(int);
		pat->mirror_io_root = *(int*)(chbuf+offt);
		offt+=sizeof(int);
		pat->replay_io = *(int*)(chbuf+offt);
		offt+=sizeof(int);
		pat->nexcls = *(int*)(chbuf+offt);
		offt+=sizeof(int);
		
		pat->excl_fnames = dtf_malloc( sizeof(char*)*pat->nexcls);
		
		for(j=0; j<pat->nexcls;j++){
            pat->excl_fnames[j] = dtf_malloc(MAX_FILE_NAME*sizeof(char));
            strcpy(pat->excl_fnames[j],(char*)(chbuf+offt));
			offt += MAX_FILE_NAME;
		}
	}
	
	pat = gl_proc.fname_ptrns;
    while(pat != NULL){
        if(pat->comp1 == gl_proc.my_comp){
            if(gl_proc.comps[ pat->comp2 ].connect_mode == DTF_UNDEFINED )
               gl_proc.comps[ pat->comp2 ].connect_mode = CONNECT_MODE_SERVER; // I am server

        } else if( gl_proc.comps[ pat->comp1 ].connect_mode == DTF_UNDEFINED){
            if(pat->comp2 == gl_proc.my_comp )
               gl_proc.comps[ pat->comp1 ].connect_mode = CONNECT_MODE_CLIENT;  //I am client
        }

        pat = pat->next;
    }
}

int load_config(const char *dtf_ini_path, const char *comp_name)
{
	int offt = 0, err;
	void *buf = NULL;
	int nprocs;
	
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	if(gl_proc.myrank == 0){
		err = parse_config(dtf_ini_path, comp_name);
		if(err) MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	}
	
	if(nprocs > 1) { 
		if(gl_proc.myrank == 0)
			pack_config(&buf, &offt);
		 err = MPI_Bcast(&offt, 1, MPI_INT, 0, MPI_COMM_WORLD );
		 CHECK_MPI(err);
		 
		 if(gl_proc.myrank != 0)
			buf = dtf_malloc(offt);
			
		 err = MPI_Bcast(buf, offt, MPI_CHAR, 0, MPI_COMM_WORLD );
		 CHECK_MPI(err);
		 
		 if(gl_proc.myrank != 0)
			unpack_config(buf, comp_name);
			
		dtf_free(buf, offt);
	}
	DTF_DBG(VERBOSE_DBG_LEVEL, "Config loaded");
	return 0;
}
