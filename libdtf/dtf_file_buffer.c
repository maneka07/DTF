#include <assert.h>
#include <string.h>
#include <unistd.h>
#include "dtf_util.h"
#include "dtf_req_match.h"
#include "dtf_file_buffer.h"
#include "dtf.h"

void unpack_file_info(MPI_Offset bufsz, void *buf, file_buffer_t *fbf)
{
    int i, varid, nvars;
    file_buffer_t *fbuf;
    int ndims;
    dtf_var_t *var;
    MPI_Offset offt = 0;
    int type;
    MPI_Datatype dtype;
    MPI_Offset *shape;
    unsigned char *chbuf = (unsigned char*)buf;
    char filename[MAX_FILE_NAME];
    DTF_DBG(VERBOSE_DBG_LEVEL, "Start unpackinf file info for of sz %d", (int)bufsz);
    /*Unpack:
       - file name
       - how many ranks created the file
       - file header size
       - header
       - number of masters
       - master list
       - number of vars
       - vars
    */
	
	if(fbf == NULL){
		/*filename*/
		memcpy(filename, chbuf, MAX_FILE_NAME);
		offt += MAX_FILE_NAME ;
		DTF_DBG(VERBOSE_DBG_LEVEL, "unpack filename %s", filename);
		fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
		assert(fbuf != NULL);
	} else {
		fbuf = fbf;
		offt += MAX_FILE_NAME ;
	}
    /*ncid*/
    fbuf->ncid = -1;
    /*how many ranks created the file*/
    fbuf->cpl_mst_info->comm_sz = (int)(*((MPI_Offset*)(chbuf+offt)));
    offt += sizeof(MPI_Offset);
    assert(fbuf->cpl_mst_info->comm_sz > 0);
    /*header size*/
    fbuf->hdr_sz = *((MPI_Offset*)(chbuf+offt));
    offt += sizeof(MPI_Offset);
    /*header*/
    fbuf->header = dtf_malloc(fbuf->hdr_sz);
    memcpy(fbuf->header, chbuf+offt, fbuf->hdr_sz);
    offt += fbuf->hdr_sz + fbuf->hdr_sz%sizeof(MPI_Offset);
	/*number of masters*/
	fbuf->cpl_mst_info->nmasters = (int)(*((MPI_Offset*)(chbuf+offt)));
	DTF_DBG(VERBOSE_DBG_LEVEL, "unpack %d masters", fbuf->cpl_mst_info->nmasters);
	offt += sizeof(MPI_Offset);
	/*list of masters*/
	fbuf->cpl_mst_info->masters = dtf_malloc(fbuf->cpl_mst_info->nmasters*sizeof(int));
	memcpy(fbuf->cpl_mst_info->masters, chbuf+offt, fbuf->cpl_mst_info->nmasters*sizeof(int));
	offt += fbuf->cpl_mst_info->nmasters*sizeof(MPI_Offset);
	//init root master
	fbuf->root_writer = fbuf->cpl_mst_info->masters[0];
	
   
    /*number of vars*/
    nvars = (int)(*((MPI_Offset*)(chbuf+offt)));
    offt += sizeof(MPI_Offset);
    DTF_DBG(VERBOSE_DBG_LEVEL, "unpack nvars %d", nvars);
    /*vars*/
    for(i = 0; i < nvars; i++){
        varid = (int)(*((MPI_Offset*)(chbuf+offt)));
        offt += sizeof(MPI_Offset);
        assert(varid >= fbuf->nvars);

        type = (int)(*((MPI_Offset*)(chbuf+offt)));
        dtype = int2mpitype(type);
        offt+=sizeof(MPI_Offset);
        ndims = (int)(*((MPI_Offset*)(chbuf+offt)));
        offt += sizeof(MPI_Offset);
        shape = (MPI_Offset*)(chbuf+offt);
        offt += sizeof(MPI_Offset)*ndims;

        var = new_var(varid, ndims, dtype, shape);
        add_var(fbuf, var);
    }
    assert(offt == bufsz);
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished unpacking");
}

static int match_str(char* pattern, const char* filename)
{
    int ret = 0;
    char *next_token, *cur_token;
    int next_pos, cur_pos, token_len;
    char *match_str;
	
    if(strlen(pattern)== 0 || strlen(filename)==0)
        return ret;
        
	next_token = strchr(pattern, '%');
	
	if(next_token == NULL){
		//not a name pattern, just check for inclusion
		if(strstr(pattern, filename)!=NULL || strstr(filename, pattern)!=NULL)
			ret = 1;
	} else {
		cur_pos = 0;
		cur_token = &pattern[cur_pos];
		
		while( next_token != NULL){
			next_pos = next_token-pattern;
			token_len = next_token - cur_token;
			pattern[next_pos] = '\0';
			
			if(next_pos - cur_pos > 0){ //check for inclusion
				match_str = strstr(filename, cur_token);
				ret = (match_str == NULL) ? 0 : 1;
			}
			pattern[next_pos] = '%';
			
			if(!ret) break;
			
			cur_pos = next_pos+1;
			cur_token = &pattern[cur_pos];
			filename = match_str + token_len;
			next_token = strchr(cur_token, '%');
		}
		if( ret && (cur_pos != strlen(pattern)))
			//check the last token
			ret = (strstr(filename, cur_token) == NULL) ? 0 : 1;
	}
	
    return ret;
}

static void init_mst_info(master_info_t* mst_info, int comm_sz)
{
    mst_info->masters = NULL;
    mst_info->my_mst = -1;
    mst_info->my_wg_sz = 0;  
    mst_info->nmasters = 0;
    mst_info->iodb = NULL;
    mst_info->nread_completed = 0;
    mst_info->comm_sz = comm_sz;
}

static int init_req_match_masters(MPI_Comm comm, master_info_t *mst_info)
{
    int wg, nranks, myrank, i;
    char* s = getenv("MAX_WORKGROUP_SIZE");
    int my_master, my_master_glob;
    int *masters;

    DTF_DBG(VERBOSE_DBG_LEVEL, "Init masters");

    if(s == NULL)
        wg = MAX_WORKGROUP_SIZE;
    else
        wg = atoi(s);
    assert(wg > 0);

    /* First, find out my master and, if I am a master, find out the size of my
       workgroup. Create the list of masters. */
    MPI_Comm_size(comm, &nranks);
    MPI_Comm_rank(comm, &myrank);

    if(nranks <= wg){
        my_master = 0;
        mst_info->my_wg_sz = nranks;
        mst_info->nmasters = 1;
    } else {
        my_master = (int)(myrank/wg) * wg;
        mst_info->my_wg_sz = wg;
        mst_info->nmasters = (int)(nranks/wg);
        if(nranks % wg > 0){
            mst_info->nmasters++;
            if(myrank >= (mst_info->nmasters-1)*wg)
                mst_info->my_wg_sz = nranks % wg;
        }
    }

    if(myrank == 0)
        DTF_DBG(VERBOSE_DBG_LEVEL, "Nmasters %d", mst_info->nmasters);
    gl_proc.stats_info.nmasters = mst_info->nmasters;

    translate_ranks(&my_master, 1, comm, gl_proc.comps[gl_proc.my_comp].comm, &my_master_glob);
    mst_info->my_mst = my_master_glob;

    mst_info->masters = (int*)dtf_malloc(mst_info->nmasters * sizeof(int));
    assert(mst_info->masters != NULL);
    masters = (int*)dtf_malloc(mst_info->nmasters * sizeof(int));
    
    masters[0] = 0;
    for(i = 1; i < mst_info->nmasters; i++){
        masters[i] = masters[i-1] + wg;
    }    
    translate_ranks(masters, mst_info->nmasters, comm, gl_proc.comps[gl_proc.my_comp].comm, mst_info->masters);

    if(myrank == 0){
        for(i = 0; i < mst_info->nmasters; i++)
            DTF_DBG(VERBOSE_ALL_LEVEL, "Rank %d is a master", mst_info->masters[i]);
    }

    dtf_free(masters, mst_info->nmasters * sizeof(int));
	DTF_DBG(VERBOSE_DBG_LEVEL, "My wg size %d", mst_info->my_wg_sz);    
    mst_info->nread_completed = 0;
    return 0;
}

static void init_iodb(ioreq_db_t *iodb)
{
    iodb->ritems = NULL;
    iodb->witems = NULL;
    iodb->nritems = 0;
    iodb->updated_flag = 0;
}

/***************************File name pattern structure**************************************/

fname_pattern_t* new_fname_pattern()
{
    fname_pattern_t *pat = dtf_malloc(sizeof(fname_pattern_t));
	pat->fname[0]='\0';
    pat->iomode = DTF_UNDEFINED;
    pat->excl_fnames = NULL;
    pat->nexcls = 0;
    pat->next = NULL;
    pat->comp1 = DTF_UNDEFINED;
    pat->comp2 = DTF_UNDEFINED;
    pat->replay_io = 0;
    pat->rdr_recorded = DTF_UNDEFINED;
    pat->wrt_recorded = DTF_UNDEFINED;
    pat->write_only = 0;
    pat->mirror_io_root=0;
    pat->io_pats = NULL;
    pat->finfo = NULL;
    pat->finfo_sz = 0;
    pat->num_sessions = 1; //default we assume file will be used only once
    
    gl_proc.stats_info.num_fpats++;
    return pat;
}


fname_pattern_t *find_fname_pattern(const char *filename)
{
	int i;
	fname_pattern_t *pat;
	int found = 0;
	
	if (strlen (filename) == 0) return NULL;
		
    pat = gl_proc.fname_ptrns;
	while(pat != NULL){
		
		if(match_str(pat->fname, filename)){
			found = 1;
		 /*First see if the filename matches against any of the exclude patterns*/
			for(i = 0; i < pat->nexcls; i++)
				if(match_str(pat->excl_fnames[i], filename)){
					found = 0;
					break;					
				}
		}
		if(found)
			break;
		pat = pat->next;
	}
	return pat;
}


/***************************File structure**************************************/

file_buffer_t* find_file_buffer(file_buffer_t* buflist, const char* file_path, int ncid)
{
	/* If a file buffer for this file exists, return it. Otherwise, check if file name
	 * matches any file name pattern. If yes, create a new file buffer and return it.
	 * Otherwise, return NULL*/
    struct file_buffer *ptr = buflist;

    while(ptr != NULL)
    {
       if((ncid >= 0) && (ptr->ncid == ncid))
                break;
       if( (file_path != NULL) &&
           (strlen(file_path)!=0 && ( (strstr(file_path, ptr->file_path)!=NULL || strstr(ptr->file_path, file_path)!=NULL)))){
		   //DTF_DBG(VERBOSE_DBG_LEVEL, "Found by name in list %s", ptr->file_path);
           break;
	   }

       ptr = ptr->next;
    }

	//DTF_DBG(VERBOSE_DBG_LEVEL, "Return fbuf %p", (void*)ptr);
    return ptr;
}

void delete_file_buffer(file_buffer_t* fbuf)
{
    int i;
    int nvars=fbuf->nvars;

    if(fbuf == NULL)
        return;
	
    assert(gl_proc.filebuf_list != NULL);
    
    DTF_DBG(VERBOSE_DBG_LEVEL, "Will delete file buffer for %s", fbuf->file_path);

    if(fbuf->header != NULL)
        dtf_free(fbuf->header, fbuf->hdr_sz);

    if(fbuf->my_mst_info != NULL){

        if(fbuf->my_mst_info->iodb != NULL){
			double t_start = MPI_Wtime();
            clean_iodb(fbuf->my_mst_info->iodb, fbuf->nvars, fbuf->cpl_mst_info->comm_sz);
            dtf_free(fbuf->my_mst_info->iodb, sizeof(ioreq_db_t));
            gl_proc.stats_info.master_time += MPI_Wtime() - t_start;
        }

        if(fbuf->my_mst_info->masters != NULL)
            dtf_free(fbuf->my_mst_info->masters, fbuf->my_mst_info->nmasters*sizeof(int));
 
        dtf_free(fbuf->my_mst_info, sizeof(master_info_t));
    }

    if(fbuf->cpl_mst_info != NULL){
        if(fbuf->cpl_mst_info->masters != NULL)
            dtf_free(fbuf->cpl_mst_info->masters, fbuf->cpl_mst_info->nmasters*sizeof(int));
        dtf_free(fbuf->cpl_mst_info, sizeof(master_info_t));
    }

    delete_ioreqs(fbuf); 

    for(i = 0; i < nvars; i++)
        delete_var(fbuf, fbuf->vars[i]);
    dtf_free(fbuf->vars, nvars*sizeof(dtf_var_t*));

	if(fbuf->ioreq_log != NULL){
		io_req_log_t *ior = fbuf->ioreq_log;
		while(ior != NULL){
			fbuf->ioreq_log = fbuf->ioreq_log->next;
			if(ior->start != NULL) dtf_free(ior->start, ior->ndims*sizeof(MPI_Offset));
			if(ior->count != NULL)dtf_free(ior->count, ior->ndims*sizeof(MPI_Offset));
			dtf_free(ior, sizeof(io_req_log_t));
			ior = fbuf->ioreq_log;
		}
	} 

    assert(fbuf->rreq_cnt == 0);
    assert(fbuf->wreq_cnt == 0);

    if(gl_proc.filebuf_list == fbuf)
		gl_proc.filebuf_list = fbuf->next;

	if(fbuf->prev != NULL)
		fbuf->prev->next = fbuf->next;
	if(fbuf->next != NULL)
		fbuf->next->prev = fbuf->prev;

    dtf_free(fbuf, sizeof(file_buffer_t));
    
    {
		DTF_DBG(VERBOSE_DBG_LEVEL, "Deleted fbuf. Left:");
		fbuf = gl_proc.filebuf_list;
		while(fbuf != NULL){
			DTF_DBG(VERBOSE_DBG_LEVEL, "%s", fbuf->file_path);
			fbuf = fbuf->next;
		}
	}
    gl_proc.stats_info.nfiles--;
}


file_buffer_t *create_file_buffer(fname_pattern_t *pat, const char* file_path, MPI_Comm comm)
{

	file_buffer_t *buf;
	int comm_sz;
	
	DTF_DBG(VERBOSE_ALL_LEVEL, "Create file buffer for file %s", file_path);

	if(strlen(file_path) > MAX_FILE_NAME){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: filename %s longer than MAX_FILE_NAME (%d)", file_path, MAX_FILE_NAME);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
	}
	buf = dtf_malloc(sizeof(struct file_buffer));
    buf->next = NULL;
    buf->prev = NULL;
    buf->is_ready = 0;
    buf->vars = NULL; //RBTreeCreate(var_cmp, var_destroy, NullFunction, var_print, NullFunction);
    buf->nvars = 0;
    buf->header = NULL;
    buf->hdr_sz = 0;
    buf->rreq_cnt = 0;
    buf->wreq_cnt = 0;
    buf->ioreq_log = NULL;
    buf->session_cnt = 0;
    buf->ncid = -1;
    buf->done_matching_flag = 0;
    buf->done_multiple_flag = 0;
    buf->fready_notify_flag = DTF_UNDEFINED;
    buf->comm = MPI_COMM_NULL;
    buf->cpl_info_shared = 0;
    buf->is_write_only = pat->write_only;
    MPI_Comm_size(comm, &comm_sz);
    buf->my_mst_info = dtf_malloc(sizeof(master_info_t));
    assert(buf->my_mst_info != NULL);
    init_mst_info(buf->my_mst_info, comm_sz);    
	buf->my_mst_info->iodb = dtf_malloc(sizeof(struct ioreq_db));
	init_iodb(buf->my_mst_info->iodb);
            
    buf->cpl_mst_info = dtf_malloc(sizeof(master_info_t));
    init_mst_info(buf->cpl_mst_info, 0);
	strcpy(buf->file_path, file_path);
	buf->reader_id = -1;
	buf->writer_id = -1;
	buf->root_writer = -1;
    buf->root_reader = -1;
    buf->is_defined = 0;
	buf->iomode = pat->iomode;
	buf->cur_transfer_epoch = 0;
	buf->is_transferring = 0;
	buf->comm = comm;
	//insert
	if(gl_proc.filebuf_list == NULL)
		gl_proc.filebuf_list = buf;
	else{
		buf->next = gl_proc.filebuf_list;
		buf->next->prev = buf;
		gl_proc.filebuf_list = buf;
	}
	
	init_req_match_masters(comm, buf->my_mst_info);
	
	gl_proc.stats_info.nfiles++;
	return buf;
}

void finalize_files()
{
    int file_cnt = 0;
    file_buffer_t *fbuf = gl_proc.filebuf_list;

	DTF_DBG(VERBOSE_DBG_LEVEL, "Finalize files");
    while(fbuf != NULL){
        DTF_DBG(VERBOSE_DBG_LEVEL, "File %s, fready_notif_flag %d", fbuf->file_path,  fbuf->fready_notify_flag);
        if(fbuf->iomode == DTF_IO_MODE_FILE){
			 if(fbuf->root_writer == gl_proc.myrank && fbuf->fready_notify_flag == DTF_UNDEFINED)
				file_cnt++;
        } 

        fbuf = fbuf->next;
    }

    DTF_DBG(VERBOSE_DBG_LEVEL, "Process has to finalize notifications for %d files (out of %d files)", file_cnt, gl_proc.stats_info.nfiles);
    assert(file_cnt <= gl_proc.stats_info.nfiles);

	file_cnt = 0;
	fbuf = gl_proc.filebuf_list;
	while(fbuf != NULL){

		if(fbuf->iomode == DTF_IO_MODE_FILE){

			//~ if(fbuf->writer_id == gl_proc.my_comp && fbuf->fready_notify_flag == RDR_NOT_NOTIFIED){
				//~ while(fbuf->root_reader == -1)
					//~ progress_comm(0);
				//~ notify_file_ready(fbuf);
			//~ }
			if(fbuf->root_writer == gl_proc.myrank && fbuf->fready_notify_flag == DTF_UNDEFINED)
				while(fbuf->fready_notify_flag != RDR_NOTIFIED)
					progress_comm(0);
			 
		} 
		file_cnt++;
		fbuf = fbuf->next;
	}
	if(file_cnt != gl_proc.stats_info.nfiles){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "Ooops, something wrong with file counting: cnt %d vs gl_proc.stats_info.nfiles %d", file_cnt, gl_proc.stats_info.nfiles);
		fbuf = gl_proc.filebuf_list;
		while(fbuf != NULL){
			DTF_DBG(VERBOSE_ERROR_LEVEL, "%s", fbuf->file_path);
			fbuf = fbuf->next;
		}
	}
	assert(file_cnt == gl_proc.stats_info.nfiles);

    //MPI_Barrier(gl_proc.comps[gl_proc.my_comp].comm);
    DTF_DBG(VERBOSE_DBG_LEVEL, "Finished finalizing notifications. Will delete file bufs");
    /*Now, delete all file buffers*/
    fbuf = gl_proc.filebuf_list;
    while(fbuf != NULL){
        delete_file_buffer(fbuf);
        fbuf = gl_proc.filebuf_list;
    }
    gl_proc.filebuf_list = NULL;
}

void clean_iodb(ioreq_db_t *iodb, int nvars, int cpl_comm_sz)
{
	int i;
	write_db_item_t *witem;

    if(iodb == NULL){
		return;
	}

	if(iodb->witems != NULL){
		for(i = 0; i < nvars; i++){

			if(iodb->witems[i] == NULL)
				continue;

			witem = iodb->witems[i];

			if(witem->ndims > 0){
				RBTreeDestroy(witem->dblocks);
				gl_proc.stats_info.malloc_size -= witem->nblocks*(sizeof(block_t)+sizeof(MPI_Offset)*2*witem->ndims);
			} else
				dtf_free(witem->dblocks, sizeof(block_t));

			dtf_free(witem, sizeof(write_db_item_t));
			iodb->witems[i] = NULL;
		}
		
        dtf_free(iodb->witems, nvars*sizeof(write_db_item_t*));
        iodb->witems = NULL;
	}
	
	assert( iodb->ritems == NULL);
			
	//~ ritem = iodb->ritems;
	//~ while(ritem != NULL){
		//~ read_dblock_t *block = ritem->dblocks;
		//~ iodb->ritems =  iodb->ritems->next;
		//~ while(block != NULL){
			//~ ritem->dblocks = ritem->dblocks->next;
			//~ dtf_free(block->start, block->ndims*sizeof(MPI_Offset));
			//~ dtf_free(block->count, block->ndims*sizeof(MPI_Offset));
			
			//~ dtf_free(block, sizeof(read_dblock_t));
			//~ block = ritem->dblocks;
			//~ ritem->nblocks--;
		//~ }
		//~ assert(ritem->nblocks == 0);
		//~ dtf_free(ritem, sizeof(read_db_item_t));
		//~ iodb->nritems--;
	//~ }
    iodb->updated_flag = 0;

	DTF_DBG(VERBOSE_DBG_LEVEL, "iodb clean");
}


void close_file(file_buffer_t *fbuf)
{
	
	//~ if(fbuf->is_transferring){
		//~ if(gl_proc.conf.iodb_build_mode == IODB_BUILD_VARID)
			//~ send_ioreqs_by_var(fbuf);
		//~ else //if(gl_proc.conf.iodb_build_mode == IODB_BUILD_BLOCK)
			//~ send_ioreqs_by_block(fbuf);
	//~ }
	
    if(fbuf->writer_id == gl_proc.my_comp){

		fbuf->is_ready = 1;

        if(fbuf->iomode == DTF_IO_MODE_FILE) {
			//Check for any incoming messages
			progress_comm(1);
            if(fbuf->fready_notify_flag == RDR_NOT_NOTIFIED){
				assert(fbuf->root_writer == gl_proc.myrank);
				while(fbuf->root_reader == -1)
					progress_comm(0);
				notify_file_ready(fbuf);
			}
		} 
       
    } 

    fbuf->is_ready = 0;  //reset flag
}

void open_file(file_buffer_t *fbuf, MPI_Comm comm)
{
    DTF_DBG(VERBOSE_DBG_LEVEL,   "Enter dtf_open %s", fbuf->file_path);

    MPI_Status status;
    int rank; //, notif_open=1;
    int err;
	
    MPI_Comm_rank(comm, &rank);
    
    if(fbuf->iomode == DTF_IO_MODE_FILE) fbuf->is_defined = 1;

    if(fbuf->reader_id == gl_proc.my_comp){
        
        if(fbuf->iomode == DTF_IO_MODE_FILE){
			DTF_DBG(VERBOSE_DBG_LEVEL, "time_stamp open file %s", fbuf->file_path);

			if(rank == 0){
				if(fbuf->root_writer == -1){
					fbuf->root_writer = inquire_root(fbuf->file_path);
					send_mst_info(fbuf, fbuf->root_writer, fbuf->writer_id);
				}
				
				DTF_DBG(VERBOSE_DBG_LEVEL, "Waiting for file to become ready");
				/*root reader rank will wait until writer finishes writing the file.
				 then it will broadcast that the file is ready to everyone else*/

				while(!fbuf->is_ready)
					progress_comm(0);
			}
			MPI_Barrier(comm);
            fbuf->is_ready = 1;			
			
			DTF_DBG(VERBOSE_DBG_LEVEL, "time_stamp file ready %s", fbuf->file_path);

        } else if(fbuf->iomode == DTF_IO_MODE_MEMORY){
			void *buf;
			int bufsz;	

            if(fbuf->header == NULL){
				fname_pattern_t *pat;
				/*Zero rank will inquire the pnetcdf header/dtf vars info/masters info
				from writer's global zero rank and then broadcast this info to other
				readers that opened the file*/
				if(rank == 0){
					if(fbuf->root_writer == -1)
						fbuf->root_writer = inquire_root(fbuf->file_path);					

					
					send_mst_info(fbuf, fbuf->root_writer, fbuf->writer_id);
					
					DTF_DBG(VERBOSE_DBG_LEVEL, "Starting to wait for file info for %s", fbuf->file_path);
					err = MPI_Probe(fbuf->root_writer, FILE_INFO_TAG, gl_proc.comps[fbuf->writer_id].comm, &status);
					CHECK_MPI(err);
					MPI_Get_count(&status, MPI_BYTE, &bufsz);

					buf = dtf_malloc(bufsz);
					err = MPI_Recv(buf, bufsz, MPI_BYTE, fbuf->root_writer, FILE_INFO_TAG, gl_proc.comps[fbuf->writer_id].comm, &status);
					CHECK_MPI(err);
				}
				DTF_DBG(VERBOSE_DBG_LEVEL, "Bcast file info to others");
				err = MPI_Bcast(&bufsz, 1, MPI_INT, 0, comm);
				CHECK_MPI(err);
				assert(bufsz > 0);

				if(rank != 0) buf = dtf_malloc(bufsz);
				
				err = MPI_Bcast(buf, bufsz, MPI_BYTE, 0, comm);
				CHECK_MPI(err);

				unpack_file_info(bufsz, buf, fbuf);
				
				pat = find_fname_pattern(fbuf->file_path);
				assert(pat != NULL);
				if(pat->replay_io){
					assert(pat->finfo_sz == 0);
					pat->finfo_sz = bufsz;
					pat->finfo = buf;
				} else 				
					dtf_free(buf, bufsz);
				
			}

			fbuf->is_ready = 1;
			fbuf->is_defined = 1;
        }
        
		fbuf->cpl_info_shared = 1;
    
    } else if(fbuf->writer_id == gl_proc.my_comp){

		if(fbuf->iomode == DTF_IO_MODE_FILE && fbuf->root_writer == gl_proc.myrank)
			fbuf->fready_notify_flag = RDR_NOT_NOTIFIED;
		assert(fbuf->is_defined);
    }
    
    DTF_DBG(VERBOSE_DBG_LEVEL,   "Exit dtf_open %s", fbuf->file_path);
}

/*Write pnetcdf header*/
void write_hdr(file_buffer_t *fbuf, MPI_Offset hdr_sz, void *header)
{
    DTF_DBG(VERBOSE_DBG_LEVEL, "Writing header (sz %d)", (int)hdr_sz);
    fbuf->hdr_sz = hdr_sz;
    fbuf->header = dtf_malloc(hdr_sz);
    memcpy(fbuf->header, header, (size_t)hdr_sz);
    return;
}

/*Read pnetcdf header*/
MPI_Offset read_hdr_chunk(file_buffer_t *fbuf, MPI_Offset offset, MPI_Offset chunk_sz, void *chunk)
{
    if(offset+chunk_sz > fbuf->hdr_sz){
        DTF_DBG(VERBOSE_ALL_LEVEL, "Warning: trying to read %llu at offt %llu but hdr sz is %llu", chunk_sz, offset, fbuf->hdr_sz);
        chunk_sz = fbuf->hdr_sz - offset;
    }

    memcpy(chunk, (unsigned char*)fbuf->header+offset, (size_t)chunk_sz);
    return chunk_sz;
}

void pack_file_info(file_buffer_t *fbuf, MPI_Offset *bufsz, void **buf)
{
    dtf_var_t *var;
    MPI_Offset sz = 0, offt = 0;
    unsigned char *chbuf;
    int i;
 //   rb_red_blk_node *var_node;
    /*Pack:
       - file name
       - how many ranks opened it
       - file header size
       - header
       - number of masters
       - master list
       - number of vars
       - vars
    */

    sz =   MAX_FILE_NAME + fbuf->hdr_sz +
           fbuf->my_mst_info->nmasters*sizeof(MPI_Offset) +
           sizeof(MPI_Offset)*4;

    sz += sz%sizeof(MPI_Offset); //padding

    for(i = 0; i < fbuf->nvars; i++){
//        tmp_var.id = i;
//        var_node = RBExactQuery(fbuf->vars, &tmp_var);
//        assert(var_node != NULL);
//        var = (dtf_var_t*)(var_node->key);
        /*sz += sizeof(var->id) + sizeof(var->el_sz) +
        sizeof(var->ndims) + sizeof(MPI_Offset)*var->ndims;*/
        sz += sizeof(MPI_Offset)*3+ sizeof(MPI_Offset)*fbuf->vars[i]->ndims;
    }

    DTF_DBG(VERBOSE_ALL_LEVEL, "Packing info: sz %lld", sz);
    *buf = dtf_malloc(sz);
    chbuf = (unsigned char*)(*buf);

    /*filename*/
    memcpy(chbuf, fbuf->file_path, MAX_FILE_NAME);
    offt += MAX_FILE_NAME;
    /*how many ranks created the file*/
    *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)fbuf->my_mst_info->comm_sz;
    offt += sizeof(MPI_Offset);
    /*header size*/
    *((MPI_Offset*)(chbuf+offt)) = fbuf->hdr_sz;
    offt += sizeof(MPI_Offset);
    /*header*/
    memcpy(chbuf+offt, fbuf->header, fbuf->hdr_sz);
    offt += fbuf->hdr_sz + fbuf->hdr_sz%sizeof(MPI_Offset);
    /*number of masters*/
    *((MPI_Offset*)(chbuf+offt)) = fbuf->my_mst_info->nmasters;
    offt += sizeof(MPI_Offset);
    if(fbuf->my_mst_info->nmasters){
        DTF_DBG(VERBOSE_DBG_LEVEL, "pack %d masters", fbuf->my_mst_info->nmasters);
        if(gl_proc.conf.single_mpirun_mode){
			MPI_Comm intercomm = fbuf->reader_id == gl_proc.my_comp ? gl_proc.comps[fbuf->writer_id].comm : gl_proc.comps[fbuf->reader_id].comm;
			//need to translate ranks from local mpi_comm_world to global
			translate_ranks(fbuf->my_mst_info->masters,fbuf->my_mst_info->nmasters, 
			                gl_proc.comps[gl_proc.my_comp].comm, intercomm, (int*)(chbuf+offt));
		} else 
			/*list of masters*/
			memcpy(chbuf+offt, fbuf->my_mst_info->masters, fbuf->my_mst_info->nmasters*sizeof(int));
        offt += fbuf->my_mst_info->nmasters*sizeof(MPI_Offset); //sizeof(int) + padding for MPI_Offset
    }

    /*number of vars*/
    *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)fbuf->nvars;
    offt += sizeof(MPI_Offset);
    DTF_DBG(VERBOSE_DBG_LEVEL, "pack %d vars", fbuf->nvars);
    /*vars*/

    for(i = 0; i < fbuf->nvars; i++){
        var = fbuf->vars[i];
        *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)var->id;
        offt += sizeof(MPI_Offset);
        *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)mpitype2int(var->dtype);
        offt += sizeof(MPI_Offset);
        *((MPI_Offset*)(chbuf+offt)) = (MPI_Offset)var->ndims;
        offt += sizeof(MPI_Offset);
        memcpy((void*)(chbuf+offt), (void*)var->shape, sizeof(MPI_Offset)*var->ndims);
        offt += sizeof(MPI_Offset)*var->ndims;
    }
    DTF_DBG(VERBOSE_ALL_LEVEL, "offt %lld", offt);
    assert(offt == sz);
    *bufsz = sz;
}
