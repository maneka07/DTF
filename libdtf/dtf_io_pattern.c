#include <stdlib.h>
#include "dtf_util.h"
#include "dtf_io_pattern.h"
#include "dtf_req_match.h"

void record_io_pat(char *filename, int rank, void *pat_data, size_t datasz,  int cur_transfer_epoch)
{
	fname_pattern_t *pat;
	io_pattern_t *iopat;
	rank_pattern_t *rpat;
	int created = 0;
	
	pat = find_fname_pattern(filename);
	assert(pat != NULL);
	
	if(!pat->replay_io)
		return;
		
	if(pat->wrt_recorded == IO_PATTERN_RECORDED)
		return;
	
	if(pat->io_pats == NULL){
		iopat = dtf_malloc(sizeof(struct io_pattern));
		assert(iopat != NULL);
		pat->io_pats = iopat;
		created = 1;
	} else {
		iopat = pat->io_pats;
		while(iopat->next != NULL){
			if(iopat->transfer_epoch == cur_transfer_epoch)
				break;
			iopat = iopat->next;
		}
		if(iopat->transfer_epoch != cur_transfer_epoch){
			io_pattern_t *tmp = dtf_malloc(sizeof(struct io_pattern));
			assert(tmp != NULL);
			iopat->next = tmp;
			iopat = iopat->next;
			created = 1;
		}
	}
	
	if(created){
		iopat->transfer_epoch = cur_transfer_epoch;
		iopat->rank_pats = NULL;
		iopat->next = NULL;
	}
	
	created = 0;
	if(iopat->rank_pats == NULL){
		rpat = dtf_malloc(sizeof(rank_pattern_t));
		assert(rpat != NULL);
		iopat->rank_pats = rpat;
		created = 1;
	} else {
		rpat = iopat->rank_pats;
		while(rpat->next != NULL){
			if(rpat->rank == rank)
				break;
			rpat = rpat->next;
		}
		if(rpat->rank != rank){
			rank_pattern_t *tmp = dtf_malloc(sizeof(rank_pattern_t));
			assert(tmp != NULL);
			rpat->next = tmp;
			rpat = rpat->next;
			created = 1;
		}
	}
	
	if(created){
		rpat->rank = rank;
		rpat->datasz = datasz;
		rpat->data = dtf_malloc(datasz);
		assert(rpat->data);
		memcpy(rpat->data, pat_data, datasz);
		rpat->next = NULL;
	} else {
		//merge
		size_t shift = sizeof(MPI_Offset) + MAX_FILE_NAME + MAX_FILE_NAME%sizeof(MPI_Offset);
		void *tmp = realloc(rpat->data, rpat->datasz + datasz - shift);
		assert(tmp != NULL);
		rpat->data = tmp;
		memcpy((unsigned char*)rpat->data + rpat->datasz, (unsigned char*)pat_data + shift, datasz - shift);
		rpat->datasz += datasz - shift;
		gl_stats.malloc_size += datasz - shift;
	}
}

void replay_io(fname_pattern_t *pat, char *filename, int epoch)
{
	io_pattern_t *iopat;
	rank_pattern_t *rpat;
	
	/*Find I/O pattern for the given epoch*/
	iopat = pat->io_pats;
	while(iopat != NULL){
		if(iopat->transfer_epoch == epoch)
			break;
		iopat = iopat->next;
	}
	
	if(iopat == NULL){
		DTF_DBG(VERBOSE_ERROR_LEVEL, "DTF Error: trying to replay I/O for epoch %d but I/O pattern for this epoch hasn't been recorded yet", epoch);
		assert(0);
	}
	
	rpat = iopat->rank_pats;
	while(rpat != NULL){
		DTF_DBG(VERBOSE_DBG_LEVEL, "Replay data for file %s for rank %d", filename, rpat->rank);
		//first replace the file name in the pattern data with the most recent filename
		memcpy( (unsigned char*)rpat->data+sizeof(MPI_Offset), filename, MAX_FILE_NAME);
		send_data(rpat->data, rpat->datasz);
		rpat = rpat->next;
	}
} 