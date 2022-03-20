/*
 © Copyright 2019 RIKEN Center for Computational Science, System Software
 *	Development Team All rights reserved.

 * The Data Transfer Framework (DTF) was designed to work with 
 * the PnetCDF library developed by Northwestern University and 
 * Argonne National Laboratory. The modified version of PnetCDF is 
 * provided with this distribution.
 * 
 * Access and use of this software shall impose the following obligations 
 * and understandings on the user. The user is granted the right, without 
 * any fee or cost, to use, copy, modify, alter, enhance and distribute 
 * this software, and any derivative works thereof, and its supporting 
 * documentation for any purpose whatsoever, provided that this entire 
 * notice appears in all copies of the software, derivative works and 
 * supporting documentation.  Further, RIKEN requests that the user credit 
 * RIKEN in any publications that result from the use of this software or 
 * in any product that includes this software.  The name RIKEN, however, 
 * may not be used in any advertising or publicity to endorse or promote 
 * any products or commercial entity unless specific written permission is 
 * obtained from RIKEN. The user also understands that RIKEN is not 
 * obligated to provide the user with any support, consulting, training or 
 * assistance of any kind with regard to the use, operation and 
 * performance of this software nor to provide the user with any updates, 
 * revisions, new versions or "bug fixes."
 * 
 * THIS SOFTWARE IS PROVIDED BY RIKEN "AS IS" AND ANY EXPRESS OR IMPLIED 
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN 
 * NO EVENT SHALL RIKEN BE LIABLE FOR ANY SPECIAL, INDIRECT OR 
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF 
 * USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR 
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE ACCESS, 
 * USE OR PERFORMANCE OF THIS SOFTWARE.
 *
*/

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
	
	DTF_DBG(VERBOSE_DBG_LEVEL, "Record io pattern");
	
	if(pat->io_pats == NULL){
		iopat = dtf_malloc(sizeof(struct io_pattern));
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
			rpat->next = tmp;
			rpat = rpat->next;
			created = 1;
		}
	}
	
	if(created){
		rpat->rank = rank;
		rpat->datasz = datasz;
		rpat->data = dtf_malloc(datasz);
		memcpy(rpat->data, pat_data, datasz);
		rpat->next = NULL;
		DTF_DBG(VERBOSE_DBG_LEVEL, "Created rank pattern for %d", rank);
	} else {
		//merge
		size_t shift = sizeof(MPI_Offset); /*rdr_rank*/
		void *tmp = realloc(rpat->data, rpat->datasz + datasz - shift);
		assert(tmp != NULL);
		rpat->data = tmp;
		memcpy((unsigned char*)rpat->data + rpat->datasz, (unsigned char*)pat_data + shift, datasz - shift);
		rpat->datasz += datasz - shift;
		gl_proc.stats_info.malloc_size += datasz - shift;
		DTF_DBG(VERBOSE_DBG_LEVEL, "Appended pattern for rank %d", rank);
	}
}

void replay_io_pat(fname_pattern_t *pat, char *filename, int epoch)
{
	io_pattern_t *iopat;
	rank_pattern_t *rpat;
	file_buffer_t *fbuf;
	
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
	fbuf = find_file_buffer(gl_proc.filebuf_list, filename, -1);
	if(fbuf == NULL)
		assert(0); 
	rpat = iopat->rank_pats;
	while(rpat != NULL){
		DTF_DBG(VERBOSE_DBG_LEVEL, "Replay data for for rank %d", rpat->rank);
		if(gl_proc.conf.use_msg_buffer)
			send_data_blocking(fbuf, (unsigned char*)(rpat->data), rpat->datasz);
		else
			send_data_nonblocking(fbuf, (unsigned char*)(rpat->data), rpat->datasz);
		rpat = rpat->next;
	}
	while(gl_proc.comps[fbuf->reader_id].out_msg_q != NULL)
		progress_send_queue();
} 
