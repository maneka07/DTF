#ifndef DTF_IO_PATTERN_H_INCLUDED
#define DTF_IO_PATTERN_H_INCLUDED

#define IO_PATTERN_RECORDING 0
#define IO_PATTERN_RECORDED  1
#include <stdio.h>

typedef struct rank_pattern{
	int rank;
	void *data;
	size_t datasz;
	
	struct rank_pattern *next;
}rank_pattern_t;

typedef struct io_pattern{
	int transfer_epoch;  /*In which transfer epoch was this pattern recorded?*/
	rank_pattern_t *rank_pats;
	struct io_pattern *next;
}io_pattern_t;


void record_io_pat(char *filename, int rank, void *pat_data, size_t datasz,  int cur_transfer_epoch);
void replay_io_pat(fname_pattern_t *pat, char *filename, int epoch);
#endif
