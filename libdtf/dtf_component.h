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


#ifndef DTF_COMPONENT_H_INCLUDED
#define DTF_COMPONENT_H_INCLUDED

#define MAX_COMP_NAME 32

#include <mpi.h>

typedef struct component{
    unsigned int    id;
    char            name[MAX_COMP_NAME];
    int             connect_mode; /*0 - I am server, 1 - I am client, -1 - undefined (no interconnection)*/
    struct dtf_msg  *in_msg_q;   /*Queue of incoming messages (those that cannot be processed right now)*/
    struct dtf_msg  *out_msg_q;
    MPI_Comm        comm;   /*intra or inter component communicator*/
    int             finalized;   /*set to one when component calls finalize function*/
}component_t;

int 	init_comp_comm();
void 	finalize_comp_comm();

#endif // DTF_COMPONENT_H_INCLUDED
