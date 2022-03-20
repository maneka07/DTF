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

#include "rb_misc.h"

/*  CONVENTIONS:  All data structures for stacks have the prefix */
/*                "stk_" to prevent name conflicts. */
/*                                                                      */
/*                Function names: Each word in a function name begins with */
/*                a capital letter.  An example funcntion name is  */
/*                CreateRedTree(a,b,c). Furthermore, each function name */
/*                should begin with a capital letter to easily distinguish */
/*                them from variables. */
/*                                                                     */
/*                Variable names: Each word in a variable name begins with */
/*                a capital letter EXCEPT the first letter of the variable */
/*                name.  For example, int newLongInt.  Global variables have */
/*                names beginning with "g".  An example of a global */
/*                variable name is gNewtonsConstant. */

/*  if DATA_TYPE is undefined then stack.h and stack.c will be code for */
/*  stacks of void *, if they are defined then they will be stacks of the */
/*  appropriate data_type */

#ifndef DATA_TYPE
#define DATA_TYPE void *
#endif

typedef struct stk_stack_node {
  DATA_TYPE info;
  struct stk_stack_node * next;
} stk_stack_node;

typedef struct stk_stack {
  stk_stack_node * top;
  stk_stack_node * tail;
} stk_stack ;

/*  These functions are all very straightforward and self-commenting so */
/*  I didn't think additional comments would be useful */
stk_stack * StackJoin(stk_stack * stack1, stk_stack * stack2);
stk_stack * StackCreate();
void StackPush(stk_stack * theStack, DATA_TYPE newInfoPointer);
void * StackPop(stk_stack * theStack);
int StackNotEmpty(stk_stack *);

