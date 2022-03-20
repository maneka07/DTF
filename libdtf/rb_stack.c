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

#include "rb_stack.h"

int StackNotEmpty(stk_stack * theStack) {
 // return( theStack ? (int) theStack->top : 0);
  return( theStack ? 1 : 0);
}

stk_stack * StackJoin(stk_stack * stack1, stk_stack * stack2) {
  if (!stack1->tail) {
    free(stack1);
    return(stack2);
  } else {
    stack1->tail->next=stack2->top;
    stack1->tail=stack2->tail;
    free(stack2);
    return(stack1);
  }
}

stk_stack * StackCreate() {
  stk_stack * newStack;

  newStack=(stk_stack *) SafeMalloc(sizeof(stk_stack));
  newStack->top=newStack->tail=NULL;
  return(newStack);
}


void StackPush(stk_stack * theStack, DATA_TYPE newInfoPointer) {
  stk_stack_node * newNode;

  if(!theStack->top) {
    newNode=(stk_stack_node *) SafeMalloc(sizeof(stk_stack_node));
    newNode->info=newInfoPointer;
    newNode->next=theStack->top;
    theStack->top=newNode;
    theStack->tail=newNode;
  } else {
    newNode=(stk_stack_node *) SafeMalloc(sizeof(stk_stack_node));
    newNode->info=newInfoPointer;
    newNode->next=theStack->top;
    theStack->top=newNode;
  }

}

DATA_TYPE StackPop(stk_stack * theStack) {
  DATA_TYPE popInfo;
  stk_stack_node * oldNode;

  if(theStack->top) {
    popInfo=theStack->top->info;
    oldNode=theStack->top;
    theStack->top=theStack->top->next;
    free(oldNode);
    if (!theStack->top) theStack->tail=NULL;
  } else {
    popInfo=NULL;
  }
  return(popInfo);
}

void StackDestroy(stk_stack * theStack,void DestFunc(void * a)) {
  stk_stack_node * x=theStack->top;
  stk_stack_node * y;

  if(theStack) {
    while(x) {
      y=x->next;
      DestFunc(x->info);
      free(x);
      x=y;
    }
    free(theStack);
  }
}

