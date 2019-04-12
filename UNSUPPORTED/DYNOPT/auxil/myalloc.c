/*
  Copyright (C) 2002, Carnegie Mellon University and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/* $Id: myalloc.c 531 2004-03-11 01:31:07Z andreasw $ */
#include <stdio.h>
#include <stdlib.h>

#ifdef AIX
#define myalloc_ myalloc
#endif

int myalloc_( int *size )
{
  void *mem;
  mem = malloc(*size);
  if( mem == NULL ) {
    fprintf(stderr,"Cannot allocate memory in myalloc.\n");
    exit(-1);
  }
  return (int)mem;
}
