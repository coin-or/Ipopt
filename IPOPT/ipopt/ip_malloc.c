/*
  Copyright (C) 2003, International Business Machines and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/* $Id: ip_malloc.c 531 2004-03-11 01:31:07Z andreasw $ */

#ifdef __cplusplus
extern "C" {
#endif

#include <config.h>

#include <stdio.h>
#if STDC_HEADERS
# include <stdlib.h>
# include <stddef.h>
#else
# if HAVE_STDLIB_H
#  include <stdlib.h>
# endif
#endif
#if HAVE_STRING_H
# if !STDC_HEADERS && HAVE_MEMORY_H
#  include <memory.h>
# endif
# include <string.h>
#endif
#if HAVE_STRINGS_H
# include <strings.h>
#endif
#if HAVE_INTTYPES_H
# include <inttypes.h>
#else
# if HAVE_STDINT_H
#  include <stdint.h>
# endif
#endif
#if HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <Ipopt.h>

fint F77_FUNC_(ip_malloc,IP_MALLOC)
           (fint *type, fint *length, fintptr *ptr)
{
  void *mem;
  size_t size;

  if( *ptr ) {
    mem = (void*)*ptr;
    free(mem);
    *ptr = 0;
  }
  if( *length != 0 ) {
    if( *type == IPMALLOC_INT )
      size=sizeof(fint);
    else if( *type == IPMALLOC_DOUBLE )
      size=sizeof(real);
    else
      return -4;
    mem = malloc((size_t)(*length)*(size));
    if( !mem )
      return -1;
    *ptr = (fintptr)mem;
  }
  return 0;
}

#ifdef __cplusplus
}
#endif
