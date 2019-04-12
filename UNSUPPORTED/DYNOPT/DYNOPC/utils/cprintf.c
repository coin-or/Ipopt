/*
  Copyright (C) 2002, Carnegie Mellon University and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/*
  $Id: cprintf.c 531 2004-03-11 01:31:07Z andreasw $

  Simply interface from FORTRAN to C's printf statement for character strings.
  This allows to perform output without newline characters at the end.

  Author:  Andreas Waechter     10-01-01

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef AIX
#define cprintf_ cprintf
#endif

#ifdef WIN32    // tell compiler that this function is to be exported
__declspec( dllexport ) void __stdcall CPRINTF(char* fformat,
					       int   lfformat,
					       char* fstring,
					       int   lfstring)
#else
void cprintf_(char* fformat,
	      char* fstring,
	      int   lfformat,
	      int   lfstring)
#endif
{
  char *string, *format;

  if( (format = (char *)malloc(lfformat+1)) &&
      (string = (char *)malloc(lfstring+1)) )
    {
      memcpy(format,fformat,lfformat);
      memcpy(string,fstring,lfstring);
      format[lfformat] = '\0';
      string[lfstring] = '\0';

      setbuf(stdout, NULL);
      fprintf(stdout,format,string);

      free(format);
      free(string);
    }
  else
    {
      fprintf(stderr,"Error in cprintf.c allocating memory!  Abort.\n");
      exit(2);
    }
}
