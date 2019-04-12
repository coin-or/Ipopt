/*
  Copyright (C) 2002, Carnegie Mellon University and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/*
  $Id: main.c 531 2004-03-11 01:31:07Z andreasw $

  Simple Wrapper for Fortran main program fmain which passes the command
  line arguments to fmain dependend of the operating system

  Author:  Andreas Waechter    09-04-01
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __sgi
void fmain_(
         char *modnam,
         char *stubnam,
         int   modnamlen,
         int   stubnamlen
#elif __linux__
void fmain_(
         char *modnam,
         char *stubnam,
         int   modnamlen,
         int   stubnamlen
#elif __sun
void fmain_(
         char *modnam,
         char *stubnam,
         int   modnamlen,
         int   stubnamlen
#elif __osf__
void fmain_(
         char *modnam,
         char *stubnam,
         int   modnamlen,
         int   stubnamlen
#elif AIX
void fmain(
         char *modnam,
         char *stubnam,
         int   modnamlen,
         int   stubnamlen        
#elif WIN32
void __stdcall FMAIN(
         char *modnam,
         int   modnamlen,
         char *stubnam,
         int   stubnamlen        
#else
#error  Need to specify how characters are passed from C to Fortran
#endif
         );


main(int argc, char *argv[])
{
#define BUFLEN 256
  char MODELDAT[] = "MODEL.DAT";
  FILE *fp;

  char *modnam, *stubnam;
  int modnamlen, stubnamlen;

  if (argc <= 1)  /* read from MODEL.DAT */
    {
      fp = fopen(MODELDAT,"r");
      if( fp == NULL )
	{
	  fprintf(stderr, "Error opening %s.\n", MODELDAT);
	  exit(1);
	}
      modnam = malloc(BUFLEN);
      stubnam = malloc(BUFLEN);
      if( !fscanf(fp,"%s%s",modnam,stubnam) )
	{
	  fprintf(stderr, "Error reading %s.\n", MODELDAT);
	  exit(1);
	}
    }
  else if (argc == 2 )   /* modnam and stubnam are the same */
    {
      modnam = stubnam = argv[1];
    }
  else if (argc == 3)
    {
      modnam  = argv[1];
      stubnam = argv[2];
    }
  else
    {
      fprintf(stderr, "Too many arguments.\n");
      exit (1);
    }

  modnamlen  = strlen(modnam);
  stubnamlen = strlen(stubnam);

#ifdef __sgi
  fmain_( modnam, stubnam, modnamlen, stubnamlen );
#elif __linux__
  fmain_( modnam, stubnam, modnamlen, stubnamlen );
#elif __sun
  fmain_( modnam, stubnam, modnamlen, stubnamlen );
#elif __osf__
  fmain_( modnam, stubnam, modnamlen, stubnamlen );
#elif AIX
  fmain( modnam, stubnam, modnamlen, stubnamlen );
#elif WIN32
  FMAIN( modnam, modnamlen, stubnam, stubnamlen );
#else
#error  Need to specify how characters are passed from C to Fortran
#endif

  return 0;
}
