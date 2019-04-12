/*
  Copyright (C) 2002, Carnegie Mellon University and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/*
  Perform changes from f2c output to file that can be used by ADOL-C
*/

/* $Id: f2c_post.c 531 2004-03-11 01:31:07Z andreasw $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFLEN 100

char* getline(char *line, int maxlen, FILE *stream)
{
  char buffer[BUFLEN+1];

  if( (fgets(line, maxlen, stream)) == NULL )
    return NULL;

  if( (strpbrk(line,"\n")) == 0 )
    {
      /*
	not yet line end reached - just ignore the rest
      */
      line[maxlen-1] = '\n';
      if( (fgets(buffer, BUFLEN, stream)) == NULL )
	return line;
      while( (strpbrk(buffer,"\n")) == 0 )
	if( (fgets(buffer, BUFLEN, stream)) == NULL )
	  return line;
    }

  return line;
}

void readerror(char *filename)
{
  fprintf(stderr, "Error reading '%s'. Abort\n", filename);
  exit(1);
}

void writeerror(char *filename)
{
  fprintf(stderr, "Error writinging '%s'. Abort\n", filename);
  exit(1);
}

main(int argc, char *argv[])
{
  FILE *infile, *outfile;
  char line[BUFLEN+1];
  const char F2CINCFILE[] = "f2c_adolc.h";

  const char KEYWORD1[] = "#ifdef __cplusplus";
  const char KEYWORD2[] = "    /* Builtin functions */";
  const char KEYWORD3[] = "static doublereal";
  const char KEYWORD4[] = "    extern /* Subroutine */ int condassign_(doublereal *,";
  const char EMPTY[] = "\n";

  char *inname, *outname;
  char *pos;

  if (argc != 3)
    {
      fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
      exit (1);
    }

  inname = argv[1];
  outname = argv[2];

  if ( (infile = fopen(argv[1],"r")) == NULL )
    {
      fprintf(stderr, "%s: Error opening '%s'. Abort\n",argv[0],inname);
      exit (1);
    }
  if ( (outfile = fopen(argv[2],"w")) == NULL )
    {
      fprintf(stderr, "%s: Error opening '%s'. Abort\n",argv[0],outname);
      exit (1);
    }

  while( strncmp((getline(line, BUFLEN, infile)),
		 KEYWORD1,strlen(KEYWORD1))!=0  );

  if( (getline(line, BUFLEN, infile)) == NULL ) readerror(inname);
  if( (getline(line, BUFLEN, infile)) == NULL ) readerror(inname);

  if( (fprintf(outfile, "#include \"%s\"\n",F2CINCFILE)) < 0 )
    writeerror(outname);
  if( (getline(line, BUFLEN, infile)) == NULL ) readerror(inname);

  while( (getline(line, BUFLEN, infile)) != NULL )
    {
      if( strncmp(line, KEYWORD1, strlen(KEYWORD1)) == 0  )
	/*
	  That's it - end of file reached
	*/
	{
	  printf("Successfully converted '%s' to '%s'.\n",inname,outname);
	  return 0;
	}

      if( strncmp(line, KEYWORD2, strlen(KEYWORD2)) == 0  )
	/*
	  Ignore the next non-empty lines
	*/
	{
	  if( (getline(line, BUFLEN, infile)) == NULL ) readerror(inname);
	  while( strncmp(line, EMPTY, strlen(EMPTY))!=0  )
	    if( (getline(line, BUFLEN, infile)) == NULL ) readerror(inname);
	}
      else if( strncmp(line, KEYWORD4, strlen(KEYWORD4)) == 0  )
	/*
	  Ignore the next line
	*/
	{
	  if( (getline(line, BUFLEN, infile)) == NULL ) readerror(inname);
	}
      else
	{
	  if( (pos = strstr(line,KEYWORD3)) != NULL )
	    /*
	      replace "static doublereal" by "static double"
	    */
	    {
	      *pos = '\0';
	      fprintf(outfile,"%s",line);
	      fprintf(outfile,"static double");
	      fprintf(outfile,"%s",(pos+strlen(KEYWORD3)));
	    }
	  else
	    fprintf(outfile,"%s",line);
	}
    }

  fprintf(stderr,"Error: Reached end of '%s' prematurely.  Abort.\n",inname);
  exit(1);
}
