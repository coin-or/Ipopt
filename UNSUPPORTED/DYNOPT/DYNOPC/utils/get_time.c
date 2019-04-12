/*
  Copyright (C) 2002, Carnegie Mellon University and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/* $Id: get_time.c 531 2004-03-11 01:31:07Z andreasw $
   Read current time into Fortran character string
   Author:    Andreas Waechter    09-27-01
*/
#include <time.h>
#include <string.h>

#ifdef WIN32
__declspec( dllexport ) void __stdcall GET_TIME
#elif AIX
void get_time
#else
void get_time_
#endif
  ( char* datestr, 
    int  ldatestr )
  
{
  char* retstr;
  time_t tp;
  int i;

  tp = time(NULL);
  retstr = ctime(&tp);
  strncpy(datestr, retstr, ldatestr);
  for( i=strlen(retstr); i < ldatestr; datestr[i++] = ' ' );
}

