/* Copyright (C) 2008   GAMS Development
 All Rights Reserved.
 This code is published under the Common Public License.

 $Id: PardisoLoader.h 330 2008-02-09 16:50:07Z stefan $

 Authors:  Stefan Vigerske
*/

#ifndef PARDISOLOADER_H_
#define PARDISOLOADER_H_

/**
 * @return Zero on success, nonzero on failure.
 */
int LSL_loadPardisoLib(const char* libname, char* msgbuf, int msglen);

/**
 * @return Zero on success, nonzero on failure.
 */
int LSL_unloadPardisoLib();

/**
 * @return Zero if not loaded, nonzero if handle is loaded
 */
int LSL_isPardisoLoaded();

#endif /*PARADISOLOADER_H_*/
