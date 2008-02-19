/* Copyright (C) 2008   GAMS Development
 All Rights Reserved.
 This code is published under the Common Public License.

 $Id: HSLLoader.h 330 2008-02-09 16:50:07Z stefan $

 Authors:  Stefan Vigerske
*/

#ifndef HSLLOADER_H_
#define HSLLOADER_H_

/**
 * @return Zero on success, nonzero on failure.
 */
int LSL_loadHSL(const char* libname, char* msgbuf, int msglen);

/**
 * @return Zero on success, nonzero on failure.
 */
int LSL_unloadHSL();

/**
 * @return Zero if not loaded, nonzero if handle is loaded
 */
int LSL_isHSLLoaded();


#endif /*HSLLOADER_H_*/
