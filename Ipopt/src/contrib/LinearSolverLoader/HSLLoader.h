/* Copyright (C) 2008 GAMS Development and others
 All Rights Reserved.
 This code is published under the Common Public License.

 $Id$

 Author: Stefan Vigerske
*/

#ifndef HSLLOADER_H_
#define HSLLOADER_H_

#ifdef __cplusplus
extern "C" {
#endif
  /** Tries to load a dynamically linked library with HSL routines.
   * Also tries to load symbols for those HSL routines that are not linked into Ipopt, i.e., HAVE_... is not defined. 
   * Return a failure if the library cannot be loaded, but not if a symbol is not found.
   * @see LSL_isMA27available
   * @see LSL_isMA28available
   * @see LSL_isMA57available
   * @see LSL_isMC19available
   * @param libname The name under which the HSL lib can be found, or NULL to use a default name (libhsl.<SHAREDLIBEXT>).
   * @param msgbuf A buffer where we can store a failure message. Assumed to be NOT NULL!
   * @param msglen Length of the message buffer.
   * @return Zero on success, nonzero on failure.
   */
  int LSL_loadHSL(const char* libname, char* msgbuf, int msglen);

  /** Unloads a loaded HSL library.
   * @return Zero on success, nonzero on failure.
   */
  int LSL_unloadHSL();

  /** Indicates whether a HSL library has been loaded.
   * @return Zero if not loaded, nonzero if handle is loaded
   */
  int LSL_isHSLLoaded();
  
  /** Indicates whether a HSL library is loaded and all symbols necessary to use MA27 have been found.
   * @return Zero if not available, nonzero if MA27 is available in the loaded library.
   */
  int LSL_isMA27available();

  /** Indicates whether a HSL library is loaded and all symbols necessary to use MA28 have been found.
   * @return Zero if not available, nonzero if MA28 is available in the loaded library.
   */
  int LSL_isMA28available();

  /** Indicates whether a HSL library is loaded and all symbols necessary to use MA57 have been found.
   * @return Zero if not available, nonzero if MA57 is available in the loaded library.
   */
  int LSL_isMA57available();
  
  /** Indicates whether a HSL library is loaded and all symbols necessary to use MA57 have been found.
   * @return Zero if not available, nonzero if MA57 is available in the loaded library.
   */
  int LSL_isMC19available();

  /** Returns name of the shared library that should contain HSL */
  char* LSL_HSLLibraryName();
  
#ifdef __cplusplus
}
#endif

#endif /*HSLLOADER_H_*/
