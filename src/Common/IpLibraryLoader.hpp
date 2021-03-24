// Copyright (C) 2021 COIN-OR Foundation
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef __IPLIBRARYLOADER_HPP__
#define __IPLIBRARYLOADER_HPP__

#include "IpReferenced.hpp"
#include "IpException.hpp"

namespace Ipopt
{

/** loading of a library at runtime
 *
 * wrapper around dlopen()/dlsym() and variants
 */
class IPOPTLIB_EXPORT LibraryLoader : public ReferencedObject
{
private:
   void* libhandle;

   /** unimplemented copy constructor */
   LibraryLoader(const LibraryLoader&);
   /** unimplemented assigment operator */
   LibraryLoader& operator=(const LibraryLoader&);

public:
   /** constructor */
   LibraryLoader()
   : libhandle(NULL)
   { }

   /** destructor */
   ~LibraryLoader()
   {
      unloadLibrary();
   }

   /** tries to load library */
   void loadLibrary(
      const std::string& libname  /**< full name of library, can include path */
      );

   /** unload library, if loaded */
   void unloadLibrary();

   /** tries to load symbol from a loaded library */
   void* loadSymbol(
      const std::string& symbolname  /**< base name of symbol */
      );
};

/** a problem occurred with a a dynamically loaded library
 *
 * e.g., library could not be loaded or a symbol could not be found
 */
DECLARE_STD_EXCEPTION(DYNAMIC_LIBRARY_FAILURE);

} // namespace Ipopt

#endif /* __IPLIBRARYLOADER_HPP__ */
