// Copyright (C) 2021 COIN-OR Foundation
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include "IpLibraryLoader.hpp"
#include "IpoptConfig.h"

#include <sstream>
#include <cstring>

#ifdef HAVE_WINDOWS_H
# include <windows.h>
#elif defined(HAVE_DLFCN_H)
# include <dlfcn.h>
#endif

namespace Ipopt
{

#ifdef IPOPT_HAS_LINEARSOLVERLOADER

#ifdef HAVE_WINDOWS_H
// add description of last error on windows to string
// see https://stackoverflow.com/questions/455434/how-should-i-use-formatmessage-properly-in-c
static
void addLastError(
   std::stringstream& s
)
{
   LPTSTR errorText = NULL;
   FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS | FORMAT_MESSAGE_ALLOCATE_BUFFER,
                 NULL, GetLastError(), 0, (LPTSTR)&errorText, 0, NULL);
   if( errorText != NULL )
   {
      s << errorText;
      LocalFree(errorText);
   }
   else
   {
      s << "see https://docs.microsoft.com/en-us/windows/win32/debug/system-error-codes)";
   }
}
#endif

void LibraryLoader::loadLibrary()
{
   if( libname.empty() )
   {
      THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, "No library name given (libname is empty)");
   }

#ifdef HAVE_WINDOWS_H
   libhandle = (void*)LoadLibrary(libname.c_str());
   if( libhandle == NULL )
   {
      std::stringstream s;
      s << "Error " << GetLastError() << " while loading DLL " << libname << ": ";
      addLastError(s);
      THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, s.str());
   }

#elif defined(HAVE_DLFCN_H)
   // ToDo switch to RTLD_LAZY for performance?
   libhandle = dlopen(libname.c_str(), RTLD_NOW);
   if( libhandle == NULL )
   {
      THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, dlerror());
   }

#else
   THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, "Do not know how to handle shared libraries on this operating system");
#endif
}

void LibraryLoader::unloadLibrary()
{
   if( libhandle == NULL )
   {
      return;
   }

#ifdef HAVE_WINDOWS_H
   if( FreeLibrary((HMODULE)libhandle) == 0 )
   {
      std::stringstream s;
      s << "Error " << GetLastError() << " while unloading " << libname << ": ";
      addLastError(s);
      THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, s.str());
   }

#elif defined(HAVE_DLFCN_H)
   if( dlclose(libhandle) != 0 )
   {
      THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, dlerror());
   }
#endif
}

/** tries to load symbol from a loaded library */
void* LibraryLoader::loadSymbol(
   const std::string& symbolname  /**< base name of symbol */
)
{
   if( libhandle == NULL )
   {
      loadLibrary();
   }
   DBG_ASSERT(libhandle != NULL);

   size_t len = symbolname.size();
   char* tripSym = new char[symbolname.size() + 2];
   void* symbol = NULL;

   for( int trip = 1; trip <= 6; trip++ )
   {
      switch( trip )
      {
         case 1: /* original */
            memcpy((void*)tripSym, (void*)symbolname.c_str(), len + 1);
            break;

         case 2: /* original_ */
            tripSym[len] = '_';
            tripSym[len + 1] = '\0';
            break;

         case 3: /* lower_ */
            for( size_t i = 0; i < len; ++i )
            {
               tripSym[i] = tolower(tripSym[i]);
            }
            break;

         case 4: /* lower */
            tripSym[len] = '\0';
            break;

         case 5: /* upper_ */
            for( size_t i = 0; i < len; ++i )
            {
               tripSym[i] = toupper(tripSym[i]);
            }
            tripSym[len] = '_';
            break;

         case 6: /* upper */
            tripSym[len] = '\0';
            break;

         default:
            ;
      }

#ifdef HAVE_WINDOWS_H
      symbol = (void*)GetProcAddress((HMODULE)libhandle, tripSym);
#elif defined(HAVE_DLFCN_H)
      symbol = dlsym(libhandle, tripSym);
#endif
      if( symbol != NULL )
      {
         break;
      }
   }

   delete[] tripSym;

   if( symbol == NULL )
   {
#ifdef HAVE_WINDOWS_H
      std::stringstream s;
      s << "Error " << GetLastError() << " while loading symbol " << symbolname << " from " << libname << ": ";
      addLastError(s);
      THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, s.str());
#elif defined(HAVE_DLFCN_H)
      THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, dlerror());
#endif
   }

   return symbol;
}

#else // IPOPT_HAS_LINEARSOLVERLOADER

void LibraryLoader::loadLibrary()
{
   THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, "Cannot load library at runtime. Ipopt has been build with --disable-linear-solver-loader.");
}

void LibraryLoader::unloadLibrary()
{
   if( libhandle == NULL )
   {
      return;
   }

   THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, "Cannot load library at runtime. Ipopt has been build with --disable-linear-solver-loader");
}

/** tries to load symbol from a loaded library */
void* LibraryLoader::loadSymbol(
   const std::string&
)
{
   THROW_EXCEPTION(DYNAMIC_LIBRARY_FAILURE, "Cannot load library at runtime. Ipopt has been build with --disable-linear-solver-loader");

   return NULL;
}

#endif // !IPOPT_HAS_LINEARSOLVERLOADER

} // namespace Ipopt
