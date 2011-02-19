/*
 * Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * file:   IpoptRJournal.hpp
 * author: Jelmer Ypma
 * date:   30 January 2011
 *
 * This file defines a C++ class that takes care of re-directing
 * output to the R terminal. Needed for Windows.
 *
 * Financial support of the UK Economic and Social Research Council 
 * through a grant (RES-589-28-0001) to the ESRC Centre for Microdata 
 * Methods and Practice (CeMMAP) is gratefully acknowledged.
 */

#ifndef __IpoptRJournal_HPP__
#define __IpoptRJournal_HPP__

#include "IpJournalist.hpp"     // ISA  Journal
#include <R.h>                  // USES Rprintf

class IpoptRJournal : public Ipopt::Journal {
    public:

    // The constructor.
    IpoptRJournal( Ipopt::EJournalLevel default_level );

    // The destructor.
    virtual ~IpoptRJournal() { };

    protected:

    // These functions override the functions in the Journal class.
    virtual void PrintImpl( 
                Ipopt::EJournalCategory category, 
                Ipopt::EJournalLevel level, 
			    const char* str);
                
    virtual void PrintfImpl(
                Ipopt::EJournalCategory category, 
                Ipopt::EJournalLevel level, 
			    const char* pformat, 
                va_list ap);
                
    virtual void FlushBufferImpl();
};

#endif
