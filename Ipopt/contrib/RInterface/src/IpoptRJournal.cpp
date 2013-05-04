/*
 * Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * file:   IpoptRJournal.cpp
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

#include "IpoptRJournal.hpp"

IpoptRJournal::IpoptRJournal( Ipopt::EJournalLevel default_level )
    : Journal("IpoptRJournal", default_level) { }

void IpoptRJournal::PrintImpl(Ipopt::EJournalCategory category, 
				 Ipopt::EJournalLevel level, const char* str) {
    
    // print string to R console
    Rprintf(str);
}

void IpoptRJournal::PrintfImpl(Ipopt::EJournalCategory category, 
				  Ipopt::EJournalLevel level, const char* pformat, 
				  va_list ap) { 
    
    // Define string
    const int MaxStrLen = 8192;
    char s[ MaxStrLen ];
    
    // R guarantees to have an implementation of vsnprintf available
    // http://www.mail-archive.com/r-devel@stat.math.ethz.ch/msg07054.html
    if ( vsnprintf( s, MaxStrLen, pformat, ap ) > MaxStrLen ) {
        Rprintf( "Warning: not all characters of next line are printed to the R console.\n" );
    }
    
    // print string to R console
    Rprintf( s );
}

void IpoptRJournal::FlushBufferImpl() {}
