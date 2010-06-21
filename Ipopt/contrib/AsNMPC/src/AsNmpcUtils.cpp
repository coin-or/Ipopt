// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Date   : 2009-05-19


#include "AsNmpcUtils.hpp"
#include <sstream>

namespace Ipopt 
{

  Index AsIndexMax(Index length, const Index* x, Index Incr)
  {
    if (length==0) {
      return 0;
    }
    Index maxval = x[0];
    for (Index i=1; i<length; i+=Incr) {
      if (x[i]>maxval) {
	maxval=x[i];
      }
    }
    return maxval;
  }

  Index AsIndexSum(Index length, const Index* x, Index Incr)
  {
    Index retval = 0;
    for (Index i=0; i<length; i+=Incr) {
      retval +=x[i];
    }
    return retval;
  }

  void append_Index(std::string& str, Index idx)
  {
    std::stringstream idx_stream;
    idx_stream << idx;
    std::string idx_string = idx_stream.str();
    str.append(idx_string);
  }

}

