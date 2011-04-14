// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-19


#include "SensUtils.hpp"
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

  SolverReturn AppReturn2SolverReturn(ApplicationReturnStatus ipopt_retval)
  {
    SolverReturn retval;
    switch (ipopt_retval) {
    case Solve_Succeeded:
      retval = SUCCESS;
      break;
    case Solved_To_Acceptable_Level:
      retval = STOP_AT_ACCEPTABLE_POINT;
      break;
    case Infeasible_Problem_Detected:
      retval = LOCAL_INFEASIBILITY;
      break;
    case Search_Direction_Becomes_Too_Small:
      retval = STOP_AT_TINY_STEP;
      break;
    case Diverging_Iterates:
      retval = DIVERGING_ITERATES;
      break;
    case User_Requested_Stop:
      retval = USER_REQUESTED_STOP;
      break;
    case Feasible_Point_Found:
      retval = FEASIBLE_POINT_FOUND;
      break;
    case Maximum_Iterations_Exceeded:
      retval = MAXITER_EXCEEDED;
      break;
    case Restoration_Failed:
      retval = RESTORATION_FAILURE;
      break;
    case Error_In_Step_Computation:
      retval = ERROR_IN_STEP_COMPUTATION;
      break;
    case Maximum_CpuTime_Exceeded:
      retval = CPUTIME_EXCEEDED;
      break;
    case Not_Enough_Degrees_Of_Freedom:
      retval = TOO_FEW_DEGREES_OF_FREEDOM;
      break;
    case Invalid_Problem_Definition:
      retval = UNASSIGNED;
      break;
    case Invalid_Option:
      retval = INVALID_OPTION;
      break;
    case Invalid_Number_Detected:
      retval = INVALID_NUMBER_DETECTED;
      break;
    case Unrecoverable_Exception:
      retval = UNASSIGNED;
      break;
    case NonIpopt_Exception_Thrown:
      retval = UNASSIGNED;
      break;
    case Insufficient_Memory:
      retval = OUT_OF_MEMORY;
      break;
    case Internal_Error:
      retval = INTERNAL_ERROR;
      break;
    default:
      retval = UNASSIGNED;
    }
    return retval;
  }
}
