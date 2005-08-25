/***********************************************************************
// Copyright (C) 2004, 2005, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
************************************************************************/

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* !!!!!!!!!!!!!!!! REMEMBER TO UPDATE IpReturnCodes.inc !!!!!!!!!!!!!!!! */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/** Return codes for the Optimize call for an application */
enum ApplicationReturnStatus
  {
    Solve_Succeeded=0,
    Solved_To_Acceptable_Level=1,
    Infeasible_Problem_Detected=2,
    Search_Direction_Becomes_Too_Small=3,

    Maximum_Iterations_Exceeded=-1,
    Restoration_Failed=-2,
    Not_Enough_Degrees_Of_Freedom=-10,
    Invalid_Problem_Definition=-11,
    Invalid_Option=-12,

    Unrecoverable_Exception=-100,
    NonIpopt_Exception_Thrown=-101,
    Insufficient_Memory=-102,
    Internal_Error=-199
  };
