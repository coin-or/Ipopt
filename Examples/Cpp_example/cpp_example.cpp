// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpIpoptApplication.hpp"
#include "MyNLP.hpp"

using namespace Ipopt;

int main(int argv, char* argc[])
{
   SmartPtr<TNLP> mynlp = new MyNLP();
   IpoptApplication* app = new IpoptApplication();
   return (int) app->OptimizeTNLP(mynlp);
}