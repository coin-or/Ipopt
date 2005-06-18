// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpIpoptData.hpp,v 1.11 2005/05/04 20:54:12 andreasw Exp $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPIPOPTTYPE_HPP__
#define __IPIPOPTTYPE_HPP__

#include <list>
#include <string>
#include <iostream>
#include "IpRegOptions.hpp"

namespace Ipopt {

  class IpoptTypeInfo;
  std::list<IpoptTypeInfo*>& IpoptTypeInfosList();
  
/** This class is for implementing "class objects" or types in IPOPT.
 *  It is primarily used so we have runtime access to each class.
 */
class IpoptTypeInfo
{
 public:
  /** Constructors / Destructors */
  //@{
  /** Standard Constructor - takes in the type_name of the derived class */
  IpoptTypeInfo(const char* type_name)
    : type_name_(type_name)  
  {
    std::cout << "Registering: " << type_name_ << std::endl;
    IpoptTypeInfosList().push_back(this);
  }

  /** Standard destructor */
  virtual ~IpoptTypeInfo() {}
  //@}

  /** Standard Access Methods */
  //@{
  std::string TypeName() { return type_name_; }
  //@}

  /** Static methods to interact with ALL class type infos */
  //@{
  /** Register all the options */
  static void RegisterAllOptions(SmartPtr<RegisteredOptions> reg_options);
  //@}

 protected:
  /** Methods to be overridded by derived classes */
  //@{
  /** Override this method in derived TypeInfo classs to register options */
  virtual void RegisterOptionsImpl(SmartPtr<RegisteredOptions> reg_options)=0;
  //@}

 private:
  /** store the type name */
  std::string type_name_;

  /** keep a static list of all IpoptTypeInfo's */
  //  static std::list<IpoptTypeInfo*> ipopt_type_infos_;
 
};

inline 
 void IpoptTypeInfo::RegisterAllOptions(SmartPtr<RegisteredOptions> reg_options)
 {
   std::list<IpoptTypeInfo*>::iterator i;
   for (i=IpoptTypeInfosList().begin(); i != IpoptTypeInfosList().end(); i++) {
     (*i)->RegisterOptionsImpl(reg_options);
   }
 }


//extern std::list<IpoptTypeInfo*> IpoptTypeInfo::ipopt_type_infos_;


////////////////////////////////////////////////////////////
// Below is what needs to be added for each IpoptType
// but I think this can be macro'ed
///////////////////////////////////////////////////////////

#define DeclareIpoptType(__class_name__)                                \
class __class_name__ ## IpoptTypeInfo : public IpoptTypeInfo            \
 {									\
  public:                                                                \
   __class_name__ ## IpoptTypeInfo()                                     \
     : IpoptTypeInfo(# __class_name__) { }	\
    \
   __class_name__ ## IpoptTypeInfo* KeepCompilerFromOptimizing();	\
  protected:                                                             \
   virtual void RegisterOptionsImpl(SmartPtr<RegisteredOptions> reg_options); \
 };                                                                      

#define DefineIpoptType(__class_name__) \
  void __class_name__ ## IpoptTypeInfo::RegisterOptionsImpl(SmartPtr<RegisteredOptions> reg_options) \
   { \
     reg_options->SetRegisteringClass(# __class_name__ );		\
     __class_name__::RegisterOptions(reg_options);				\
     reg_options->SetRegisteringClass("Unknown Class");			\
   } \
 \
__class_name__ ## IpoptTypeInfo _ ## __class_name__ ## IpoptTypeInfo_; \
 \
__class_name__ ## IpoptTypeInfo* \
__class_name__ ## IpoptTypeInfo::KeepCompilerFromOptimizing() \
   { return &_ ## __class_name__ ## IpoptTypeInfo_; }

// class ClassAIpoptTypeInfo : public IpoptTypeInfo
// {
//  public:
//   ClassAIpoptTypeInfo()
//     : IpoptTypeInfo("ClassA") {}

//  protected:
//   virtual void RegisterOptionsImpl()
//   { 
//     //  ClassA::RegisterOptions(); 
//   }
// };
// ClassAIpoptTypeInfo _ClassAIpoptTypeInfo_;

} // end namespace Ipopt

#endif
