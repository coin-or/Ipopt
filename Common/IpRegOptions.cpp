// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOptionsList.hpp,v 1.2 2005/01/12 02:58:15 andreasw Exp $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpRegOptions.hpp"
#include <set>

namespace Ipopt {

  void RegisteredOption::OutputDescription(const Journalist& jnlst) const
  {
    std::string type_str = "Unknown";
    if (type_ ==OT_Number) { type_str = "Real Number"; }
    else if (type_ ==OT_Integer) { type_str = "Integer"; }
    else if (type_ ==OT_String) { type_str = "String"; }
    
    jnlst.Printf(J_ERROR, J_MAIN, "\n### %s (%s) ###\nClass: %s\nDescription: %s\n", 
		  name_.c_str(), type_str.c_str(), registering_class_.c_str(), description_.c_str()); 

    if (type_ ==OT_Number) { 
      if (has_lower_) { jnlst.Printf(J_ERROR, J_MAIN, "%g", lower_); } 
      else { jnlst.Printf(J_ERROR, J_MAIN, "-inf"); }

      if (lower_strict_) { jnlst.Printf(J_ERROR, J_MAIN, " < "); }
      else { jnlst.Printf(J_ERROR, J_MAIN, " <= "); }
      
      jnlst.Printf(J_ERROR, J_MAIN, "(%g)", default_number_);

      if (has_upper_ && upper_strict_) { jnlst.Printf(J_ERROR, J_MAIN, " < ");}
      else { jnlst.Printf(J_ERROR, J_MAIN, " <= "); }

      if (has_upper_) { jnlst.Printf(J_ERROR, J_MAIN, "%g\n", upper_); }
      else { jnlst.Printf(J_ERROR, J_MAIN, "+inf\n"); }
    }
    else if (type_ ==OT_Integer) {
      if (has_lower_) {	jnlst.Printf(J_ERROR, J_MAIN, "%d", (Index)lower_); }
      else { jnlst.Printf(J_ERROR, J_MAIN, "-inf"); }

      jnlst.Printf(J_ERROR, J_MAIN, " <= (%d) <= ", (Index)default_number_);

      if (has_upper_) {	jnlst.Printf(J_ERROR, J_MAIN, "%d\n", (Index)upper_); }
      else { jnlst.Printf(J_ERROR, J_MAIN, "+inf\n"); }
   }
    else if (type_ ==OT_String) { 
      std::vector<string_entry>::const_iterator i;
      jnlst.Printf(J_ERROR, J_MAIN, "Valid Settings:\n");      
      for (i = valid_strings_.begin(); i != valid_strings_.end(); i++) {
	jnlst.Printf(J_ERROR, J_MAIN, "\t%s (%s)\n", (*i).value_.c_str(), (*i).description_.c_str());      
      }
      jnlst.Printf(J_ERROR, J_MAIN, "Default: \"%s\"\n", default_string_.c_str());      
    }
  }

  void RegisteredOption::OutputShortDescription(const Journalist& jnlst) const
  {
    jnlst.Printf(J_ERROR, J_MAIN, "\n# %s\n", description_.c_str());
    jnlst.Printf(J_ERROR, J_MAIN, "%s: ",  name_.c_str());


    if (type_ ==OT_Number) { 
      if (has_lower_) { jnlst.Printf(J_ERROR, J_MAIN, "%g", lower_); } 
      else { jnlst.Printf(J_ERROR, J_MAIN, "-inf"); }

      if (lower_strict_) { jnlst.Printf(J_ERROR, J_MAIN, " < "); }
      else { jnlst.Printf(J_ERROR, J_MAIN, " <= "); }
      
      jnlst.Printf(J_ERROR, J_MAIN, "(%g)", default_number_);

      if (has_upper_ && upper_strict_) { jnlst.Printf(J_ERROR, J_MAIN, " < ");}
      else { jnlst.Printf(J_ERROR, J_MAIN, " <= "); }

      if (has_upper_) { jnlst.Printf(J_ERROR, J_MAIN, "%g\n", upper_); }
      else { jnlst.Printf(J_ERROR, J_MAIN, "+inf\n"); }
    }
    else if (type_ ==OT_Integer) {
      if (has_lower_) {	jnlst.Printf(J_ERROR, J_MAIN, "%d", (Index)lower_); }
      else { jnlst.Printf(J_ERROR, J_MAIN, "-inf"); }

      jnlst.Printf(J_ERROR, J_MAIN, " <= (%d) <= ", (Index)default_number_);

      if (has_upper_) {	jnlst.Printf(J_ERROR, J_MAIN, "%d\n", (Index)upper_); }
      else { jnlst.Printf(J_ERROR, J_MAIN, "+inf\n"); }
   }
    else if (type_ ==OT_String) { 
      std::vector<string_entry>::const_iterator i;
      jnlst.Printf(J_ERROR, J_MAIN, "(\"%s\")\n", default_string_.c_str());      
      for (i = valid_strings_.begin(); i != valid_strings_.end(); i++) {
	jnlst.Printf(J_ERROR, J_MAIN, "\t- %s (%s)\n", (*i).value_.c_str(), (*i).description_.c_str());      
      }
    }
  }

  void RegisteredOptions::AddNumberOption(const std::string& name, const std::string& description, 
					   Number default_value)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_Number);
    option->SetDefaultNumber(default_value);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddLowerBoundedNumberOption(const std::string& name, const std::string& description,  
						      Number lower, bool strict, Number default_value)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_Number);
    option->SetDefaultNumber(default_value);
    option->SetLowerNumber(lower, strict);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }
  void RegisteredOptions::AddUpperBoundedNumberOption(const std::string& name, const std::string& description,  
						      Number upper, bool strict, Number default_value)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_Number);
    option->SetDefaultNumber(default_value);
    option->SetUpperNumber(upper, strict);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddBoundedNumberOption(const std::string& name, const std::string& description,  
						 Number lower, bool lower_strict, 
						 Number upper, bool upper_strict, Number default_value)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_Number);
    option->SetDefaultNumber(default_value);
    option->SetLowerNumber(lower, lower_strict);
    option->SetUpperNumber(upper, upper_strict);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddIntegerOption(const std::string& name, const std::string& description,  Index default_value)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_Integer);
    option->SetDefaultInteger(default_value);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddLowerBoundedIntegerOption(const std::string& name, const std::string& description,  Index lower, Index default_value)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_Integer);
    option->SetDefaultInteger(default_value);
    option->SetLowerInteger(lower);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddUpperBoundedIntegerOption(const std::string& name, const std::string& description,  Index upper, Index default_value)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_Integer);
    option->SetDefaultInteger(default_value);
    option->SetUpperInteger(upper);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddBoundedIntegerOption(const std::string& name, const std::string& description,  Index lower, Index upper, Index default_value)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_Integer);
    option->SetDefaultInteger(default_value);
    option->SetLowerInteger(lower);
    option->SetUpperInteger(upper);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddStringOption(const std::string& name, const std::string& description,  const std::string& default_value, 
					  const std::vector<std::string>& settings, const std::vector<std::string>& descriptions)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    DBG_ASSERT(settings.size() == descriptions.size());
    for (int i=0; i<settings.size(); i++) {
      option->AddValidStringSetting(settings[i], descriptions[i]);
    }
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddStringOption1(const std::string& name, const std::string& description,  const std::string& default_value, 
					   const std::string& setting1, const std::string& description1)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

    void RegisteredOptions::AddStringOption2(const std::string& name, const std::string& description,  const std::string& default_value, 
			  const std::string& setting1, const std::string& description1,
			  const std::string& setting2, const std::string& description2)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    option->AddValidStringSetting(setting2, description2);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddStringOption3(const std::string& name, const std::string& description,  const std::string& default_value, 
			  const std::string& setting1, const std::string& description1,
			  const std::string& setting2, const std::string& description2,
			  const std::string& setting3, const std::string& description3)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    option->AddValidStringSetting(setting2, description2);
    option->AddValidStringSetting(setting3, description3);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddStringOption4(const std::string& name, const std::string& description,  const std::string& default_value, 
			  const std::string& setting1, const std::string& description1,
			  const std::string& setting2, const std::string& description2,
			  const std::string& setting3, const std::string& description3,
			  const std::string& setting4, const std::string& description4)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    option->AddValidStringSetting(setting2, description2);
    option->AddValidStringSetting(setting3, description3);
    option->AddValidStringSetting(setting4, description4);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddStringOption5(const std::string& name, const std::string& description,  const std::string& default_value, 
			  const std::string& setting1, const std::string& description1,
			  const std::string& setting2, const std::string& description2,
			  const std::string& setting3, const std::string& description3,
			  const std::string& setting4, const std::string& description4,
			  const std::string& setting5, const std::string& description5)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    option->AddValidStringSetting(setting2, description2);
    option->AddValidStringSetting(setting3, description3);
    option->AddValidStringSetting(setting4, description4);
    option->AddValidStringSetting(setting5, description5);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddStringOption6(const std::string& name, const std::string& description,  const std::string& default_value, 
			  const std::string& setting1, const std::string& description1,
			  const std::string& setting2, const std::string& description2,
			  const std::string& setting3, const std::string& description3,
			  const std::string& setting4, const std::string& description4,
			  const std::string& setting5, const std::string& description5,
			  const std::string& setting6, const std::string& description6)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    option->AddValidStringSetting(setting2, description2);
    option->AddValidStringSetting(setting3, description3);
    option->AddValidStringSetting(setting4, description4);
    option->AddValidStringSetting(setting5, description5);
    option->AddValidStringSetting(setting6, description6);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void RegisteredOptions::AddStringOption7(const std::string& name, const std::string& description,  const std::string& default_value, 
			  const std::string& setting1, const std::string& description1,
			  const std::string& setting2, const std::string& description2,
			  const std::string& setting3, const std::string& description3,
			  const std::string& setting4, const std::string& description4,
			  const std::string& setting5, const std::string& description5,
		          const std::string& setting7, const std::string& description7,
			  const std::string& setting6, const std::string& description6)
  {
    SmartPtr<RegisteredOption> option = new RegisteredOption(name, description, current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    option->AddValidStringSetting(setting2, description2);
    option->AddValidStringSetting(setting3, description3);
    option->AddValidStringSetting(setting4, description4);
    option->AddValidStringSetting(setting5, description5);
    option->AddValidStringSetting(setting6, description6);
    option->AddValidStringSetting(setting6, description7);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED, 
	       std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  SmartPtr<const RegisteredOption> RegisteredOptions::GetOption(const std::string& name)
  {
    SmartPtr<const RegisteredOption> option;
    std::map< std::string, SmartPtr<RegisteredOption> >::iterator reg_option = registered_options_.find(name);
    if (reg_option == registered_options_.end()) {
      option = NULL;
    }
    else {
      option = ConstPtr(reg_option->second);
    }

    return option;
  }
  
  void RegisteredOptions::OutputOptionDocumentation(const Journalist& jnlst)
  {
    // create a set to print sorted output
    std::set<std::string> classes;
    std::map< std::string, SmartPtr<RegisteredOption> >::iterator option;
    for (option = registered_options_.begin(); option != registered_options_.end(); option++) {
      classes.insert(option->second->RegisteringClass());
    }

    std::set<std::string>::iterator i;
    for (i = classes.begin(); i != classes.end(); i++) {
      jnlst.Printf(J_ERROR, J_MAIN, "\n\n### %s ###\n", (*i).c_str());
      for (option = registered_options_.begin(); option != registered_options_.end(); option++) {
	if (option->second->RegisteringClass() == (*i)) {
	  option->second->OutputShortDescription(jnlst);
	}
      }
    }

//     // This output should be sorted by registered class at some point, not by name
//     std::map< std::string, SmartPtr<RegisteredOption> >::iterator option;
//     for (option = registered_options_.begin(); option != registered_options_.end(); option++) {
//       option->second->OutputDescription(jnlst);
//     }
  }

} // namespace Ipopt
