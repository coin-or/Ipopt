// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOptionsList.hpp,v 1.2 2005/01/12 02:58:15 andreasw Exp $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpRegOptions.hpp"
#include <set>

#ifdef OLD_C_HEADERS
# include <ctype.h>
#else
# include <cctype>
#endif

namespace Ipopt
{

  Index RegisteredOption::next_counter_ = 1;

  void RegisteredOption::OutputDescription(const Journalist& jnlst) const
  {
    std::string type_str = "Unknown";
    if (type_ ==OT_Number) {
      type_str = "Real Number";
    }
    else if (type_ ==OT_Integer) {
      type_str = "Integer";
    }
    else if (type_ ==OT_String) {
      type_str = "String";
    }

    jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                 "\n### %s (%s) ###\nClass: %s\nDescription: %s\n",
                 name_.c_str(), type_str.c_str(),
                 registering_class_.c_str(), short_description_.c_str());

    if (type_ ==OT_Number) {
      if (has_lower_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%g", lower_);
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "-inf");
      }

      if (lower_strict_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " < ");
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " <= ");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "(%g)", default_number_);

      if (has_upper_ && upper_strict_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " < ");
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " <= ");
      }

      if (has_upper_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%g\n", upper_);
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "+inf\n");
      }
    }
    else if (type_ ==OT_Integer) {
      if (has_lower_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%d", (Index)lower_);
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "-inf");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " <= (%d) <= ", (Index)default_number_);

      if (has_upper_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%d\n", (Index)upper_);
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "+inf\n");
      }
    }
    else if (type_ ==OT_String) {
      std::vector<string_entry>::const_iterator i;
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "Valid Settings:\n");
      for (i = valid_strings_.begin(); i != valid_strings_.end(); i++) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\t%s (%s)\n",
                     (*i).value_.c_str(), (*i).description_.c_str());
      }
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "Default: \"%s\"\n",
                   default_string_.c_str());
    }
  }

  void RegisteredOption::OutputShortDescription(const Journalist& jnlst) const
  {
    jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%-30s",  name_.c_str());


    if (type_ == OT_Number) {
      if (has_lower_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%10g", lower_);
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%10s", "-inf");
      }

      if (has_lower_ && !lower_strict_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " <= ");
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " <  ");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "(%11g)", default_number_);

      if (has_upper_ && !upper_strict_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " <= ");
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " <  ");
      }

      if (has_upper_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%-10g\n", upper_);
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%-10s\n", "+inf");
      }
    }
    else if (type_ == OT_Integer) {
      if (has_lower_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%10d <= ", (Index)lower_);
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%10s <  ", "-inf");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "(%11d)",
                   (Index)default_number_);

      if (has_upper_) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " <= %-10d\n", (Index)upper_);
      }
      else {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " <  %-10s\n", "+inf");
      }
    }
    else if (type_ == OT_String) {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "(\"%s\")\n",
                   default_string_.c_str());
    }
    jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "   ");
    jnlst.PrintStringOverLines(J_SUMMARY, J_DOCUMENTATION, 3, 76,
                               short_description_.c_str());
    if (long_description_ != "") {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n     ");
      jnlst.PrintStringOverLines(J_SUMMARY, J_DOCUMENTATION, 5, 74,
                                 long_description_.c_str());
    }
    if (type_ == OT_String) {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n   Possible values:\n");
      for (std::vector<string_entry>::const_iterator
           i = valid_strings_.begin();
           i != valid_strings_.end(); i++) {
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "    - %-23s [",
                     (*i).value_.c_str());

        jnlst.PrintStringOverLines(J_SUMMARY, J_DOCUMENTATION, 31, 48,
                                   (*i).description_.c_str());
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "]\n");
      }
    }
    else {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n");
    }
    jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n");
  }

  bool RegisteredOption::IsValidStringSetting(const std::string& value) const
  {
    DBG_ASSERT(type_ == OT_String);

    std::vector<string_entry>::const_iterator i;
    for (i = valid_strings_.begin(); i != valid_strings_.end(); i++) {
      if (i->value_ == "*" || string_equal_insensitive(i->value_, value)) {
        return true;
      }
    }
    return false;
  }

  std::string
  RegisteredOption::MapStringSetting(const std::string& value) const
  {
    DBG_ASSERT(type_ == OT_String);

    std::string matched_setting = "";

    std::vector<string_entry>::const_iterator i;
    for (i = valid_strings_.begin(); i != valid_strings_.end(); i++) {
      if (i->value_ == "*") {
        matched_setting = value;
      }
      else if (string_equal_insensitive(i->value_, value)) {
        matched_setting = i->value_;
      }
    }
    return matched_setting;
  }

  Index
  RegisteredOption::MapStringSettingToEnum(const std::string& value) const
  {
    DBG_ASSERT(type_ == OT_String);

    Index matched_setting = -1;

    Index cnt = 0;
    std::vector<string_entry>::const_iterator i;
    for (i = valid_strings_.begin(); i != valid_strings_.end(); i++) {
      ASSERT_EXCEPTION(i->value_ != "*", IpoptException,
                       "Cannot map a wildcard setting to an enumeration");
      if (string_equal_insensitive(i->value_, value)) {
        matched_setting = cnt;
        break;
      }
      cnt++;
    }

    ASSERT_EXCEPTION(matched_setting != -1, ERROR_CONVERTING_STRING_TO_ENUM,
                     std::string("Could not find a match for setting ") + value +
                     " in option: " + name_);
    return matched_setting;
  }

  bool
  RegisteredOption::string_equal_insensitive(const std::string& s1,
      const std::string& s2) const
  {
    using namespace std;

    if (s1.size()!=s2.size())
      return false;

    string::const_iterator i1 = s1.begin();
    string::const_iterator i2 = s2.begin();

    while(i1!=s1.end()) {
      if (toupper(*i1)!=toupper(*i2))
        return false;
      i1++;
      i2++;
    }
    return true;
  }

  void
  RegisteredOptions::AddNumberOption(const std::string& name,
                                     const std::string& short_description,
                                     Number default_value,
                                     const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_Number);
    option->SetDefaultNumber(default_value);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(),
                     OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddLowerBoundedNumberOption(const std::string& name,
      const std::string& short_description,
      Number lower, bool strict,
      Number default_value,
      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_Number);
    option->SetDefaultNumber(default_value);
    option->SetLowerNumber(lower, strict);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddUpperBoundedNumberOption(const std::string& name,
      const std::string& short_description,
      Number upper, bool strict,
      Number default_value,
      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_Number);
    option->SetDefaultNumber(default_value);
    option->SetUpperNumber(upper, strict);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddBoundedNumberOption(const std::string& name,
      const std::string& short_description,
      Number lower, bool lower_strict,
      Number upper, bool upper_strict,
      Number default_value,
      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_Number);
    option->SetDefaultNumber(default_value);
    option->SetLowerNumber(lower, lower_strict);
    option->SetUpperNumber(upper, upper_strict);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddIntegerOption(const std::string& name,
                                      const std::string& short_description,
                                      Index default_value,
                                      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_Integer);
    option->SetDefaultInteger(default_value);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddLowerBoundedIntegerOption(const std::string& name,
      const std::string& short_description,
      Index lower, Index default_value,
      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_Integer);
    option->SetDefaultInteger(default_value);
    option->SetLowerInteger(lower);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddUpperBoundedIntegerOption(const std::string& name,
      const std::string& short_description,
      Index upper, Index default_value,
      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_Integer);
    option->SetDefaultInteger(default_value);
    option->SetUpperInteger(upper);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddBoundedIntegerOption(const std::string& name,
      const std::string& short_description,
      Index lower, Index upper,
      Index default_value,
      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_Integer);
    option->SetDefaultInteger(default_value);
    option->SetLowerInteger(lower);
    option->SetUpperInteger(upper);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddStringOption(const std::string& name,
                                     const std::string& short_description,
                                     const std::string& default_value,
                                     const std::vector<std::string>& settings,
                                     const std::vector<std::string>& descriptions,
                                     const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    DBG_ASSERT(settings.size() == descriptions.size());
    for (int i=0; i<(int)settings.size(); i++) {
      option->AddValidStringSetting(settings[i], descriptions[i]);
    }
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddStringOption1(const std::string& name,
                                      const std::string& short_description,
                                      const std::string& default_value,
                                      const std::string& setting1,
                                      const std::string& description1,
                                      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddStringOption2(const std::string& name,
                                      const std::string& short_description,
                                      const std::string& default_value,
                                      const std::string& setting1,
                                      const std::string& description1,
                                      const std::string& setting2,
                                      const std::string& description2,
                                      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    option->AddValidStringSetting(setting2, description2);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddStringOption3(const std::string& name,
                                      const std::string& short_description,
                                      const std::string& default_value,
                                      const std::string& setting1,
                                      const std::string& description1,
                                      const std::string& setting2,
                                      const std::string& description2,
                                      const std::string& setting3,
                                      const std::string& description3,
                                      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    option->AddValidStringSetting(setting2, description2);
    option->AddValidStringSetting(setting3, description3);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  void
  RegisteredOptions::AddStringOption4(const std::string& name,
                                      const std::string& short_description,
                                      const std::string& default_value,
                                      const std::string& setting1,
                                      const std::string& description1,
                                      const std::string& setting2,
                                      const std::string& description2,
                                      const std::string& setting3,
                                      const std::string& description3,
                                      const std::string& setting4,
                                      const std::string& description4,
                                      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
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

  void
  RegisteredOptions::AddStringOption5(const std::string& name,
                                      const std::string& short_description,
                                      const std::string& default_value,
                                      const std::string& setting1,
                                      const std::string& description1,
                                      const std::string& setting2,
                                      const std::string& description2,
                                      const std::string& setting3,
                                      const std::string& description3,
                                      const std::string& setting4,
                                      const std::string& description4,
                                      const std::string& setting5,
                                      const std::string& description5,
                                      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
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

  void
  RegisteredOptions::AddStringOption6(const std::string& name,
                                      const std::string& short_description,
                                      const std::string& default_value,
                                      const std::string& setting1,
                                      const std::string& description1,
                                      const std::string& setting2,
                                      const std::string& description2,
                                      const std::string& setting3,
                                      const std::string& description3,
                                      const std::string& setting4,
                                      const std::string& description4,
                                      const std::string& setting5,
                                      const std::string& description5,
                                      const std::string& setting6,
                                      const std::string& description6,
                                      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
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

  void
  RegisteredOptions::AddStringOption7(const std::string& name,
                                      const std::string& short_description,
                                      const std::string& default_value,
                                      const std::string& setting1,
                                      const std::string& description1,
                                      const std::string& setting2,
                                      const std::string& description2,
                                      const std::string& setting3,
                                      const std::string& description3,
                                      const std::string& setting4,
                                      const std::string& description4,
                                      const std::string& setting5,
                                      const std::string& description5,
                                      const std::string& setting6,
                                      const std::string& description6,
                                      const std::string& setting7,
                                      const std::string& description7,
                                      const std::string& long_description)
  {
    SmartPtr<RegisteredOption> option =
      new RegisteredOption(name, short_description, long_description,
                           current_registering_class_);
    option->SetType(OT_String);
    option->SetDefaultString(default_value);
    option->AddValidStringSetting(setting1, description1);
    option->AddValidStringSetting(setting2, description2);
    option->AddValidStringSetting(setting3, description3);
    option->AddValidStringSetting(setting4, description4);
    option->AddValidStringSetting(setting5, description5);
    option->AddValidStringSetting(setting6, description6);
    option->AddValidStringSetting(setting7, description7);
    ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                     std::string("The option: ") + option->Name() + " has already been registered by someone else");
    registered_options_[name] = option;
  }

  SmartPtr<const RegisteredOption> RegisteredOptions::GetOption(const std::string& name)
  {
    std::string tag_only = name;
    std::string::size_type pos = name.rfind(".", name.length());
    if (pos != std::string::npos) {
      tag_only = name.substr(pos+1, name.length()-pos);
    }
    SmartPtr<const RegisteredOption> option;
    std::map< std::string, SmartPtr<RegisteredOption> >::iterator reg_option = registered_options_.find(tag_only);
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
    std::set
      <std::string> classes;
    std::map <std::string, SmartPtr<RegisteredOption> >::iterator option;
    for (option = registered_options_.begin(); option != registered_options_.end(); option++) {
      classes.insert(option->second->RegisteringClass());
    }

    std::set
      <std::string>::iterator i;
    for (i = classes.begin(); i != classes.end(); i++) {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "\n### %s ###\n\n", (*i).c_str());
      std::map<Index, SmartPtr<RegisteredOption> > class_options;
      for (option = registered_options_.begin();
           option != registered_options_.end(); option++) {
        if (option->second->RegisteringClass() == (*i)) {

          class_options[option->second->Counter()] = option->second;
        }
      }
      std::map<Index, SmartPtr<RegisteredOption> >::const_iterator co;
      for (co = class_options.begin(); co != class_options.end(); co++) {
        co->second->OutputShortDescription(jnlst);
      }
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n");
    }
  }

} // namespace Ipopt
