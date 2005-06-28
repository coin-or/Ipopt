// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpOptionsList.hpp"
#ifdef OLD_C_HEADERS
# include <ctype.h>
#else
# include <cctype>
#endif

namespace Ipopt
{

  void OptionsList::SetValue(const std::string& tag, const std::string& value)
  {
    if (IsValid(reg_options_)) { 
      SmartPtr<const RegisteredOption> option = reg_options_->GetOption(tag);
      
      if (IsNull(option)) {
	std::string msg = "Tried to set Option: " + tag;
	msg += ". It is not a valid option. Please check the list of available options.";
	THROW_EXCEPTION(OPTION_NOT_REGISTERED, msg);
      }

      if (option->Type() != OT_String) {
	std::string msg = "Tried to set Option: " + tag;
	msg += ". It is a valid option, but it is of type ";
	if (option->Type() == OT_Number) {
	  msg += " Number";
	}
	else if (option->Type() == OT_Integer) {
	  msg += " Integer";
	}
	else {
	  msg += " Unknown";
	}
	msg += ", not of type String. Please check the documentation for options.";
	if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	THROW_EXCEPTION(OPTION_VALUE_IS_INCORRECT_TYPE, msg);
      }

      if (!option->IsValidStringSetting(value)) {
	std::string msg = "Setting: " + value;
	msg += " is not a valid setting for Option: ";
	msg += tag;
	msg += ". Check the option documentation.";
	if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	THROW_EXCEPTION(OPTION_OUT_OF_RANGE, msg);
      }
    }

    OptionsList::OptionValue optval(lowercase(value));
    options_[lowercase(tag)] = optval;
  }

  void OptionsList::SetNumericValue(const std::string& tag, Number value)
  {
    char buffer[256];
    sprintf(buffer, "%g", value);

    if (IsValid(reg_options_)) { 
      SmartPtr<const RegisteredOption> option = reg_options_->GetOption(tag);
      
      if (IsNull(option)) {
	std::string msg = "Tried to set Option: " + tag;
	msg += ". It is not a valid option. Please check the list of available options.";
	THROW_EXCEPTION(OPTION_NOT_REGISTERED, msg);
      }

      if (option->Type() != OT_Number) {
	std::string msg = "Tried to set Option: " + tag;
	msg += ". It is a valid option, but it is of type ";
	if (option->Type() == OT_String) {
	  msg += " String";
	}
	else if (option->Type() == OT_Integer) {
	  msg += " Integer";
	}
	else {
	  msg += " Unknown";
	}
	msg += ", not of type Number. Please check the documentation for options.";
	if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	THROW_EXCEPTION(OPTION_VALUE_IS_INCORRECT_TYPE, msg);
      }

      if (!option->IsValidNumberSetting(value)) {
	std::string msg = "Setting: ";
	msg += buffer;
	msg += " is not a valid setting for Option: ";
	msg += tag;
	msg += ". Check the option documentation.";
	if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	THROW_EXCEPTION(OPTION_OUT_OF_RANGE, msg);
      }
    }

    OptionsList::OptionValue optval(buffer);
    options_[lowercase(tag)] = optval;
  }

  void OptionsList::SetIntegerValue(const std::string& tag, Index value)
  {
    char buffer[256];
    sprintf(buffer, "%d", value);

    if (IsValid(reg_options_)) { 
      SmartPtr<const RegisteredOption> option = reg_options_->GetOption(tag);
      
      if (IsNull(option)) {
	std::string msg = "Tried to set Option: " + tag;
	msg += ". It is not a valid option. Please check the list of available options.";
	THROW_EXCEPTION(OPTION_NOT_REGISTERED, msg);
      }

      if (option->Type() != OT_Integer) {
	std::string msg = "Tried to set Option: " + tag;
	msg += ". It is a valid option, but it is of type ";
	if (option->Type() == OT_String) {
	  msg += " String";
	}
	else if (option->Type() == OT_Number) {
	  msg += " Number";
	}
	else {
	  msg += " Unknown";
	}
	msg += ", not of type Integer. Please check the documentation for options.";
	if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	THROW_EXCEPTION(OPTION_VALUE_IS_INCORRECT_TYPE, msg);
      }

      if (!option->IsValidIntegerSetting(value)) {
	std::string msg = "Setting: ";
	msg += buffer;
	msg += " is not a valid setting for Option: ";
	msg += tag;
	msg += ". Check the option documentation.";
	if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	THROW_EXCEPTION(OPTION_OUT_OF_RANGE, msg);
      }
    }

    OptionsList::OptionValue optval(buffer);
    options_[lowercase(tag)] = optval;
  }

  bool OptionsList::GetValue(const std::string& tag, std::string& value,
			     const std::string& prefix) const
  {
    SmartPtr<const RegisteredOption> option = NULL;

    bool found = find_tag(tag, prefix, value);

    if (IsValid(reg_options_)) {
      option = reg_options_->GetOption(tag);
      if (IsNull(option)) {
	std::string msg = "IPOPT tried to get the value of Option: " + tag;
	msg += ". It is not a valid registered option.";
	  THROW_EXCEPTION(OPTION_NOT_REGISTERED, msg);
      }

      if (option->Type() != OT_String) {
	std::string msg = "IPOPT tried to get the value of Option: " + tag;
	msg += ". It is a valid option, but it is of type ";
	if (option->Type() == OT_Integer) {
	  msg += " Integer";
	}
	else if (option->Type() == OT_Number) {
	  msg += " Number";
	}
	else {
	  msg += " Unknown";
	}
	msg += ", not of type String. Please check the documentation for options.";
	if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	THROW_EXCEPTION(OPTION_VALUE_IS_INCORRECT_TYPE, msg);
      }

      if (found) {
	value = option->MapStringSetting(value);
      }
      else {
	value = option->DefaultString();
      }
    }

    return found;
  }

  bool OptionsList::GetEnumValue(const std::string& tag, Index& value,
				 const std::string& prefix) const
  {
    std::string str;
    SmartPtr<const RegisteredOption> option = NULL;

    bool found = find_tag(tag, prefix, str);

    if (IsValid(reg_options_)) {
      option = reg_options_->GetOption(tag);
      if (IsNull(option)) {
	std::string msg = "IPOPT tried to get the value of Option: " + tag;
	msg += ". It is not a valid registered option.";
	  THROW_EXCEPTION(OPTION_NOT_REGISTERED, msg);
      }

      if (option->Type() != OT_String) {
	std::string msg = "IPOPT tried to get the value of Option: " + tag;
	msg += ". It is a valid option, but it is of type ";
	if (option->Type() == OT_Integer) {
	  msg += " Integer";
	}
	else if (option->Type() == OT_Number) {
	  msg += " Number";
	}
	else {
	  msg += " Unknown";
	}
	msg += ", not of type String. Please check the documentation for options.";
	if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	THROW_EXCEPTION(OPTION_VALUE_IS_INCORRECT_TYPE, msg);
      }

      if (found) {
	value = option->MapStringSettingToEnum(str);
      }
      else {
	value = option->DefaultStringAsEnum();
      }
    }

    return found;
  }

  bool OptionsList::GetBoolValue(const std::string& tag, bool& value,
				 const std::string& prefix) const
  {
    std::string str;
    bool ret = GetValue(tag, str, prefix);
    if (str == "no" || str == "false" || str == "off") {
      value = false;
    }
    else if (str == "yes" || str == "true" || str == "on") {
      value = true;
    }
    else {
      THROW_EXCEPTION(OPTION_OUT_OF_RANGE, "Tried to get a boolean from an option and failed.");
      ret = false;
    }

    return ret;
  }

  bool OptionsList::GetNumericValue(const std::string& tag, Number& value,
				    const std::string& prefix) const
  { 
    SmartPtr<const RegisteredOption> option = NULL;

    if (IsValid(reg_options_)) {
      option = reg_options_->GetOption(tag);
      if (IsNull(option)) {
	std::string msg = "IPOPT tried to get the value of Option: " + tag;
	msg += ". It is not a valid registered option.";
	  THROW_EXCEPTION(OPTION_NOT_REGISTERED, msg);
      }

      if (option->Type() != OT_Number) {
	std::string msg = "IPOPT tried to get the value of Option: " + tag;
	msg += ". It is a valid option, but it is of type ";
	if (option->Type() == OT_Integer) {
	  msg += " Integer";
	}
	else if (option->Type() == OT_String) {
	  msg += " String";
	}
	else {
	  msg += " Unknown";
	}
	msg += ", not of type Number. Please check the documentation for options.";
	if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	THROW_EXCEPTION(OPTION_VALUE_IS_INCORRECT_TYPE, msg);
      }
    }

    std::string strvalue;
    if (find_tag(tag, prefix, strvalue)) {
      char* p_end;
      Number retval = strtod(strvalue.c_str(), &p_end);
      if (*p_end!='\0' && !isspace(*p_end)) {
	std::string msg = "Option \"" + tag +
	  "\": Double value expected, but non-numeric value \"" +
	  strvalue+"\" found.\n";
	THROW_EXCEPTION(OPTION_VALUE_IS_NONNUMERIC, msg);
      }
      value = retval;
      return true;
    }
    else if (IsValid(option)) {
      value = option->DefaultNumber();
      return false;
    }
    return false;
  }

  bool OptionsList::GetIntegerValue(const std::string& tag, Index& value,
				    const std::string& prefix) const
  {
    SmartPtr<const RegisteredOption> option = NULL;

    if (IsValid(reg_options_)) {
      option = reg_options_->GetOption(tag);
      if (IsNull(option)) {
	std::string msg = "IPOPT tried to get the value of Option: " + tag;
	msg += ". It is not a valid registered option.";
	  THROW_EXCEPTION(OPTION_NOT_REGISTERED, msg);
      }

      if (option->Type() != OT_Integer) {
	std::string msg = "IPOPT tried to get the value of Option: " + tag;
	msg += ". It is a valid option, but it is of type ";
	if (option->Type() == OT_Number) {
	  msg += " Number";
	}
	else if (option->Type() == OT_String) {
	  msg += " String";
	}
	else {
	  msg += " Unknown";
	}
	msg += ", not of type Integer. Please check the documentation for options.";
	if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	THROW_EXCEPTION(OPTION_VALUE_IS_INCORRECT_TYPE, msg);
      }
    }

    std::string strvalue;
    if (find_tag(tag, prefix, strvalue)) {
      char* p_end;
      Index retval = strtol(strvalue.c_str(), &p_end, 10);
      if (*p_end!='\0' && !isspace(*p_end)) {
	std::string msg = "Option \"" + tag +
	  "\": Integer value expected, but non-integer value \"" +
	  strvalue+"\" found.\n";
	THROW_EXCEPTION(OPTION_VALUE_IS_NONINTEGER, msg);
      }
      value = retval;
      return true;
    }
    else if (IsValid(option)) {
      value = option->DefaultInteger();
      return false;
    }

    return false;
  }

  const std::string& OptionsList::lowercase(const std::string tag) const
  {
    lowercase_buffer_ = tag;
    for(Index i=0; i<(Index)tag.length(); i++) {
      lowercase_buffer_[i] = tolower(tag[i]);
    }
    return lowercase_buffer_;
  }

  void OptionsList::PrintList(std::string& list) const
  {
    list.clear();
    char buffer[256];
    sprintf(buffer, "%40s   %-20s %s\n", "Name", "Value", "# times used");
    list += buffer;
    for(std::map< std::string, OptionValue >::const_iterator p = options_.begin();
	p != options_.end();
	p++ )
      {
	sprintf(buffer, "%40s = %-20s %6d\n", p->first.c_str(),
		p->second.Value().c_str(), p->second.Counter());
	list += buffer;
      }
  }

  bool OptionsList::ReadFromFile(const Journalist& jnlst,
				 FILE* fp)
  {
    DBG_ASSERT(fp);

    jnlst.Printf(J_DETAILED, J_MAIN, "Start reading options from file.\n");

    while (true) {
      std::string tag;
      std::string value;

      if (!readnexttoken(fp, tag)) {
	// That's it - end of file reached.
	jnlst.Printf(J_DETAILED, J_MAIN,
		     "Finished reading options from file.\n");
	return true;
      }

      if (!readnexttoken(fp, value)) {
	// Can't read value for a given tag
	jnlst.Printf(J_ERROR, J_MAIN,
		     "Error reading value for tag %s from file.\n",
		     tag.c_str());
	return false;
      }

      // Now add the value for the options list
      jnlst.Printf(J_DETAILED, J_MAIN,
		   "Adding option \"%s\" with value \"%s\" to OptionsList.\n",
		   tag.c_str(), value.c_str());

      if (IsValid(reg_options_)) {
	SmartPtr<const RegisteredOption> option = reg_options_->GetOption(tag);
	if (IsNull(option)) {
	  std::string msg = "Read Option: ";
	  msg += tag;
	  msg += ". It is not a valid option. Check the list of available options.";
	  THROW_EXCEPTION(OPTION_NOT_REGISTERED, msg);
	}
	
	if (option->Type() == OT_String) {
	  SetValue(tag, value);
	}
	else if (option->Type() == OT_Number) {
	  char* p_end;
	  Number retval = strtod(value.c_str(), &p_end);
	  if (*p_end!='\0' && !isspace(*p_end)) {
	    std::string msg = "Option \"" + tag +
	      "\": Double value expected, but non-numeric value \"" +
	      value + "\" found.\n";
	    THROW_EXCEPTION(OPTION_VALUE_IS_NONNUMERIC, msg);
	  }
	  SetNumericValue(tag, retval);
	}
	else if (option->Type() == OT_Integer) {
	  char* p_end;
	  Index retval = strtol(value.c_str(), &p_end, 10);
	  if (*p_end!='\0' && !isspace(*p_end)) {
	    std::string msg = "Option \"" + tag +
	      "\": Integer value expected, but non-integer value \"" +
	      value + "\" found.\n";
	    if (IsValid(jnlst_)) { option->OutputDescription(*jnlst_); }
	    THROW_EXCEPTION(OPTION_VALUE_IS_NONINTEGER, msg);
	  }
	  SetIntegerValue(tag, retval);
	}
	else {
	  DBG_ASSERT(false && "Option Type: Unknown");
	}
      }
      else {
	SetValue(tag, value);
      }
    }
  }

  bool OptionsList::find_tag(const std::string& tag,
			     const std::string& prefix,
			     std::string& value) const
  {
    bool found=false;
    std::map< std::string, OptionValue >::const_iterator p;

    if (prefix != "") {
      p = options_.find(lowercase(prefix+tag));
      if (p != options_.end()) {
	found = true;
      }
    }

    if (!found) {
      p = options_.find(lowercase(tag));
      if (p != options_.end()) {
	found = true;
      }
    }

    if (found) {
      value = p->second.GetValue();
    }

    return found;
  }

  bool OptionsList::readnexttoken(FILE* fp, std::string& token)
  {
    token.clear();
    int c = fgetc(fp);

    // First get rid of all comments and white spaces
    while (c!=EOF && (isspace(c) || c=='#') ) {
      if (c=='#') {
	for (c=fgetc(fp);
	     c!='\n' && c!=EOF;
	     c=fgetc(fp));
      }
      c=fgetc(fp);
    }

    // Now read the token
    while (c!=EOF && !isspace(c)) {
      token += c;
      c = fgetc(fp);
    }

  return (c!=EOF);
  }

} // namespace Ipopt

