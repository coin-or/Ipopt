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
    OptionsList::OptionValue optval(lowercase(value));
    options_[lowercase(tag)] = optval;
  }

  void OptionsList::SetNumericValue(const std::string& tag, Number value)
  {
    char buffer[256];
    sprintf(buffer, "%g", value);
    OptionsList::OptionValue optval(buffer);
    options_[lowercase(tag)] = optval;
  }

  void OptionsList::SetIntegerValue(const std::string& tag, Index value)
  {
    char buffer[256];
    sprintf(buffer, "%d", value);
    OptionsList::OptionValue optval(buffer);
    options_[lowercase(tag)] = optval;
  }

  bool OptionsList::GetValue(const std::string& tag, std::string& value,
			     const std::string& prefix) const
  {
    return find_tag(tag, prefix, value);
  }

  bool OptionsList::GetNumericValue(const std::string& tag, Number& value,
				    const std::string& prefix) const
  {
    std::string strvalue;
    if (!find_tag(tag, prefix, strvalue)) {
      return false;
    }

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

  bool OptionsList::GetIntegerValue(const std::string& tag, Index& value,
				    const std::string& prefix) const
  {
    std::string strvalue;
    if (!find_tag(tag, prefix, strvalue)) {
      return false;
    }

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
      SetValue(tag, value);
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
    char c = fgetc(fp);

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

