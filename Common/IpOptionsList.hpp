// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPOPTLIST_HPP__
#define __IPOPTLIST_HPP__

#include "IpUtils.hpp"
#include "IpReferenced.hpp"
#include "IpException.hpp"
#include "IpRegOptions.hpp"
#include <map>

#ifdef OLD_C_HEADERS
# include <stdio.h>
#else
# include <cstdio>
#endif

namespace Ipopt
{
  /** This class stores a list of user set options.  Each options is
   *  identified by a case-insensitive keyword (tag).  Its value is
   *  stored internally as a string (always lower case), but for
   *  convenience set and get methods are provided to obtain Index and
   *  Number type values.  For each keyword we also keep track of how
   *  often the value of an option has been requested by a get method.
   */
  class OptionsList : public ReferencedObject
  {
    /** Class for storing the value and counter for each option in
     *  OptionsList. */
    class OptionValue
    {
    public:
      /**@name Constructors/Destructors */
      //@{
      /** Default constructor (needed for the map) */
      OptionValue()
          :
          initialized_(false)
      {}

      /** Constructor given the value */
      OptionValue(std::string value)
          :
          value_(value),
          counter_(0),
          initialized_(true)
      {}

      /** Copy Constructor */
      OptionValue(const OptionValue& copy)
          :
          value_(copy.value_),
          counter_(copy.counter_),
          initialized_(copy.initialized_)
      {}

      /** Equals operator */
      void operator=(const OptionValue& copy)
      {
        value_=copy.value_;
        counter_=copy.counter_;
        initialized_=copy.initialized_;
      }

      /** Default Destructor */
      ~OptionValue()
      {}
      //@}

      /** Method for retrieving the value of an option.  Calling this
       *  method will increase the counter by one. */
      std::string GetValue() const
      {
        DBG_ASSERT(initialized_);
        counter_++;
        return value_;
      }

      /** Method for retrieving the value without increasing the
       *  counter */
      std::string Value() const
      {
        DBG_ASSERT(initialized_);
        return value_;
      }

      /** Method for accessing current value of the request counter */
      Index Counter() const
      {
        DBG_ASSERT(initialized_);
        return counter_;
      }

    private:
      /** Value for this option */
      std::string value_;

      /** Counter for requests */
      mutable Index counter_;

      /** for debugging */
      bool initialized_;
    };

  public:
    /**@name Constructors/Destructors */
    //@{
    OptionsList(SmartPtr<RegisteredOptions> reg_options, SmartPtr<Journalist> jnlst)
        : reg_options_(reg_options), jnlst_(jnlst)
    {}

    OptionsList()
    {}

    /** Copy Constructor */
    OptionsList(const OptionsList& copy)
    {
      options_ = copy.options_;
    }

    /** Default destructor */
    virtual ~OptionsList()
    {}

    /** Overloaded Equals Operator */
    void operator=(const OptionsList& source)
    {
      options_ = source.options_;
    }
    //@}

    /** @name Exceptions that can be used to indicate errors with
    options */
    //@{
    DECLARE_STD_EXCEPTION(OPTION_NOT_REGISTERED);
    DECLARE_STD_EXCEPTION(OPTION_VALUE_IS_INCORRECT_TYPE);
    DECLARE_STD_EXCEPTION(OPTION_OUT_OF_RANGE);
    DECLARE_STD_EXCEPTION(OPTION_VALUE_IS_NONINTEGER);
    DECLARE_STD_EXCEPTION(OPTION_VALUE_IS_NONNUMERIC);
    //@}

    /** @name Get / Set Methods */
    //@{
    void SetRegisteredOptions(const SmartPtr<RegisteredOptions> reg_options)
    {
      reg_options_ = reg_options;
    }
    void SetJournalist(const SmartPtr<Journalist> jnlst)
    {
      jnlst_ = jnlst;
    }
    //@}
    /** @name Methods for setting options */
    //@{
    void SetValue(const std::string& tag, const std::string& value);
    void SetNumericValue(const std::string& tag, Number value);
    void SetIntegerValue(const std::string& tag, Index value);
    //@}

    /** @name Method for retrieving values from the options list.  If
     *  a tag is not found, the methods return false. */
    //@{
    bool GetValue(const std::string& tag, std::string& value,
                  const std::string& prefix) const;
    bool GetEnumValue(const std::string& tag, Index& value,
                      const std::string& prefix) const;
    bool GetBoolValue(const std::string& tag, bool& value,
                      const std::string& prefix) const;
    bool GetNumericValue(const std::string& tag, Number& value,
                         const std::string& prefix) const;
    bool GetIntegerValue(const std::string& tag, Index& value,
                         const std::string& prefix) const;
    //@}

    /** Get a string with the list of all options (tag, value, counter) */
    void PrintList(std::string& list) const;

    /** Read options from a file with name filename.  Returns false if
     *  an error was encountered. */
    bool ReadFromFile(const Journalist& jnlst, FILE* fp);

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    //    OptionsList();

    //@}

    /** map for storing the options */
    std::map< std::string, OptionValue > options_;

    /** list of all the registered options to validate against */
    SmartPtr<RegisteredOptions> reg_options_;

    /** Journalist for writing error messages, etc. */
    SmartPtr<Journalist> jnlst_;

    /** auxilliary method for converting sting to all lower-case
     *  letters */
    const std::string& lowercase(const std::string tag) const;

    /** auxilliary method for finding the value for a tag in the
     *  options list.  This method first looks for the concatenated
     *  string prefix+tag (if prefix is not ""), and if this is not
     *  found, it looks for tag.  The return value is true iff
     *  prefix+tag or tag is found.  In that case, the corresponding
     *  string value is copied into value. */
    bool find_tag(const std::string& tag, const std::string& prefix,
                  std::string& value) const;

    /** read the next token from stream fp.  Returns false, if EOF was
     *  reached before a tokens was ecountered. */
    bool readnexttoken(FILE* fp, std::string& token);

    /** auxilliary string set by lowercase method */
    mutable std::string lowercase_buffer_;
  };

} // namespace Ipopt

#endif
