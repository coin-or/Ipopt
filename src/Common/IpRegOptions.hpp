// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-06-18

#ifndef __IPREGOPTIONS_HPP__
#define __IPREGOPTIONS_HPP__

#include "IpUtils.hpp"
#include "IpReferenced.hpp"
#include "IpException.hpp"
#include "IpSmartPtr.hpp"

#include <map>
#include <set>
#include <list>

namespace Ipopt
{

enum RegisteredOptionType
{
   OT_Number,
   OT_Integer,
   OT_String,
   OT_Unknown
};

class RegisteredOption;

/** A category of registered options. */
class IPOPTLIB_EXPORT RegisteredCategory: public ReferencedObject
{
   friend class RegisteredOptions;
public:
   /// Constructor
   ///
   /// Use negative value for priority to suppress it being included in documentation.
   RegisteredCategory(
      std::string name,
      int         priority
      )
   : name_(name),
     priority_(priority)
   { }

   /// name of category
   const std::string& Name() const
   {
      return name_;
   }

   /// priority of category
   int Priority() const
   {
      return priority_;
   }

   /// gives list of options in this category
   const std::list<SmartPtr<RegisteredOption> >& RegisteredOptions() const
   {
      return regoptions_;
   }

   // class comparing two categories by priority
   class ComparePriority
   {
   public:
      bool operator()(
         const SmartPtr<RegisteredCategory>& lhs,
         const SmartPtr<RegisteredCategory>& rhs
         ) const
      {
         DBG_ASSERT(IsValid(lhs));
         DBG_ASSERT(IsValid(rhs));
         return lhs->priority_ < rhs->priority_;
      }
   };

private:
   /// unimplemented default constructor
   RegisteredCategory();
   /// unimplemented copy constructor
   RegisteredCategory(const RegisteredCategory&);
   /// unimplemented assignment operator
   RegisteredCategory& operator=(const RegisteredCategory&);

   /// name of category
   std::string name_;

   /// priority of category (used to decide whether to print and printing order)
   int priority_;

   /// options of this category
   std::list<SmartPtr<RegisteredOption> > regoptions_;
};

/** Option that has been registered. */
class IPOPTLIB_EXPORT RegisteredOption: public ReferencedObject
{
   friend class RegisteredOptions;
public:
   /** class to hold the valid string settings for a string option */
   class string_entry
   {
   public:
      string_entry(
         const std::string& value,
         const std::string& description
      )
         : value_(value),
           description_(description)
      { }

      std::string value_;
      std::string description_;
   };

   /** Constructors / Destructors */
   ///@{
   RegisteredOption(
      Index counter
   )
      : type_(OT_Unknown),
        advanced_(false),
        has_lower_(false),
        has_upper_(false),
        counter_(counter)
   { }

   RegisteredOption(
      const std::string& name,
      const std::string& short_description,
      const std::string& long_description,
      const SmartPtr<RegisteredCategory>& registering_category,
      Index counter,
      bool advanced = false
   )
      : name_(name),
        short_description_(short_description),
        long_description_(long_description),
        registering_category_(registering_category),
        type_(OT_Unknown),
        advanced_(advanced),
        has_lower_(false),
        has_upper_(false),
        counter_(counter)
   { }

   RegisteredOption(
      const RegisteredOption& copy
   )
      : name_(copy.name_),
        short_description_(copy.short_description_),
        long_description_(copy.long_description_),
        registering_category_(copy.registering_category_),
        type_(copy.type_),
        advanced_(copy.advanced_),
        has_lower_(copy.has_lower_),
        lower_(copy.lower_),
        has_upper_(copy.has_upper_),
        upper_(copy.upper_),
        valid_strings_(copy.valid_strings_),
        counter_(copy.counter_)
   { }

   virtual ~RegisteredOption()
   { }
   ///@}

   DECLARE_STD_EXCEPTION(ERROR_CONVERTING_STRING_TO_ENUM);

   /** Standard Get / Set Methods */
   ///@{
   /** Get the option's name (tag in the input file) */
   virtual const std::string& Name() const
   {
      return name_;
   }
   /** Set the option's name (tag in the input file) */
   virtual void SetName(
      const std::string& name
   )
   {
      name_ = name;
   }
   /** Get the short description */
   virtual const std::string& ShortDescription() const
   {
      return short_description_;
   }

   /** Get the long description */
   virtual const std::string& LongDescription() const
   {
      return long_description_;
   }

   /** Set the short description */
   virtual void SetShortDescription(
      const std::string& short_description
   )
   {
      short_description_ = short_description;
   }

   /** Set the long description */
   virtual void SetLongDescription(
      const std::string& long_description
   )
   {
      long_description_ = long_description;
   }

   /** Get the registering class */
   virtual SmartPtr<const RegisteredCategory> RegisteringCategory() const
   {
      return ConstPtr(registering_category_);
   }

   /** Get the Option's type */
   virtual const RegisteredOptionType& Type() const
   {
      return type_;

   }
   /** Set the Option's type */
   virtual void SetType(
      const RegisteredOptionType& type
   )
   {
      type_ = type;
   }

   /** Get the advanced flag */
   virtual bool Advanced() const
   {
      return advanced_;
   }
   /** Set the advanced flag */
   virtual void SetAdvanced(
      bool advanced = true
   )
   {
      advanced_ = advanced;
   }

   /** Counter */
   virtual Index Counter() const
   {
      return counter_;
   }
   ///@}

   /** @name Get / Set methods valid for specific types
    *
    * @note The Type must be set before calling these methods.
    */
   ///@{
   /** check if the option has a lower bound
    *
    * can be called for OT_Number & OT_Integer
    */
   virtual const bool& HasLower() const
   {
      DBG_ASSERT(type_ == OT_Number || type_ == OT_Integer);
      return has_lower_;
   }

   /** check if the lower bound is strict
    *
    * can be called for OT_Number
    */
   virtual const bool& LowerStrict() const
   {
      DBG_ASSERT(type_ == OT_Number && has_lower_ == true);
      return lower_strict_;
   }

   /** get the Number version of the lower bound
    *
    * can be called for OT_Number
    */
   virtual Number LowerNumber() const
   {
      DBG_ASSERT(has_lower_ == true && type_ == OT_Number);
      return lower_;
   }

   /** set the Number version of the lower bound
    *
    * can be called for OT_Number
    */
   virtual void SetLowerNumber(
      const Number& lower,
      const bool&   strict
   )
   {
      DBG_ASSERT(type_ == OT_Number);
      lower_ = lower;
      lower_strict_ = strict, has_lower_ = true;
   }

   /** get the Integer version of the lower bound
    *
    * can be called for OT_Integer
    */
   virtual Index LowerInteger() const
   {
      DBG_ASSERT(has_lower_ == true && type_ == OT_Integer);
      return (Index) lower_;
   }

   /** set the Integer version of the lower bound
    *
    * can be called for OT_Integer
    */
   virtual void SetLowerInteger(
      const Index& lower
   )
   {
      DBG_ASSERT(type_ == OT_Integer);
      lower_ = (Number) lower;
      has_lower_ = true;
   }

   /** check if the option has an upper bound
    *
    * can be called for OT_Number & OT_Integer
    */
   virtual const bool& HasUpper() const
   {
      DBG_ASSERT(type_ == OT_Number || type_ == OT_Integer);
      return has_upper_;
   }

   /** check if the upper bound is strict
    *
    * can be called for OT_Number
    */
   virtual const bool& UpperStrict() const
   {
      DBG_ASSERT(type_ == OT_Number && has_upper_ == true);
      return upper_strict_;
   }

   /** get the Number version of the upper bound
    *
    * can be called for OT_Number
    */
   virtual Number UpperNumber() const
   {
      DBG_ASSERT(has_upper_ == true && type_ == OT_Number);
      return upper_;
   }

   /** set the Number version of the upper bound
    *
    * can be called for OT_Number
    */
   virtual void SetUpperNumber(
      const Number& upper,
      const bool&   strict
   )
   {
      DBG_ASSERT(type_ == OT_Number);
      upper_ = upper;
      upper_strict_ = strict;
      has_upper_ = true;
   }

   /** get the Integer version of the upper bound
    *
    * can be called for OT_Integer
    */
   virtual Index UpperInteger() const
   {
      DBG_ASSERT(has_upper_ == true && type_ == OT_Integer);
      return (Index) upper_;
   }

   /** set the Integer version of the upper bound
    *
    * can be called for OT_Integer
    */
   virtual void SetUpperInteger(
      const Index& upper
   )
   {
      DBG_ASSERT(type_ == OT_Integer);
      upper_ = (Number) upper;
      has_upper_ = true;
   }

   /** method to add valid string entries
    *
    * can be called for OT_String
    */
   virtual void AddValidStringSetting(
      const std::string value,
      const std::string description)
   {
      DBG_ASSERT(type_ == OT_String);
      valid_strings_.push_back(string_entry(value, description));
   }

   /** get the default as a Number
    *
    * can be called for OT_Number
    */
   virtual Number DefaultNumber() const
   {
      DBG_ASSERT(type_ == OT_Number);
      return default_number_;
   }

   /** Set the default as a Number
    *
    * can be called for OT_Number
    */
   virtual void SetDefaultNumber(
      const Number& default_value
   )
   {
      DBG_ASSERT(type_ == OT_Number);
      default_number_ = default_value;
   }

   /** get the default as an Integer
    *
    * can be called for OT_Integer
    */
   virtual Index DefaultInteger() const
   {
      DBG_ASSERT(type_ == OT_Integer);
      return (Index) default_number_;
   }

   /** Set the default as an Integer
    *
    * can be called for OT_Integer
    */
   virtual void SetDefaultInteger(
      const Index& default_value
   )
   {
      DBG_ASSERT(type_ == OT_Integer);
      default_number_ = (Number) default_value;
   }

   /** get the default as a string
    *
    * can be called for OT_String
    */
   virtual std::string DefaultString() const
   {
      DBG_ASSERT(type_ == OT_String);
      return default_string_;
   }

   /** get the default as a string, but as the index of the string in the list
    *
    *  helps map from a string to an enum
    *
    *  can be called for OT_String
    */
   virtual Index DefaultStringAsEnum() const
   {
      DBG_ASSERT(type_ == OT_String);
      return MapStringSettingToEnum(default_string_);
   }

   /** Set the default as a string
    *
    * can be called for OT_String
    */
   virtual void SetDefaultString(
      const std::string& default_value
   )
   {
      DBG_ASSERT(type_ == OT_String);
      default_string_ = default_value;
   }

   /** get the valid string settings
    *
    * can be called for OT_String
    */
   virtual std::vector<string_entry> GetValidStrings() const
   {
      DBG_ASSERT(type_ == OT_String);
      return valid_strings_;
   }

   /** Check if the Number value is a valid setting
    *
    * can be called for OT_Number
    * */
   virtual bool IsValidNumberSetting(
      const Number& value
   ) const
   {
      DBG_ASSERT(type_ == OT_Number);
      if( has_lower_ && ((lower_strict_ == true && value <= lower_) || (lower_strict_ == false && value < lower_)) )
      {
         return false;
      }
      if( has_upper_ && ((upper_strict_ == true && value >= upper_) || (upper_strict_ == false && value > upper_)) )
      {
         return false;
      }
      return true;
   }

   /** Check if the Integer value is a valid setting
    *
    * can be called for OT_Integer
    */
   virtual bool IsValidIntegerSetting(
      const Index& value
   ) const
   {
      DBG_ASSERT(type_ == OT_Integer);
      if( has_lower_ && value < lower_ )
      {
         return false;
      }
      if( has_upper_ && value > upper_ )
      {
         return false;
      }
      return true;
   }

   /** Check if the String value is a valid setting
    *
    * can be called for OT_String
    */
   virtual bool IsValidStringSetting(
      const std::string& value
   ) const;

   /** Map a user setting (allowing any case) to the case used when
    *  the setting was registered.
    */
   virtual std::string MapStringSetting(
      const std::string& value
   ) const;

   /** Map a user setting (allowing any case) to the index of the
    *  matched setting in the list of string settings.
    *
    *  Helps map a string setting to an enumeration.
    */
   virtual Index MapStringSettingToEnum(
      const std::string& value
   ) const;
   ///@}

   /** output a description of the option */
   virtual void OutputDescription(
      const Journalist& jnlst
   ) const;

   /** output a more concise version */
   virtual void OutputShortDescription(
      const Journalist& jnlst
   ) const;

   /** output a latex version */
   virtual void OutputLatexDescription(
      const Journalist& jnlst
   ) const;

   /** output a doxygen version */
   virtual void OutputDoxygenDescription(
      const Journalist& jnlst
   ) const;

private:
   std::string name_;
   std::string short_description_;
   std::string long_description_;
   SmartPtr<RegisteredCategory> registering_category_;
   RegisteredOptionType type_;
   bool advanced_;

   bool has_lower_;
   bool lower_strict_;
   Number lower_;
   bool has_upper_;
   bool upper_strict_;
   Number upper_;
   Number default_number_;

   std::vector<string_entry> valid_strings_;
   std::string default_string_;

   /** Has the information as how many-th option this one was
    *  registered. */
   const Index counter_;

   void MakeValidLatexString(
      std::string source,
      std::string& dest
   ) const;

   std::string MakeValidLatexNumber(
      Number value
   ) const;

   std::string MakeValidHTMLNumber(
      Number value
   ) const;

   /** Compare two strings and return true if they are equal (case insensitive comparison) */
   bool string_equal_insensitive(
      const std::string& s1,
      const std::string& s2
   ) const;
};

/** Class for storing registered options.
 *
 * Used for validation and documentation.
 */
class IPOPTLIB_EXPORT RegisteredOptions: public ReferencedObject
{
public:
   typedef std::map<std::string, SmartPtr<RegisteredOption> > RegOptionsList;
   typedef std::map<std::string, SmartPtr<RegisteredCategory> > RegCategoriesList;
   typedef std::set<SmartPtr<RegisteredCategory>, RegisteredCategory::ComparePriority> RegCategoriesByPriority;

   /** Constructors / Destructors */
   ///@{
   /** Default Constructor */
   RegisteredOptions()
      : next_counter_(0)
   { }

   /** Destructor */
   virtual ~RegisteredOptions()
   {
      // break circular reference between registered options and registered categories
      for( RegCategoriesList::iterator it(registered_categories_.begin()); it != registered_categories_.end(); ++it )
         it->second->regoptions_.clear();
   }
   ///@}

   DECLARE_STD_EXCEPTION(OPTION_ALREADY_REGISTERED);

   /** Methods to interact with registered options */
   ///@{
   /** set the registering class
    *
    * If nonempty name, then all subsequent options will be added with the registered category.
    * If empty name, then all subsequent options will not be added to any registered category.
    *
    * If the category doesn't exist yet, it will be created with given data.
    * If it exists already, given priority and undocumented flag are ignored.
    */
   virtual void SetRegisteringCategory(
      const std::string& registering_category,
      int                priority = 0
   )
   {
      if( registering_category.empty() )
      {
         current_registering_category_ = NULL;
         return;
      }

      SmartPtr<RegisteredCategory>& reg_categ = registered_categories_[registering_category];
      if( !IsValid(reg_categ) )
         reg_categ = new RegisteredCategory(registering_category, priority);
      current_registering_category_ = reg_categ;
   }

   /** set the registering class
    *
    * If not NULL, then all subsequent options will be added with the registered category.
    * If NULL, then all subsequent options will not be added to any registered category.
    */
   virtual void SetRegisteringCategory(
      SmartPtr<RegisteredCategory> registering_category
   )
   {
      current_registering_category_ = registering_category;
      if( !IsValid(registering_category) )
         return;

      SmartPtr<RegisteredCategory>& reg_categ = registered_categories_[registering_category->Name()];
      if( !IsValid(reg_categ) )
         reg_categ = registering_category;
      else
      {
         // if we already had a category under this name, then it should be the same as the given one
         DBG_ASSERT(reg_categ == registering_category);
      }
   }

   /** retrieve the value of the current registering category */
   virtual SmartPtr<RegisteredCategory> RegisteringCategory()
   {
      return current_registering_category_;
   }

   /** Add a Number option (with no restrictions) */
   virtual void AddNumberOption(
      const std::string& name,
      const std::string& short_description,
      Number             default_value,
      const std::string& long_description = "",
      bool               advanced = false
   );

   /** Add a Number option (with a lower bound) */
   virtual void AddLowerBoundedNumberOption(
      const std::string& name,
      const std::string& short_description,
      Number             lower,
      bool               strict,
      Number             default_value,
      const std::string& long_description = "",
      bool               advanced = false
   );

   /** Add a Number option (with a upper bound) */
   virtual void AddUpperBoundedNumberOption(
      const std::string& name,
      const std::string& short_description,
      Number             upper,
      bool               strict,
      Number             default_value,
      const std::string& long_description = "",
      bool               advanced = false
   );

   /** Add a Number option (with a both bounds) */
   virtual void AddBoundedNumberOption(
      const std::string& name,
      const std::string& short_description,
      Number             lower,
      bool               lower_strict,
      Number             upper,
      bool               upper_strict,
      Number             default_value,
      const std::string& long_description = "",
      bool               advanced = false
   );

   /** Add a Integer option (with no restrictions) */
   virtual void AddIntegerOption(
      const std::string& name,
      const std::string& short_description,
      Index              default_value,
      const std::string& long_description = "",
      bool               advanced = false
   );

   /** Add a Integer option (with a lower bound) */
   virtual void AddLowerBoundedIntegerOption(
      const std::string& name,
      const std::string& short_description,
      Index              lower,
      Index              default_value,
      const std::string& long_description = "",
      bool               advanced = false
   );

   /** Add a Integer option (with a upper bound) */
   virtual void AddUpperBoundedIntegerOption(
      const std::string& name,
      const std::string& short_description,
      Index              upper,
      Index              default_value,
      const std::string& long_description = "",
      bool               advanced = false
   );

   /** Add a Integer option (with a both bounds) */
   virtual void AddBoundedIntegerOption(
      const std::string& name,
      const std::string& short_description,
      Index              lower,
      Index              upper,
      Index              default_value,
      const std::string& long_description = "",
      bool               advanced = false
   );

   /** Add a String option (with no restrictions) */
   virtual void AddStringOption(
      const std::string&              name,
      const std::string&              short_description,
      const std::string&              default_value,
      const std::vector<std::string>& settings,
      const std::vector<std::string>& descriptions,
      const std::string&              long_description = "",
      bool               advanced = false
   );

   /** Methods that make adding string options with only a few entries easier */
   virtual void AddStringOption1(
      const std::string& name,
      const std::string& short_description,
      const std::string& default_value,
      const std::string& setting1,
      const std::string& description1,
      const std::string& long_description = "",
      bool               advanced = false
   );

   virtual void AddStringOption2(
      const std::string& name,
      const std::string& short_description,
      const std::string& default_value,
      const std::string& setting1,
      const std::string& description1,
      const std::string& setting2,
      const std::string& description2,
      const std::string& long_description = "",
      bool               advanced = false
   );

   virtual void AddStringOption3(
      const std::string& name,
      const std::string& short_description,
      const std::string& default_value,
      const std::string& setting1,
      const std::string& description1,
      const std::string& setting2,
      const std::string& description2,
      const std::string& setting3,
      const std::string& description3,
      const std::string& long_description = "",
      bool               advanced = false
   );

   virtual void AddStringOption4(
      const std::string& name,
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
      const std::string& long_description = "",
      bool               advanced = false
   );

   virtual void AddStringOption5(
      const std::string& name,
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
      const std::string& long_description = "",
      bool               advanced = false
   );

   virtual void AddStringOption6(
      const std::string& name,
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
      const std::string& long_description = "",
      bool               advanced = false
   );

   virtual void AddStringOption7(
      const std::string& name,
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
      const std::string& long_description = "",
      bool               advanced = false
   );

   virtual void AddStringOption8(
      const std::string& name,
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
      const std::string& setting8,
      const std::string& description8,
      const std::string& long_description = "",
      bool               advanced = false
   );

   virtual void AddStringOption9(
      const std::string& name,
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
      const std::string& setting8,
      const std::string& description8,
      const std::string& setting9,
      const std::string& description9,
      const std::string& long_description = "",
      bool               advanced = false
   );

   virtual void AddStringOption10(
      const std::string& name,
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
      const std::string& setting8,
      const std::string& description8,
      const std::string& setting9,
      const std::string& description9,
      const std::string& setting10,
      const std::string& description10,
      const std::string& long_description = "",
      bool               advanced = false
   );

   /** Get a registered option
    *
    * @return NULL, if the option does not exist
    */
   virtual SmartPtr<const RegisteredOption> GetOption(
      const std::string& name
   );

   /** Giving access to iteratable representation of the registered options */
   const RegOptionsList& RegisteredOptionsList() const
   {
      return registered_options_;
   }

   /** Giving access to registered categories */
   const RegCategoriesList& RegisteredCategories() const
   {
      return registered_categories_;
   }

   /** Giving access to registered categories ordered by (increasing) priority
    *
    * Result is stored in given set.
    */
   void RegisteredCategoriesByPriority(
      RegCategoriesByPriority& categories
      ) const;

   /** Output documentation in text format
    *
    * If categories is empty, then all options are printed.
    */
   virtual void OutputOptionDocumentation(
      const Journalist&             jnlst,
      const std::list<std::string>& categories = std::list<std::string>(),
      bool                          output_advanced = false
   ) const;

   /** Output documentation in Latex format to include in a latex file
    *
    * If options_to_print is empty, then all options are printed.
    */
   virtual void OutputLatexOptionDocumentation(
      const Journalist&             jnlst,
      const std::list<std::string>& options_to_print = std::list<std::string>(),
      bool                          output_advanced = false
   ) const;

   /** Output documentation in Doxygen format to include in doxygen documentation
    *
    * If options_to_print is empty, then all options are printed.
    */
   virtual void OutputDoxygenOptionDocumentation(
      const Journalist&             jnlst,
      const std::list<std::string>& options_to_print = std::list<std::string>(),
      bool                          output_advanced = false
   ) const;
   ///@}

private:
   void AddOption(
      const SmartPtr<RegisteredOption>& option
      );

   RegOptionsList registered_options_;
   RegCategoriesList registered_categories_;

   Index next_counter_;
   SmartPtr<RegisteredCategory> current_registering_category_;
};

} // namespace Ipopt

#endif
