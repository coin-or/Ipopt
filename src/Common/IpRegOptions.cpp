// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-06-18

#include "IpoptConfig.h"
#include "IpRegOptions.hpp"
#include "IpOptionsList.hpp"

#include <cstdio>
#include <cctype>

namespace Ipopt
{

void RegisteredOption::OutputDescription(
   const Journalist& jnlst
) const
{
   std::string type_str = "Unknown";
   if( type_ == OT_Number )
   {
      type_str = "Real Number";
   }
   else if( type_ == OT_Integer )
   {
      type_str = "Integer";
   }
   else if( type_ == OT_String )
   {
      type_str = "String";
   }

   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                "\n### %s (%s) %s ###\nCategory: %s\nDescription: %s\n", name_.c_str(), type_str.c_str(),
                advanced_ ? "(advanced)" : "",
                IsValid(registering_category_) ? registering_category_->Name().c_str() : "n/a", short_description_.c_str());

   if( type_ == OT_Number )
   {
      if( has_lower_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%g", lower_);
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "-inf");
      }

      if( lower_strict_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " < ");
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <= ");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "(%g)", default_number_);

      if( has_upper_ && upper_strict_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " < ");
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <= ");
      }

      if( has_upper_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%g\n", upper_);
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "+inf\n");
      }
   }
   else if( type_ == OT_Integer )
   {
      if( has_lower_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%" IPOPT_INDEX_FORMAT, (Index) lower_);
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "-inf");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   " <= (%" IPOPT_INDEX_FORMAT ") <= ", (Index) default_number_);

      if( has_upper_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%" IPOPT_INDEX_FORMAT "\n", (Index) upper_);
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "+inf\n");
      }
   }
   else if( type_ == OT_String )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "Valid Settings:\n");
      for( std::vector<string_entry>::const_iterator i = valid_strings_.begin(); i != valid_strings_.end(); ++i )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "\t%s (%s)\n", (*i).value_.c_str(), (*i).description_.c_str());
      }
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "Default: \"%s\"\n", default_string_.c_str());
   }
}

void RegisteredOption::OutputLatexDescription(
   const Journalist& jnlst
) const
{
   std::string latex_name;
   MakeValidLatexString(name_, latex_name);
   std::string latex_desc;
   MakeValidLatexString(short_description_, latex_desc);
   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                "\\paragraph{%s:}\\label{opt:%s} ", latex_name.c_str(), name_.c_str());
   if( advanced_ )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "(advanced) ");
   }
   if( short_description_.length() == 0 )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "~");
   }
   else
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%s",
                   latex_desc.c_str());
   }
   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                " \\\\\n");

   //    Index length = name_.length() + short_description_.length();
   //    DBG_ASSERT(length <= 80);
   //    jnlst.PrintStringOverLines(J_SUMMARY, J_DOCUMENTATION, 0, 50,
   //                               latex_desc.c_str());

   if( long_description_ != "" )
   {
      latex_desc = "";
      MakeValidLatexString(long_description_,
                           latex_desc);
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   " %s", latex_desc.c_str());
   }

   if( type_ == OT_Number )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   " The valid range for this real option is \n$");
      std::string buff;
      if( has_lower_ )
      {
         buff = MakeValidLatexNumber(lower_);
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%s", buff.c_str());
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "{\\tt -inf}");
      }

      if( has_lower_ && !lower_strict_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " \\le ");
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <  ");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "{\\tt %s }", latex_name.c_str());

      if( has_upper_ && !upper_strict_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " \\le ");
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <  ");
      }

      if( has_upper_ )
      {
         buff = MakeValidLatexNumber(upper_);
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%s", buff.c_str());
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "{\\tt +inf}");
      }

      buff = MakeValidLatexNumber(default_number_);
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "$\nand its default value is $%s$.\n\n", buff.c_str());

   }
   else if( type_ == OT_Integer )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   " The valid range for this integer option is\n$");
      if( has_lower_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%" IPOPT_INDEX_FORMAT " \\le ", (Index) lower_);
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%s <  ", "{\\tt -inf}");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "{\\tt %s }", latex_name.c_str());

      if( has_upper_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " \\le %" IPOPT_INDEX_FORMAT "", (Index) upper_);
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <  %s", "{\\tt +inf}");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "$\nand its default value is $%" IPOPT_INDEX_FORMAT "$.\n\n", (Index) default_number_);
   }
   else if( type_ == OT_String )
   {
      std::string buff;
      MakeValidLatexString(default_string_, buff);
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   " The default value for this string option is \"%s\".\n", buff.c_str());

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "\\\\ \nPossible values:\n");
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "\\begin{itemize}\n");
      for( std::vector<string_entry>::const_iterator i = valid_strings_.begin(); i != valid_strings_.end(); ++i )
      {
         std::string latex_value;
         MakeValidLatexString((*i).value_, latex_value);
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "   \\item %s", latex_value.c_str());

         if( (*i).description_.length() > 0 )
         {
            MakeValidLatexString((*i).description_, latex_desc);
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                         ": %s", latex_desc.c_str());
         }
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "\n");
      }
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "\\end{itemize}\n");
   }
   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, ""
                "\n");
}

void RegisteredOption::MakeValidLatexString(
   const std::string& source,
   std::string&       dest
) const
{
   for( std::string::const_iterator c = source.begin(); c != source.end(); ++c )
   {
      if( *c == '_' )
      {
         dest.append("\\_");
      }
      else if( *c == '^' )
      {
         dest.append("\\^");
      }
      else
      {
         dest += *c;
      }
   }
}

std::string RegisteredOption::MakeValidLatexNumber(
   Number value
) const
{
   char buffer[256];
   Snprintf(buffer, 255, "%g", value);
   std::string source = buffer;
   std::string dest;

   bool found_e = false;
   for( std::string::iterator c = source.begin(); c != source.end(); ++c )
   {
      if( *c == 'e' )
      {
         found_e = true;
         dest.append(" \\cdot 10^{");
      }
      else
      {
         dest += *c;
      }
   }
   if( found_e )
   {
      dest.append("}");
   }

   return dest;
}

void RegisteredOption::OutputDoxygenDescription(
   const Journalist& jnlst
) const
{
   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                "\\anchor OPT_%s\n<strong>%s</strong>", name_.c_str(), name_.c_str());
   if( advanced_ )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " (<em>advanced</em>)");
   }

   if( short_description_.length() > 0 )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   ": %s", short_description_.c_str());
   }
   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                "\n<blockquote>\n");

   if( long_description_ != "" )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   " %s", long_description_.c_str());
   }

   if( type_ == OT_Number )
   {
      std::string buff;
      if( has_lower_ || has_upper_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " The valid range for this real option is ");
         if( has_lower_ )
         {
            buff = MakeValidHTMLNumber(lower_);
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%s",
                         buff.c_str());

            if( !lower_strict_ )
            {
               jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                            " &le; ");
            }
            else
            {
               jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                            " < ");
            }
         }
         //else
         //{
         //   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
         //       "-&infin; < ");
         //}

         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%s",
                      name_.c_str());

         if( has_upper_ )
         {
            if( !upper_strict_ )
            {
               jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                            " &le; ");
            }
            else
            {
               jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                            " < ");
            }

            buff = MakeValidHTMLNumber(upper_);
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%s",
                         buff.c_str());
         }
         // else
         // {
         //   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
         //       "< &infin;");
         //}
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " The valid range for this real option is unrestricted");
      }

      buff = MakeValidHTMLNumber(default_number_);
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   " and its default value is %s.\n", buff.c_str());

   }
   else if( type_ == OT_Integer )
   {
      if( has_lower_ || has_upper_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " The valid range for this integer option is ");
         if( has_lower_ )
         {
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                         "%" IPOPT_INDEX_FORMAT " &le; ", (Index) lower_);
         }
         //else
         //{
         //   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
         //        "-&infin; < ");
         //}

         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "%s",
                      name_.c_str());

         if( has_upper_ )
         {
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                         " &le; %" IPOPT_INDEX_FORMAT "", (Index) upper_);
         }
         //else
         //{
         //   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
         //       " < &infin;");
         //}
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " The valid range for this integer option is unrestricted");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   " and its default value is %" IPOPT_INDEX_FORMAT ".\n", (Index) default_number_);
   }
   else if( type_ == OT_String )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   " The default value for this string option is \"%s\".\n", default_string_.c_str());

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "\nPossible values:");

      bool havedescr = false;
      for( std::vector<string_entry>::const_iterator i = valid_strings_.begin(); i != valid_strings_.end() && !havedescr; ++i )
         if( (*i).description_.length() > 0 )
         {
            havedescr = true;
         }

      if( havedescr )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n");
         for( std::vector<string_entry>::const_iterator i = valid_strings_.begin(); i != valid_strings_.end(); ++i )
         {
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " - %s", i->value_.c_str());
            if( (*i).description_.length() > 0 )
            {
               jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, ": %s", i->description_.c_str());
            }
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n");
         }
      }
      else
      {
         for( std::vector<string_entry>::const_iterator i = valid_strings_.begin(); i != valid_strings_.end(); ++i )
         {
            if( i != valid_strings_.begin() )
            {
               jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, ",");
            }
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, " %s", i->value_.c_str());
         }
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n");
      }

      /*
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
          "Possible values:\n");
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
           "|Value|Description|\n");
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
           "|:----|:----------|\n");
      for( std::vector<string_entry>::const_iterator i = valid_strings_.begin(); i != valid_strings_.end(); i++ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
              "|%s|%s|\n", i->value_.c_str(), i->description_.c_str());
      }
      */
   }
   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                "</blockquote>\n\n");
}

std::string RegisteredOption::MakeValidHTMLNumber(
   Number value
) const
{
   char buffer[256];
   Snprintf(buffer, 255, "%g", value);
   std::string source = buffer;
   std::string dest;

   bool found_e = false;
   for( std::string::iterator c = source.begin(); c != source.end(); ++c )
   {
      if( *c == 'e' )
      {
         found_e = true;
         if( dest == "1" )
         {
            dest = "";
         }
         else if( dest == "-1" )
         {
            dest = "-";
         }
         else
         {
            dest.append(" &middot; ");
         }
         dest += "10<sup>";
      }
      else
      {
         dest += *c;
      }
   }
   if( found_e )
   {
      dest.append("</sup>");
   }

   return dest;
}

void RegisteredOption::OutputShortDescription(
   const Journalist& jnlst
) const
{
   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                "%-30s", name_.c_str());

   if( type_ == OT_Number )
   {
      if( has_lower_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%10g", lower_);
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%10s", "-inf");
      }

      if( has_lower_ && !lower_strict_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <= ");
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <  ");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "(%11g)", default_number_);

      if( has_upper_ && !upper_strict_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <= ");
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <  ");
      }

      if( has_upper_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%-10g\n", upper_);
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%-10s\n", "+inf");
      }
   }
   else if( type_ == OT_Integer )
   {
      if( has_lower_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%10" IPOPT_INDEX_FORMAT " <= ", (Index) lower_);
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "%10s <  ", "-inf");
      }

      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "(%11" IPOPT_INDEX_FORMAT ")", (Index) default_number_);

      if( has_upper_ )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <= %-10" IPOPT_INDEX_FORMAT "\n", (Index) upper_);
      }
      else
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      " <  %-10s\n", "+inf");
      }
   }
   else if( type_ == OT_String )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "(\"%s\")\n", default_string_.c_str());
   }
   if( advanced_ )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "   Advanced option for expert users.\n");
   }
   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                "   ");
   jnlst.PrintStringOverLines(J_SUMMARY, J_DOCUMENTATION, 3, 76, short_description_);
   if( long_description_ != "" )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "\n     ");
      jnlst.PrintStringOverLines(J_SUMMARY, J_DOCUMENTATION, 5, 74, long_description_);
   }
   if( type_ == OT_String )
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "\n   Possible values:\n");
      for( std::vector<string_entry>::const_iterator i = valid_strings_.begin(); i != valid_strings_.end(); ++i )
      {
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "    - %-23s", (*i).value_.c_str());

         if( (*i).description_.length() > 0 )
         {
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                         " [");
            jnlst.PrintStringOverLines(J_SUMMARY, J_DOCUMENTATION, 31, 48, (*i).description_);
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                         "]");
         }
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                      "\n");
      }
   }
   else
   {
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                   "\n");
   }
   jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                "\n");
}

bool RegisteredOption::IsValidStringSetting(
   const std::string& value
) const
{
   DBG_ASSERT(type_ == OT_String);

   for( std::vector<string_entry>::const_iterator i = valid_strings_.begin(); i != valid_strings_.end(); ++i )
   {
      if( i->value_ == "*" || string_equal_insensitive(i->value_, value) )
      {
         return true;
      }
   }
   return false;
}

std::string RegisteredOption::MapStringSetting(
   const std::string& value
) const
{
   DBG_ASSERT(type_ == OT_String);

   std::string matched_setting = "";

   for( std::vector<string_entry>::const_iterator i = valid_strings_.begin(); i != valid_strings_.end(); ++i )
   {
      if( i->value_ == "*" )
      {
         matched_setting = value;
      }
      else if( string_equal_insensitive(i->value_, value) )
      {
         matched_setting = i->value_;
      }
   }
   return matched_setting;
}

Index RegisteredOption::MapStringSettingToEnum(
   const std::string& value
) const
{
   DBG_ASSERT(type_ == OT_String);

   Index matched_setting = -1;

   Index cnt = 0;
   for( std::vector<string_entry>::const_iterator i = valid_strings_.begin(); i != valid_strings_.end(); ++i )
   {
      ASSERT_EXCEPTION(i->value_ != "*", IpoptException, "Cannot map a wildcard setting to an enumeration");
      if( string_equal_insensitive(i->value_, value) )
      {
         matched_setting = cnt;
         break;
      }
      cnt++;
   }

   ASSERT_EXCEPTION(matched_setting != -1, ERROR_CONVERTING_STRING_TO_ENUM,
                    std::string("Could not find a match for setting ") + value + " in option: " + name_);
   return matched_setting;
}

bool RegisteredOption::string_equal_insensitive(
   const std::string& s1,
   const std::string& s2
) const
{
   using namespace std;

   if( s1.size() != s2.size() )
   {
      return false;
   }

   string::const_iterator i1 = s1.begin();
   string::const_iterator i2 = s2.begin();

   while( i1 != s1.end() )
   {
      if( toupper(*i1) != toupper(*i2) )
      {
         return false;
      }
      ++i1;
      ++i2;
   }
   return true;
}

void RegisteredOptions::SetRegisteringCategory(
   const std::string& registering_category,
   int                priority
)
{
   if( registering_category.empty() )
   {
      current_registering_category_ = NULL;
      return;
   }

   SmartPtr<RegisteredCategory>& reg_categ = registered_categories_[registering_category];
   if( !IsValid(reg_categ) )
   {
      reg_categ = new RegisteredCategory(registering_category, priority);
   }
   current_registering_category_ = reg_categ;
}

void RegisteredOptions::SetRegisteringCategory(
   SmartPtr<RegisteredCategory> registering_category
)
{
   current_registering_category_ = registering_category;
   if( !IsValid(registering_category) )
   {
      return;
   }

   SmartPtr<RegisteredCategory>& reg_categ = registered_categories_[registering_category->Name()];
   if( !IsValid(reg_categ) )
   {
      reg_categ = registering_category;
   }
   else
   {
      // if we already had a category under this name, then it should be the same as the given one
      DBG_ASSERT(reg_categ == registering_category);
   }
}

void RegisteredOptions::AddOption(
   const SmartPtr<RegisteredOption>& option
)
{
   ASSERT_EXCEPTION(registered_options_.find(option->Name()) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
                    std::string("The option: ") + option->Name() + " has already been registered by someone else");
   registered_options_[option->Name()] = option;

   if( IsValid(option->registering_category_) )
   {
      option->registering_category_->regoptions_.push_back(option);
   }
}

void RegisteredOptions::AddNumberOption(
   const std::string& name,
   const std::string& short_description,
   Number             default_value,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_Number);
   option->SetDefaultNumber(default_value);
   AddOption(option);
}

void RegisteredOptions::AddLowerBoundedNumberOption(
   const std::string& name,
   const std::string& short_description,
   Number             lower,
   bool               strict,
   Number             default_value,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_Number);
   option->SetDefaultNumber(default_value);
   option->SetLowerNumber(lower, strict);
   AddOption(option);
}

void RegisteredOptions::AddUpperBoundedNumberOption(
   const std::string& name,
   const std::string& short_description,
   Number             upper,
   bool               strict,
   Number             default_value,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_Number);
   option->SetDefaultNumber(default_value);
   option->SetUpperNumber(upper, strict);
   AddOption(option);
}

void RegisteredOptions::AddBoundedNumberOption(
   const std::string& name,
   const std::string& short_description,
   Number             lower,
   bool               lower_strict,
   Number             upper,
   bool               upper_strict,
   Number             default_value,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_Number);
   option->SetDefaultNumber(default_value);
   option->SetLowerNumber(lower, lower_strict);
   option->SetUpperNumber(upper, upper_strict);
   AddOption(option);
}

void RegisteredOptions::AddIntegerOption(
   const std::string& name,
   const std::string& short_description,
   Index              default_value,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_Integer);
   option->SetDefaultInteger(default_value);
   AddOption(option);
}

void RegisteredOptions::AddLowerBoundedIntegerOption(
   const std::string& name,
   const std::string& short_description,
   Index              lower,
   Index              default_value,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_Integer);
   option->SetDefaultInteger(default_value);
   option->SetLowerInteger(lower);
   AddOption(option);
}

void RegisteredOptions::AddUpperBoundedIntegerOption(
   const std::string& name,
   const std::string& short_description,
   Index              upper,
   Index              default_value,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_Integer);
   option->SetDefaultInteger(default_value);
   option->SetUpperInteger(upper);
   AddOption(option);
}

void RegisteredOptions::AddBoundedIntegerOption(
   const std::string& name,
   const std::string& short_description,
   Index              lower,
   Index              upper,
   Index              default_value,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_Integer);
   option->SetDefaultInteger(default_value);
   option->SetLowerInteger(lower);
   option->SetUpperInteger(upper);
   AddOption(option);
}

void RegisteredOptions::AddStringOption(
   const std::string&              name,
   const std::string&              short_description,
   const std::string&              default_value,
   const std::vector<std::string>& settings,
   const std::vector<std::string>& descriptions,
   const std::string&              long_description,
   bool                            advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   DBG_ASSERT(settings.size() == descriptions.size());
   for( size_t i = 0; i < settings.size(); i++ )
   {
      option->AddValidStringSetting(settings[i], descriptions[i]);
   }
   AddOption(option);
}

void RegisteredOptions::AddStringOption1(
   const std::string& name,
   const std::string& short_description,
   const std::string& default_value,
   const std::string& setting1,
   const std::string& description1,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   option->AddValidStringSetting(setting1, description1);
   AddOption(option);
}

void RegisteredOptions::AddStringOption2(
   const std::string& name,
   const std::string& short_description,
   const std::string& default_value,
   const std::string& setting1,
   const std::string& description1,
   const std::string& setting2,
   const std::string& description2,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   option->AddValidStringSetting(setting1, description1);
   option->AddValidStringSetting(setting2, description2);
   AddOption(option);
}

void RegisteredOptions::AddStringOption3(
   const std::string& name,
   const std::string& short_description,
   const std::string& default_value,
   const std::string& setting1,
   const std::string& description1,
   const std::string& setting2,
   const std::string& description2,
   const std::string& setting3,
   const std::string& description3,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   option->AddValidStringSetting(setting1, description1);
   option->AddValidStringSetting(setting2, description2);
   option->AddValidStringSetting(setting3, description3);
   AddOption(option);
}

void RegisteredOptions::AddStringOption4(
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
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   option->AddValidStringSetting(setting1, description1);
   option->AddValidStringSetting(setting2, description2);
   option->AddValidStringSetting(setting3, description3);
   option->AddValidStringSetting(setting4, description4);
   AddOption(option);
}

void RegisteredOptions::AddStringOption5(
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
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   option->AddValidStringSetting(setting1, description1);
   option->AddValidStringSetting(setting2, description2);
   option->AddValidStringSetting(setting3, description3);
   option->AddValidStringSetting(setting4, description4);
   option->AddValidStringSetting(setting5, description5);
   AddOption(option);
}

void RegisteredOptions::AddStringOption6(
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
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   option->AddValidStringSetting(setting1, description1);
   option->AddValidStringSetting(setting2, description2);
   option->AddValidStringSetting(setting3, description3);
   option->AddValidStringSetting(setting4, description4);
   option->AddValidStringSetting(setting5, description5);
   option->AddValidStringSetting(setting6, description6);
   AddOption(option);
}

void RegisteredOptions::AddStringOption7(
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
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   option->AddValidStringSetting(setting1, description1);
   option->AddValidStringSetting(setting2, description2);
   option->AddValidStringSetting(setting3, description3);
   option->AddValidStringSetting(setting4, description4);
   option->AddValidStringSetting(setting5, description5);
   option->AddValidStringSetting(setting6, description6);
   option->AddValidStringSetting(setting7, description7);
   AddOption(option);
}

void RegisteredOptions::AddStringOption8(
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
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   option->AddValidStringSetting(setting1, description1);
   option->AddValidStringSetting(setting2, description2);
   option->AddValidStringSetting(setting3, description3);
   option->AddValidStringSetting(setting4, description4);
   option->AddValidStringSetting(setting5, description5);
   option->AddValidStringSetting(setting6, description6);
   option->AddValidStringSetting(setting7, description7);
   option->AddValidStringSetting(setting8, description8);
   AddOption(option);
}

void RegisteredOptions::AddStringOption9(
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
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   option->AddValidStringSetting(setting1, description1);
   option->AddValidStringSetting(setting2, description2);
   option->AddValidStringSetting(setting3, description3);
   option->AddValidStringSetting(setting4, description4);
   option->AddValidStringSetting(setting5, description5);
   option->AddValidStringSetting(setting6, description6);
   option->AddValidStringSetting(setting7, description7);
   option->AddValidStringSetting(setting8, description8);
   option->AddValidStringSetting(setting9, description9);
   AddOption(option);
}

void RegisteredOptions::AddStringOption10(
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
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value);
   option->AddValidStringSetting(setting1, description1);
   option->AddValidStringSetting(setting2, description2);
   option->AddValidStringSetting(setting3, description3);
   option->AddValidStringSetting(setting4, description4);
   option->AddValidStringSetting(setting5, description5);
   option->AddValidStringSetting(setting6, description6);
   option->AddValidStringSetting(setting7, description7);
   option->AddValidStringSetting(setting8, description8);
   option->AddValidStringSetting(setting9, description9);
   option->AddValidStringSetting(setting10, description10);
   AddOption(option);
}

/** Create a string value with two possible settings: yes and no */
void RegisteredOptions::AddBoolOption(
   const std::string& name,
   const std::string& short_description,
   bool               default_value,
   const std::string& long_description,
   bool               advanced
)
{
   SmartPtr<RegisteredOption> option = new RegisteredOption(name, short_description, long_description,
         current_registering_category_, next_counter_++, advanced);
   option->SetType(OT_String);
   option->SetDefaultString(default_value ? "yes" : "no");
   option->AddValidStringSetting("yes", "");
   option->AddValidStringSetting("no", "");
   AddOption(option);
}

SmartPtr<const RegisteredOption> RegisteredOptions::GetOption(
   const std::string& name
)
{
   std::string tag_only = name;
   std::string::size_type pos = name.rfind(".", name.length());
   if( pos != std::string::npos )
   {
      tag_only = name.substr(pos + 1, name.length() - pos);
   }
   SmartPtr<const RegisteredOption> option;
   std::map<std::string, SmartPtr<RegisteredOption> >::iterator reg_option = registered_options_.find(tag_only);
   if( reg_option == registered_options_.end() )
   {
      option = NULL;
   }
   else
   {
      option = ConstPtr(reg_option->second);
   }

   return option;
}

/** Giving access to registered categories ordered by priority (decreasing) */
void RegisteredOptions::RegisteredCategoriesByPriority(
   RegCategoriesByPriority& categories
) const
{
   for( RegCategoriesList::const_iterator it = registered_categories_.begin(); it != registered_categories_.end(); ++it )
   {
      categories.insert(it->second);
   }
}

/** Output documentation
 *
 * Format is decided according to print_options_mode parameter.
 */
void RegisteredOptions::OutputOptionDocumentation(
   const Journalist&             jnlst,
   SmartPtr<OptionsList>         options,
   int                           minpriority
) const
{
   OutputMode printmode;
   Index enum_int;
   options->GetEnumValue("print_options_mode", enum_int, "");
   printmode = OutputMode(enum_int);

   bool printadvanced;
   options->GetBoolValue("print_advanced_options", printadvanced, "");

   RegCategoriesByPriority cats;
   RegisteredCategoriesByPriority(cats);
   for( RegCategoriesByPriority::const_iterator cat_it = cats.begin(); cat_it != cats.end(); ++cat_it )
   {
      if( (*cat_it)->Priority() < minpriority )
      {
         break;
      }

      bool firstopt = true;
      for( std::list<SmartPtr<RegisteredOption> >::const_iterator opt_it = (*cat_it)->RegisteredOptions().begin(); opt_it != (*cat_it)->RegisteredOptions().end(); ++opt_it )
      {
         if( !printadvanced && (*opt_it)->Advanced() )
         {
            continue;
         }

         if( firstopt )
         {
            const std::string& catname = (*cat_it)->Name();
            switch( printmode )
            {
               case OUTPUTTEXT :
               {
                  jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n### %s ###\n\n", catname.c_str());
                  break;
               }

               case OUTPUTLATEX:
               {
                  jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\\subsection{%s}\n\n", catname.c_str());
                  break;
               }

               case OUTPUTDOXYGEN:
               {
                  std::string anchorname = catname;
                  for( std::string::iterator it = anchorname.begin(); it != anchorname.end(); ++it )
                     if( !isalnum(*it) )
                     {
                        *it = '_';
                     }

                  jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\\subsection OPT_%s %s\n\n", anchorname.c_str(), catname.c_str());
               }
            }

            firstopt = false;
         }

         switch( printmode )
         {
            case OUTPUTTEXT :
               (*opt_it)->OutputShortDescription(jnlst);
               break;

            case OUTPUTLATEX:
               (*opt_it)->OutputLatexDescription(jnlst);
               break;

            case OUTPUTDOXYGEN:
               (*opt_it)->OutputDoxygenDescription(jnlst);
               break;
         }
      }
      jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n");
   }
}

void RegisteredOptions::OutputOptionDocumentation(
   const Journalist&             jnlst,
   const std::list<std::string>& categories
) const
{
   if( !categories.empty() )
   {
      for( std::list<std::string>::const_iterator i = categories.begin(); i != categories.end(); ++i )
      {
         RegCategoriesList::const_iterator cat_it = registered_categories_.find(*i);
         // skip nonexisting category
         if( cat_it == registered_categories_.end() )
         {
            continue;
         }

         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n### %s ###\n\n", i->c_str());

         for( std::list<SmartPtr<RegisteredOption> >::const_iterator opt_it = cat_it->second->RegisteredOptions().begin(); opt_it != cat_it->second->RegisteredOptions().end(); ++opt_it )
         {
            (*opt_it)->OutputShortDescription(jnlst);
         }
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n");
      }
   }
   else
   {
      for( RegCategoriesList::const_iterator cat_it = registered_categories_.begin(); cat_it != registered_categories_.end(); ++cat_it )
      {
         if( cat_it->second->Priority() < 0 )
         {
            continue;
         }

         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n### %s ###\n\n", cat_it->first.c_str());

         for( std::list<SmartPtr<RegisteredOption> >::const_iterator opt_it = cat_it->second->RegisteredOptions().begin(); opt_it != cat_it->second->RegisteredOptions().end(); ++opt_it )
         {
            if( (*opt_it)->Advanced() )
            {
               continue;
            }

            (*opt_it)->OutputShortDescription(jnlst);
         }
         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n");
      }
   }
}

void RegisteredOptions::OutputLatexOptionDocumentation(
   const Journalist&             jnlst,
   const std::list<std::string>& options_to_print
) const
{
   if( !options_to_print.empty() )
   {
      for( std::list<std::string>::const_iterator coption = options_to_print.begin(); coption != options_to_print.end(); ++coption )
      {
         if( coption->c_str()[0] == '#' )
         {
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\\subsection{%s}\n\n", &coption->c_str()[1]);
            continue;
         }

         SmartPtr<RegisteredOption> option = registered_options_.at(*coption);
         DBG_ASSERT(IsValid(option));

         option->OutputLatexDescription(jnlst);
      }
   }
   else
   {
      RegCategoriesByPriority cats;
      RegisteredCategoriesByPriority(cats);
      for( RegCategoriesByPriority::const_iterator cat_it = cats.begin(); cat_it != cats.end(); ++cat_it )
      {
         if( (*cat_it)->Priority() < 0 )
         {
            break;
         }

         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\\subsection{%s}\n\n", (*cat_it)->Name().c_str());

         for( std::list<SmartPtr<RegisteredOption> >::const_iterator opt_it = (*cat_it)->RegisteredOptions().begin(); opt_it != (*cat_it)->RegisteredOptions().end(); ++opt_it )
         {
            if( (*opt_it)->Advanced() )
            {
               continue;
            }

            (*opt_it)->OutputLatexDescription(jnlst);
         }
      }
   }
}

void RegisteredOptions::OutputDoxygenOptionDocumentation(
   const Journalist&             jnlst,
   const std::list<std::string>& options_to_print
) const
{
   if( !options_to_print.empty() )
   {
      for( std::list<std::string>::const_iterator coption = options_to_print.begin(); coption != options_to_print.end(); ++coption )
      {
         if( (*coption)[0] == '#' )
         {
            std::string anchorname = &coption->c_str()[1];
            for( std::string::iterator it = anchorname.begin(); it != anchorname.end(); ++it )
               if( !isalnum(*it) )
               {
                  *it = '_';
               }
            jnlst.Printf(J_SUMMARY, J_DOCUMENTATION,
                         "\\subsection OPT_%s %s\n\n", anchorname.c_str(), &coption->c_str()[1]);

            continue;
         }

         SmartPtr<RegisteredOption> option = registered_options_.at(*coption);
         DBG_ASSERT(IsValid(option));

         option->OutputDoxygenDescription(jnlst);
      }
   }
   else
   {
      RegCategoriesByPriority cats;
      RegisteredCategoriesByPriority(cats);
      for( RegCategoriesByPriority::const_iterator cat_it = cats.begin(); cat_it != cats.end(); ++cat_it )
      {
         if( (*cat_it)->Priority() < 0 )
         {
            break;
         }

         std::string anchorname = (*cat_it)->Name();
         for( std::string::iterator it = anchorname.begin(); it != anchorname.end(); ++it )
            if( !isalnum(*it) )
            {
               *it = '_';
            }

         jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\\subsection OPT_%s %s\n\n", anchorname.c_str(), (*cat_it)->Name().c_str());

         for( std::list<SmartPtr<RegisteredOption> >::const_iterator opt_it = (*cat_it)->RegisteredOptions().begin(); opt_it != (*cat_it)->RegisteredOptions().end(); ++opt_it )
         {
            if( (*opt_it)->Advanced() )
            {
               continue;
            }

            (*opt_it)->OutputDoxygenDescription(jnlst);
         }
      }
   }
}

void RegisteredOptions::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   roptions->SetRegisteringCategory("Output");
   roptions->AddStringOption3(
      "print_options_mode",
      "format in which to print options documentation",
      "text",
      "text", "Ordinary text",
      "latex", "LaTeX formatted",
      "doxygen", "Doxygen (markdown) formatted");

   roptions->AddBoolOption(
      "print_advanced_options",
      "whether to print also advanced options",
      false,
      "",
      true);

}

} // namespace Ipopt
