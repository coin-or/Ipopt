// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPJOURNALIST_HPP__
#define __IPJOURNALIST_HPP__

#ifdef OLD_C_HEADERS
# include <stdarg.h>
#else
# include <cstdarg>
#endif
#include <string>
#include <vector>
#include "IpTypes.hpp"
#include "IpReferenced.hpp"

namespace Ipopt
{

  // forward declarations
  class Journal;
  class Matrix;
  class Vector;

  /**@name Journalist Enumerations. */
  //@{
  /** Print Level Enum. */
  enum EJournalLevel {
    J_NONE=0,
    J_ERROR,
    J_WARNING,
    J_SUMMARY,
    J_DETAILED,
    J_MOREDETAILED,
    J_VECTOR,
    J_MOREVECTOR,
    J_MATRIX,
    J_ALL,
    J_LAST_LEVEL
  };

  /** Category Selection Enum. */
  enum EJournalCategory {
    J_DBG=0,
    J_STATISTICS,
    J_MAIN,
    J_INITIALIZATION,
    J_BARRIER_UPDATE,
    J_SOLVE_PD_SYSTEM,
    J_FRAC_TO_BOUND,
    J_LINEAR_ALGEBRA,
    J_LINE_SEARCH,
    J_SOLUTION,
    J_LAST_CATEGORY
  };
  //@}

  /** Class responsible for all message output.
   * This class is responsible for all messaging and output.
   * The "printing" code or "author" should send ALL messages to the
   * Journalist, indicating an appropriate category and print level.
   * The journalist then decides, based on reader specified
   * acceptance criteria, which message is actually printed in which 
   * journals.
   * This allows the printing code to send everything, while the 
   * "reader" can decide what they really want to see.
   * 
   * Authors:
   * Authors use the 
   * Journals: You can add as many Journals as you like to the
   * Journalist with the AddJournal method. Each one represents 
   * a different printing location (or file).  Then, you can 
   * call the "print" methods of the Journalist to output
   * information to each of the journals.
   * 
   * Acceptance Criteria: Each print message should be flagged 
   * appropriately with an EJournalCategory and EJournalLevel.
   * 
   * The AddJournal
   * method returns a pointer to the newly created Journal object
   * (if successful) so you can set Acceptance criteria for that
   * particular location.
   * 
   */
  class Journalist : public ReferencedObject
  {
  public:
    /**@name Constructor / Desructor. */
    //@{
    /** Constructor. */
    Journalist();

    /** Destructor... */
    virtual ~Journalist();
    //@}

    /**@name Author Methods.
     * These methods are used by authoring code, or code that wants
     * to report some information.
     */
    //@{
    /** Method to print a formatted string */
    void Printf(EJournalLevel level, EJournalCategory category,
                const char* format, ...) const;

    /** Method to print a formatted string with indentation */
    void PrintfIndented(EJournalLevel level,
                        EJournalCategory category,
                        Index indent_level,
                        const char* format, ...) const;

    /** Method for printing a vector.  This calls the Vector print
     *  methods for each appropriate Journal.
     */
    void PrintVector(EJournalLevel level,
                     EJournalCategory category,
                     const std::string& name,
                     const Vector& vector,
                     Index indent=0,
                     std::string prefix="") const;

    /** Method for printing a matrix.  This calls the Matrix print
     *  methods for each appropriate Journal.
     */
    void PrintMatrix(EJournalLevel level,
                     EJournalCategory category,
                     const std::string& name,
                     const Matrix& matrix,
                     Index indent=0,
                     std::string prefix="") const;

    /** Method to print a formatted string
     * using the va_list argument. */
    void VPrintf(EJournalLevel level,
                 EJournalCategory category,
                 const char* pformat,
                 va_list ap) const;

    /** Method to print a formatted string with indentation,
     * using the va_list argument. */
    void VPrintfIndented(EJournalLevel level,
                         EJournalCategory category,
                         Index indent_level,
                         const char* pformat,
                         va_list ap) const;

    /** Method that returns true if there is a Journal that would
     *  write output for the given JournalLevel and JournalCategory.
     *  This is useful if expensive computation would be required for
     *  a particular output.  The author code can check with this
     *  method if the computations are indeed required.
     */
    bool ProduceOutput(EJournalLevel level,
                       EJournalCategory category) const;

    //@}

    /**@name Reader Methods.
     * These methods are used by the reader. The reader will setup the 
     * journalist with each output file and the acceptance
     * criteria for that file.
     *
     * Use these methods to setup the wanted journals (files).
     * These are the internal objects that keep track of the print levels 
     * for each category. Then use the internal Journal objects to
     * set specific print levels for each category (or keep defaults).
     *  
     * Note: the lifetime of these internal objects is controlled by
     * the Journalist class. Do not try to delete and do not try to keep the
     * pointer outside of the current scope. 
     */
    //@{
    /** Add a new journal.  The location_name is a string identifier,
     *  which can be used to obtain the pointer to the new Journal at
     *  a later point using the GetJournal method.  fname is the name
     *  of the * file to which this Journal corresponds.  Use
     *  fname="stdout" * for stdout, and use fname="stderr" for
     *  stderr.  This method * returns the Journal pointer so you can
     *  set specific acceptance * criteria.  It returns NULL if there
     *  was a problem creating a * new Journal.  The default_level is
     *  used to initialize the * printing level for all categories.
     */
    Journal* AddJournal(
      const std::string& location_name,    /** identifier */
      const std::string& fname,
      EJournalLevel default_level = J_WARNING
    );

    /** Get an existing journal.  You can use this method to change
     *  the acceptance criteria at runtime.
     */
    Journal* GetJournal(const std::string& location_name);
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    Journalist(const Journalist&);

    /** Overloaded Equals Operator */
    void operator=(const Journalist&);
    //@}

    //** Private Data Members. */
    //@{
    std::vector<Journal*> journals_;
    //@}
  };

  /** Journal class (part of the Journalist implementation.).  This
   * class stores the specific information for one journal, the file
   * location (stdout, stderr, or disk), and the acceptance criteria.
   */
  class Journal
  {
  public:
    /** Constructor. */
    Journal(const std::string& name, EJournalLevel default_level);

    /** Destructor. */
    ~Journal();

    /** Get the name of the Journal */
    std::string Name();

    /** Set the print level for a particular category. */
    void SetPrintLevel(
      EJournalCategory category, EJournalLevel level
    );

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
    Journal();

    /** Copy Constructor */
    Journal(const Journal&);

    /** Overloaded Equals Operator */
    void operator=(const Journal&);
    //@}

    /** Name of the output location */
    std::string name_;

    /** vector of integers indicating the level for each category */
    Index print_levels_[J_LAST_CATEGORY];

    /** FILE pointer for the output destination */
    FILE* file_;

    /**@name Private methods for the Journalist (friend) to call*/
    //@{
    /** Open a new file for the output location.
     *  Special Names: stdout means stdout,
     *               : stderr means stderr.
     *
     *  Return code is false only if the file with the given name
     *  could not be opened.
     */
    bool Open(const char* fname);

    /** Ask if a particular print level/category is accepted by the
     * journal.
     */
    bool IsAccepted(
      EJournalCategory category, EJournalLevel level
    ) const;

    /** Print to the designated output location */
    void Print(const char* str);

    /** Printf to the designated output location */
    void Printf(const char* pformat, va_list ap);

    /** Print vector to the designated output location */
    void PrintVector(std::string name, const Vector& vector, Index indent, std::string prefix);

    /** Print matrix to the designated output location */
    void PrintMatrix(const std::string name, const Matrix& matrix, Index indent, std::string prefix);
    //@}

    friend class Journalist;
  };
}

#endif
