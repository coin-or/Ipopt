// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpJournalist.hpp"
#include "IpDebug.hpp"

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

#ifdef HAVE_CSTRING
# include <cstring>
#else
# ifdef HAVE_STRING_H
#  include <string.h>
# else
#  error "don't have header file for string"
# endif
#endif

#ifdef HAVE_MPI
# include "IpMpi.hpp"
#endif

namespace Ipopt
{

  Journalist::Journalist()
      :
      my_rank_(0),
      collecting_output_(false)
  {
#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
#endif
  }

  Journalist::~Journalist()
  {
    DBG_ASSERT(!collecting_output_);
    journals_.clear();
  }

  void Journalist::Printf( EJournalLevel level, EJournalCategory category,
                           const char* pformat, ... ) const
  {
    // wrap the arguments and pass to VPrintf
    va_list ap;
    va_start(ap, pformat);

    VPrintf(level, category, pformat, ap);

    va_end(ap);
  }

  void Journalist::PrintStringOverLines(EJournalLevel level,
                                        EJournalCategory category,
                                        Index indent_spaces, Index max_length,
                                        const std::string& line) const
  {
    DBG_ASSERT(indent_spaces + max_length + 1 < 1024);
    char buffer[1024];
    std::string::size_type last_line_pos = 0;
    std::string::size_type last_word_pos = 0;
    bool first_line = true;
    Index buffer_pos = 0;

    while (last_line_pos < line.length()) {
      std::string::size_type line_pos = last_line_pos;
      Index curr_length = 0;
      while (curr_length < max_length && line_pos < line.length()) {
        buffer[buffer_pos] = line[line_pos];
        if (line[line_pos] == ' ') {
          last_word_pos = line_pos+1;
        }
        curr_length++;
        buffer_pos++;
        line_pos++;
      }
      if (line_pos == line.length()) {
        // This is the last line to be printed.
        buffer[buffer_pos] = '\0';
        Printf(level, category, "%s", buffer);
        break;
      }
      if (last_word_pos == last_line_pos) {
        if (line[line_pos]==' ') {
          buffer[buffer_pos] = '\0';
          last_word_pos = line_pos+1;
          last_line_pos = line_pos+1;
        }
        else {
          // The current word is too long to fit into one line
          // split word over two lines
          buffer[buffer_pos-1] = '-';
          buffer[buffer_pos] = '\0';
          last_word_pos = line_pos-1;
          last_line_pos = last_word_pos;
        }
      }
      else {
        // insert '\0' character after last complete word
        buffer[buffer_pos-(line_pos-last_word_pos)-1] = '\0';
        last_line_pos = last_word_pos;
      }

      Printf(level, category, "%s\n", buffer);
      if (first_line) {
        for (Index i=0; i<indent_spaces; i++) {
          buffer[i] = ' ';
        }
        first_line = false;
      }
      buffer_pos = indent_spaces;
    }
  }

  void Journalist::PrintfIndented( EJournalLevel level,
                                   EJournalCategory category, Index indent_level,
                                   const char* pformat, ... ) const
  {
    // wrap the arguments and pass to VPrintfIndented
    va_list ap;
    va_start(ap, pformat);

    VPrintfIndented(level, category, indent_level, pformat, ap);

    va_end(ap);
  }

  //   void Journalist::PrintVector(EJournalLevel level,
  //                                EJournalCategory category,
  //                                const std::string& name,
  //                                const Vector& vector,
  //                                Index indent,
  //                                const std::string prefix) const
  //   {
  //     // print the msg on every journal that accepts
  //     // the category and output level
  //     for (Index i=0; i<(Index)journals_.size(); i++) {
  //       if (journals_[i]->IsAccepted(category, level)) {
  //         // print the message
  //         journals_[i]->PrintVector(name, vector, indent, prefix);
  //       }
  //     }
  //   }

  //   void Journalist::PrintMatrix(EJournalLevel level,
  //                                EJournalCategory category,
  //                                const std::string& name,
  //                                const Matrix& matrix,
  //                                Index indent /*=0*/,
  //                                std::string prefix /*=""*/) const
  //   {
  //     // print the msg on every journal that accepts
  //     // the category and output level
  //     for (Index i=0; i<(Index)journals_.size(); i++) {
  //       if (journals_[i]->IsAccepted(category, level)) {
  //         // print the message
  //         journals_[i]->PrintMatrix(name, matrix, indent, prefix);
  //       }
  //     }
  //   }

  void Journalist::VPrintf(
    EJournalLevel level,
    EJournalCategory category,
    const char* pformat, va_list ap) const
  {
    // print the msg on every journal that accepts
    // the category and output level
    for (Index i=0; i<(Index)journals_.size(); i++) {
      if (journals_[i]->IsAccepted(category, level)) {
        if (collecting_output_) {
#ifdef HAVE_VA_COPY
          va_list apcopy;
          va_copy(apcopy, ap);
          vsnprintf(buffer_, 32767, pformat, apcopy);
          va_end(apcopy);
#else

          vsnprintf(buffer_, 32767, pformat, ap);
#endif

          journals_[i]->CollectionBuffer() += buffer_;
        }
        else {
          if (my_rank_ != 0) return;
          // print the message
#ifdef HAVE_VA_COPY
          va_list apcopy;
          va_copy(apcopy, ap);
          journals_[i]->Printf(category, level, pformat, apcopy);
          va_end(apcopy);
#else

          journals_[i]->Printf(category, level, pformat, ap);
#endif
        }
      }
    }
  }

  void Journalist::VPrintfIndented(
    EJournalLevel level,
    EJournalCategory category,
    Index indent_level,
    const char* pformat, va_list ap) const
  {
    // print the msg on every journal that accepts
    // the category and output level
    for (Index i=0; i<(Index)journals_.size(); i++) {
      if (journals_[i]->IsAccepted(category, level)) {

        if (collecting_output_) {
          for (Index s=0; s<indent_level; s++) {
            journals_[i]->CollectionBuffer() += "  ";
          }
#ifdef HAVE_VA_COPY
          va_list apcopy;
          va_copy(apcopy, ap);
          vsnprintf(buffer_, 32767, pformat, apcopy);
          va_end(apcopy);
#else

          vsnprintf(buffer_, 32767, pformat, ap);
#endif

          journals_[i]->CollectionBuffer() += buffer_;
        }
        else {
          if (my_rank_ != 0) return;
          // indent the appropriate amount
          for (Index s=0; s<indent_level; s++) {
            journals_[i]->Print(category, level, "  ");
          }

          // print the message
#ifdef HAVE_VA_COPY
          va_list apcopy;
          va_copy(apcopy, ap);
          journals_[i]->Printf(category, level, pformat, apcopy);
          va_end(apcopy);
#else

          journals_[i]->Printf(category, level, pformat, ap);
#endif

        }
      }
    }
  }

  bool Journalist::ProduceOutput(EJournalLevel level,
                                 EJournalCategory category) const
  {
    for (Index i=0; i<(Index)journals_.size(); i++) {
      if (journals_[i]->IsAccepted(category, level)) {
        return true;
      }
    }
    return false;
  }

  bool Journalist::AddJournal(const SmartPtr<Journal> jrnl)
  {
    DBG_ASSERT(!collecting_output_);
    DBG_ASSERT(IsValid(jrnl));
    std::string name = jrnl->Name();

    SmartPtr<Journal> temp = GetJournal(name);
    DBG_ASSERT(IsNull(temp));
    if (IsValid(temp)) {
      return false;
    }

    journals_.push_back(jrnl);
    return true;
  }

  SmartPtr<Journal> Journalist::AddFileJournal(
    const std::string& journal_name,
    const std::string& fname,
    EJournalLevel default_level
  )
  {
    SmartPtr<FileJournal> temp = new FileJournal(journal_name, default_level);

    // Open the file (Note:, a fname of "stdout" is handled by the
    // Journal class to mean stdout, etc.
    // We don't want to open the journal if we are not the root node
    if ( (my_rank_!=0 || temp->Open(fname.c_str())) &&
         AddJournal(GetRawPtr(temp))) {
      return GetRawPtr(temp);
    }
    return NULL;
  }

  void Journalist::FlushBuffer() const
  {
    DBG_ASSERT(!collecting_output_);
    if (my_rank_ != 0) return;
    for (Index i=0; i<(Index)journals_.size(); i++) {
      journals_[i]->FlushBuffer();
    }
  }

  SmartPtr<Journal> Journalist::GetJournal(
    const std::string& journal_name
  )
  {
    SmartPtr<Journal> retValue = NULL;

    // try to find the journal
    for (Index i=0; i<(Index)journals_.size(); i++) {
      SmartPtr<Journal> tmp = journals_[i];
      if (tmp->Name() == journal_name) {
        retValue = tmp;
        break;
      }
    }

    return retValue;
  }

  void Journalist::StartDistributedOutput() const
  {
    DBG_ASSERT(!collecting_output_);
#ifdef HAVE_MPI
    collecting_output_ = true;
#endif
  }

  void Journalist::FinishDistributedOutput() const
  {
#ifdef HAVE_MPI
    DBG_ASSERT(collecting_output_);
    collecting_output_ = false;

    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int* recvcounts = NULL;
    int* displs = NULL;
    char* total_buffer = NULL;
    if (my_rank_==0) {
      recvcounts = new int[num_procs];
      displs = new int[num_procs];
    }

    for (Index i=0; i<(Index)journals_.size(); i++) {
      int local_buflen = journals_[i]->CollectionBuffer().length()+1;

      // Send buffer sizes to root process for memory allocation
      MPI_Gather(&local_buflen, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (my_rank_==0) {
        displs[0] = 0;
        for (int p=1; p<num_procs; p++) {
          displs[p] = displs[p-1] + recvcounts[p-1];
        }
        const int total_buflen = displs[num_procs-1] + recvcounts[num_procs-1];
        total_buffer = new char[total_buflen];
      }

      // Send all buffer to the root process
      const char* local_buffer = journals_[i]->CollectionBuffer().c_str();
      MPI_Gatherv(const_cast<char*>(local_buffer), local_buflen, MPI_CHAR, total_buffer, recvcounts, displs, MPI_CHAR, 0, MPI_COMM_WORLD);

      // On root process, do the output
      if (my_rank_==0) {
        for (int p=0; p<num_procs; p++) {
          journals_[i]->Print(J_ANY, J_NONE, &total_buffer[displs[p]]);
        }
        delete [] total_buffer;
      }

      journals_[i]->CollectionBuffer().clear();
    }

    delete [] recvcounts;
    delete [] displs;
#endif
  }

  void Journalist::DeleteAllJournals()
  {
    for (Index i=0; i<(Index)journals_.size(); i++) {
      journals_[i]=NULL;
    }
    journals_.resize(0);
  }

  ///////////////////////////////////////////////////////////////////////////
  //                 Implementation of the Journal class                   //
  ///////////////////////////////////////////////////////////////////////////

  Journal::Journal(
    const std::string& name,
    EJournalLevel default_level
  )
      :
      name_(name)
  {
    for (Index i=0; i<J_LAST_CATEGORY; i++) {
      print_levels_[i] = default_level;
    }
  }

  Journal::~Journal()
  {}

  std::string Journal::Name()
  {
    return name_;
  }

  bool Journal::IsAccepted(
    EJournalCategory category,
    EJournalLevel level
  ) const
  {
    if (category == J_ANY) return true;
    if (print_levels_[(Index)category] >= (Index) level) {
      return true;
    }

    return false;
  }

  void Journal::SetPrintLevel(
    EJournalCategory category,
    EJournalLevel level)
  {
    DBG_ASSERT((Index)category >= 0);
    print_levels_[(Index)category] = (Index) level;
  }

  void Journal::SetAllPrintLevels(
    EJournalLevel level)
  {
    for (Index category=(Index)J_DBG;
         category<(Index)J_USER_APPLICATION;
         category++) {
      print_levels_[category] = (Index) level;
    }
  }


  ///////////////////////////////////////////////////////////////////////////
  //                 Implementation of the FileJournal class               //
  ///////////////////////////////////////////////////////////////////////////

  FileJournal::FileJournal(
    const std::string& name,
    EJournalLevel default_level
  )
      :
      Journal(name, default_level),
      file_(NULL)
  {}

  FileJournal::~FileJournal()
  {
    if (file_ && file_ != stdout && file_ != stderr) {
      // close the file
      fclose(file_);
    }
    file_ = NULL;
  }


  bool FileJournal::Open(const char* fname)
  {
    if (file_ && file_ != stdout && file_ != stderr) {
      // file already opened, close it
      fclose(file_);
    }
    file_ = NULL;

    if (strcmp("stdout", fname)==0) {
      file_=stdout;
      return true;
    }
    else if (strcmp("stderr", fname)==0) {
      file_=stderr;
      return true;
    }
    else {
      // open the file on disk
      file_ = fopen(fname, "w+");
      if (file_) {
        return true;
      }
    }

    return false;
  }


  void FileJournal::PrintImpl(EJournalCategory category, EJournalLevel level,
                              const char* str)
  {
    DBG_START_METH("Journal::Print", 0);
    if (file_) {
      fprintf(file_, str);
      DBG_EXEC(0, fflush(file_));
    }
  }

  void FileJournal::PrintfImpl(EJournalCategory category, EJournalLevel level,
                               const char* pformat, va_list ap)
  {
    DBG_START_METH("Journal::Printf", 0);
    if (file_) {
      vfprintf(file_, pformat, ap);
      DBG_EXEC(0, fflush(file_));
    }
  }

  void FileJournal::FlushBufferImpl()
  {
    if (file_) {
      fflush(file_);
    }
  }


  ///////////////////////////////////////////////////////////////////////////
  //                 Implementation of the StreamJournal class               //
  ///////////////////////////////////////////////////////////////////////////

  StreamJournal::StreamJournal(
    const std::string& name,
    EJournalLevel default_level
  )
      :
      Journal(name, default_level),
      os_(NULL)
  {}

  void StreamJournal::SetOutputStream(std::ostream* os)
  {
    os_ = os;
  }

  void StreamJournal::PrintImpl(EJournalCategory category, EJournalLevel level,
                                const char* str)
  {
    DBG_START_METH("StreamJournal::PrintImpl", 0);
    if (os_) {
      *os_ << str;
      DBG_EXEC(0, *os_ << std::flush);
    }
  }

  void StreamJournal::PrintfImpl(EJournalCategory category, EJournalLevel level,
                                 const char* pformat, va_list ap)
  {
    DBG_START_METH("StreamJournal::PrintfImpl", 0);
    if (os_) {
      vsnprintf(buffer_, 32767, pformat, ap);
      *os_ << buffer_;
      DBG_EXEC(0, *os_ << std::flush);
    }
  }

  void StreamJournal::FlushBufferImpl()
  {
    if (os_) {
      *os_ << std::flush;
    }
  }

} // namespace Ipopt
