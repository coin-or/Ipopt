// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter, Sanjeeb Dash     IBM    2009-06-17
//             (based on IpTripletHelper.hpp rev 1312)

#ifndef __IPPARTRIPLETHELPER_HPP__
#define __IPPARTRIPLETHELPER_HPP__

#include "IpTypes.hpp"
#include "IpException.hpp"

namespace Ipopt
{

  /** forward declarations */
  class Matrix;
  class ParGenMatrix;
  class ParSymMatrix;
  class DiagMatrix;
  class IdentityMatrix;
  class ExpansionMatrix;
  class ParExpansionMatrix;
  class ScaledMatrix;
  class SymScaledMatrix;
  class SumMatrix;
  class SumSymMatrix;
  class ZeroMatrix;
  class CompoundMatrix;
  class CompoundSymMatrix;
  class TransposeMatrix;
  class Vector;
  class CompoundVector;
  class ParVector;

  class ParTripletHelper
  {
  public:
    /**@name A set of recursive routines that help with the Triplet format. */
    //@{
    /** find the total number of triplet entries of a Matrix */
    static Index GetNumberEntries(const Matrix& matrix);

    /** fill the irows, jcols structure for the triplet format from the matrix */
    static void FillRowCol(Index n_entries, const Matrix& matrix, Index* iRow, Index* jCol, Index row_offset=0, Index col_offset=0);

    /** fill the values for the triplet format from the matrix */
    static void FillValues(Index n_entries, const Matrix& matrix, Number* values);

    /** fill the values from a possibly distributed vector into one
     *  single array */
    static void FillAllValuesFromVector(Index n_entries, const Vector& vector, Number* values);

    /** fill the values into a possibly distributed vector from one
     *  single array */
    static void PutAllValuesInVector(Index n_entries, const Number* values, Vector& vector);

    /** find the number of local entries induced by Vector for the
     *  elements residing in this process. */
    static Index GetLocalNumberEntries(const Vector& vector);

    /** Fill array with global positions of local entries.  global_pos
     *  must have length as least GetLocalNumberEntries. */
    static void GetGlobalPos(Index n_local_entries, const Vector& vector,
			     Index* global_pos, Index offset = 0);

    /** fill local_values with the local elements of vector */
    static void FillLocalValuesFromVector(Index n_local_entries,
					  const Vector& vector,
					  Number* local_values);

    /** fill the local elements of vector with the elements in local_values */
    static void PutLocalValuesInVector(Index n_local_entries,
				       const Number* local_values,
				       Vector& vector);

    //@}

  private:
    /** find the total number of triplet entries for the SumMatrix */
    static Index GetNumberEntries_(const SumMatrix& matrix);

    /** find the total number of triplet entries for the SumSymMatrix */
    static Index GetNumberEntries_(const SumSymMatrix& matrix);

    /** find the total number of triplet entries for the CompoundMatrix */
    static Index GetNumberEntries_(const CompoundMatrix& matrix);

    /** find the total number of triplet entries for the CompoundSymMatrix */
    static Index GetNumberEntries_(const CompoundSymMatrix& matrix);

    /** find the total number of triplet entries induced by Vector */
    static Index GetNumberEntries_(const Vector& vector);

    static void FillRowCol_(Index n_entries, const ParGenMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const ParGenMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const ParSymMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const ParSymMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const DiagMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const DiagMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const IdentityMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const IdentityMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const ParExpansionMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const ParExpansionMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const ExpansionMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const ExpansionMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const Vector& vector, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const Vector& vector, Number* values);

    static void FillRowCol_(Index n_entries, const SumMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const SumMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const SumSymMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const SumSymMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const CompoundMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const CompoundMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const CompoundSymMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const CompoundSymMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const ScaledMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const ScaledMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const SymScaledMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const SymScaledMatrix& matrix, Number* values);

    static void FillRowCol_(Index n_entries, const TransposeMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol);

    static void FillValues_(Index n_entries, const TransposeMatrix& matrix, Number* values);

    static void FillAllValuesFromVector_(Index n_entries, const CompoundVector& vector, Number* values);

    static void FillAllValuesFromVector_(Index n_entries, const ParVector& vector, Number* values);

    static void PutAllValuesInVector_(Index n_entries, const Number* values, CompoundVector& vector);

    static void PutAllValuesInVector_(Index n_entries, const Number* values, ParVector& vector);

    static Index GetLocalNumberEntries_(const CompoundVector& vector);

    static Index GetLocalNumberEntries_(const ParVector& vector);

    static void GetGlobalPos_(Index n_local_entries, const CompoundVector& vector,
			      Index* global_pos, Index offset);

    static void GetGlobalPos_(Index n_local_entries, const ParVector& vector,
			      Index* global_pos, Index offset);

    static void FillLocalValuesFromVector_(Index n_local_entries,
					  const CompoundVector& vector,
					  Number* local_values);

    static void FillLocalValuesFromVector_(Index n_local_entries,
					  const ParVector& vector,
					  Number* local_values);

    static void PutLocalValuesInVector_(Index n_local_entries,
					const Number* local_values,
					CompoundVector& vector);

    static void PutLocalValuesInVector_(Index n_local_entries,
					const Number* local_values,
					ParVector& vector);
  };
} // namespace Ipopt

#endif
