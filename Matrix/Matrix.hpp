/*

(C) Matthew Swanson

This file is part of _PROJECT_NAME_.

_PROJECT_NAME_ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

_PROJECT_NAME_ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with _PROJECT_NAME_.  If not, see <http://www.gnu.org/licenses/>.

*/
#pragma once

#include <cstddef>
#include <cstring>
#include <cmath>
#include <utility>
#include <string>
#include <sstream>
#include <array>
#include <memory>
#include <iomanip>


namespace cg {

const float PI = 3.14159265359f;
const float DegreeToRadian = (PI / 180);
const float RadianToDegree = (180 / PI);


namespace rmat {


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////RAW MATRIX//////////
///////////////////////////////////////////////////////////////////////////////////////////////////


template<std::size_t Columns, typename DataType>
using RowMatrix = std::array<DataType, Columns>;

template<std::size_t Rows, typename DataType>
using ColumnMatrix = std::array<RowMatrix<1, DataType>, Rows>;

template<std::size_t Rows, std::size_t Columns, typename DataType>
using RawMatrix = std::array<RowMatrix<Columns, DataType>, Rows>;

template<std::size_t Rows, typename DataType>
using RawVec = ColumnMatrix<Rows, DataType>;

template<typename DataType>
using Vector2 = RawVec<2, DataType>;

template<typename DataType>
using Vector3 = RawVec<3, DataType>;

template<typename DataType>
using Vector4 = RawVec<4, DataType>;

//////////////////////////////////////////////////////////////////////////////////RAW MATRIX///////
///////////////////////////////////////////////////////////////////////////////Conversions/////////
///////////////////////////////////////////////////////////////////////////////////////////////////


/**View a line array as a RawMatrix.
\param m The line array.
\return A RawMatrix reference to the line array.*/
template<std::size_t _R, std::size_t _C, typename _T>
RawMatrix<_R, _C, _T>*
ViewAsMatrix(std::array<_T, _R*_C>& m)
{
	return ((RawMatrix<_R, _C, _T>*)&m);
}
/**View a line array as a RawMatrix.
\param m The line array.
\return A RawMatrix reference to the line array.*/
template<std::size_t _R, std::size_t _C, typename _T>
const RawMatrix<_R, _C, _T>*
ViewAsMatrix(const std::array<_T, _R*_C>& m)
{
	return ((const RawMatrix<_R, _C, _T>*)&m);
}

/**Get a pointer to a matrix as a single line matrix. Example: A 3x3 matrix
will be viewed as a 1x9 matrix. The 1x9 matrix will reflect changes made to
the 3x3 matrix and vice versa.
\param m The RawMatrix<R,C,T>
\return A reference to the RawMatrix as a std::array<T,R*C>.*/
template<std::size_t Rows, std::size_t Columns, typename DataType>
std::array<DataType, Rows*Columns>*
ViewAsArray(RawMatrix<Rows, Columns, DataType>& m)
{
	return ((std::array<DataType, Rows*Columns>*)&m);
}
/**Get a pointer to a matrix as a single line matrix. Example: A 3x3 matrix
will be viewed as a 1x9 matrix. The 1x9 matrix will reflect changes made to
the 3x3 matrix and vice versa.
\param m The RawMatrix<R,C,T>
\return A reference to the RawMatrix as a std::array<T,R*C>.*/
template<std::size_t Rows, std::size_t Columns, typename DataType>
const std::array<DataType, Rows*Columns>*
ViewAsArray(const RawMatrix<Rows, Columns, DataType>& m)
{
	return ((const std::array<DataType, Rows*Columns>*)&m);
}
//////////////////////////////////////////////////////////////////////////////////RAW MATRIX///////
///////////////////////////////////////////////////////////////////////////////Creations///////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/**Create a raw matrix.
\param nums A std::array of numbers, which will be copied. EAch _R numbers will
be a row.
\return The raw matrix.*/
template<std::size_t _R, std::size_t _C, typename _T>
RawMatrix<_R, _C, _T> Create(std::array<_T, _R*_C>&& nums)
{
	RawMatrix<_R, _C, _T> ret;
	auto& maskRet = *rmat::ViewAsArray(ret);
	maskRet = std::move(nums);
	return ret;
}
/**Create a raw matrix.
\param nums A std::array of numbers, which will be copied. EAch _R numbers will
be a row.
\return The raw matrix.*/
template<std::size_t _R, std::size_t _C, typename _T>
RawMatrix<_R, _C, _T> Create(const std::array<_T, _R*_C>& nums)
{
	RawMatrix<_R, _C, _T> ret;
	std::memcpy(ret.data(), nums.data(), _R*_C * sizeof(_T));
	return ret;
}


//////////////////////////////////////////////////////////////////////////////////RAW MATRIX///////
///////////////////////////////////////////////////////////////////////////////Accessors///////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/**Print a matrix to a string.
\param m The matrix.
\return The matrix as a string.*/
template<std::size_t _R, std::size_t _C, typename _T>
std::string ToString(const RawMatrix<_R, _C, _T>& m)
{
	std::size_t digitSize = 1;
	for (std::size_t i = 0; i < _R; ++i)
		for (std::size_t j = 0; j < _C; ++j)
		{
			auto s = std::to_string(m[i][j]).size();
			if (digitSize < s)
				digitSize = s;
		}
	std::stringstream ss;
	ss << std::string((_C - 1) * (3 + digitSize) + (2 + digitSize), '-')
		<< std::endl << "|";
	for (std::size_t i = 0; i < _R; ++i)
		for (std::size_t j = 0; j < _C; ++j)
		{
			if (j == 0 && i != 0)
				ss << "|" << std::endl << "|";
			if (j != 0)
				ss << " , ";
			ss << std::setfill(' ') << std::setw(digitSize);
			ss << std::to_string(m[i][j]);
		}
	ss << "|" << std::endl
		<< std::string((_C - 1) * (3 + digitSize) + (2 + digitSize), '-');
	return ss.str();
}

/**Extract a column from a matrix.
\param c The column to get.
\param m The matrix for which to extract the column.
\return A ColumnMatrix that is the column (not a reference though).*/
template<std::size_t _R, std::size_t _C, typename _T>
ColumnMatrix<_R, _T> ColumnOf(std::size_t c, const RawMatrix<_R, _C, _T>& m)
{
	auto ret = Create<_R, 1, _T>({});
	for (std::size_t i = 0; i < _R; ++i)
		ret[i][0] = m[i][c];
	return ret;
}
/**Extract a row from a matrix.
\param r The row to get.
\param m The matrix for which to extract the row.
\return A RowMatrix that is the row (not a reference though).*/
template<std::size_t _R, std::size_t _C, typename _T>
RowMatrix<_C, _T> RowOf(std::size_t r, const RawMatrix<_R, _C, _T>& m)
{
	return m[r];
}

/**Set a column in a matrix.
\param c The ColumnMatrix to set.
\param i The index of `m` to set.
\param m The matrix to alter.*/
template<std::size_t _R, std::size_t _C, typename _T>
void SetColumn(std::size_t i, const ColumnMatrix<_R, _T>& c,
	RawMatrix<_R, _C, _T>& m)
{
	for (std::size_t a = 0; a < _R; ++a)
		m[a][i] = c[a][0];
}
/**Set a row in a matrix.
\param r The RowMatrix to set.
\param i The index of `m` to set.
\param m The matrix to alter.*/
template<std::size_t _R, std::size_t _C, typename _T>
void SetRow(const RowMatrix<_R, _T>& r, RawMatrix<_R, _C, _T>& m, std::size_t i)
{
	m[i] = r;
}

//////////////////////////////////////////////////////////////////////////////////RAW MATRIX///////
///////////////////////////////////////////////////////////////////////////////Operators///////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/**Get the dot product of a single row and column.
\param row The row in the op.
\param colp The column in the op.
\return The single number that is the answer to the DOT operation.*/
template<std::size_t S, typename T>
T Mul(const RowMatrix<S, T>& row, const ColumnMatrix<S, T>& colp)
{
	auto& col = *rmat::ViewAsArray<S, 1, T>(colp);
	T result = 0;
	for (std::size_t i = 0; i < S; ++i)
	{
		result += row[i] * col[i];
	}
	return result;
}

/**Operate on raw matrices.
\param lhs The first matrix.
\param rhs The second matrix.
\return A matrix that is the result of the operation.*/
template<std::size_t _Rw, std::size_t _Cl, std::size_t _Cl2, typename _T>
RawMatrix<_Rw, _Cl2, _T> Mul(const RawMatrix<_Rw, _Cl, _T>& lhs,
	const RawMatrix<_Cl, _Cl2, _T>& rhs)
{
	auto ret = Create<_Rw, _Cl2, _T>({});
	for (std::size_t i = 0; i < _Rw; ++i)
		for (std::size_t j = 0; j < _Cl2; ++j)
			ret[i][j] = Mul(lhs[i], ColumnOf(j, rhs));
	return ret;
}

/**Scale a matrix by a scalar.
\param s The scalar.
\param m The matrix. It will be modified.*/
template<std::size_t _R, std::size_t _C, typename _T>
void Scale(RawMatrix<_R, _C, _T>& m, const _T& s)
{
	auto& retV = *rmat::ViewAsArray(m);
	for (_T& t : retV)
		t *= s;
}
/**Scale a matrix by a scalar.
\param s The scalar.
\param m The matrix. It will be modified.*/
template<std::size_t _C, typename _T>
void Scale(RowMatrix<_C, _T>& m, const _T& s)
{
	for (_T& t : m)
		t *= s;
}
/**Scale a matrix by a scalar.
\param s The scalar.
\param m The matrix. It will be modified.*/
template<std::size_t _R, typename _T>
void Scale(ColumnMatrix<_R, _T>& m, const _T& s)
{
	auto& ret = *rmat::ViewAsArray(m);
	for (_T& t : ret)
		t *= s;
}
/**Scale a matrix by a scalar.
\param s The scalar.
\param m The matrix.
\return A copy of `m` that was scaled.*/
template<std::size_t _R, std::size_t _C, typename _T>
RawMatrix<_R, _C, _T> Mul(const RawMatrix<_R, _C, _T>& m, const _T& s)
{
	auto copy = m;
	Scale(copy, s);
	return copy;
}

/**Add two matrices together. Must be the same size.
\param m1 The first matrix.
\param m2 The second matrix.
\return The result matrix.*/
template<std::size_t _R, std::size_t _C, typename _T>
RawMatrix<_R, _C, _T> Add(const RawMatrix<_R, _C, _T>& m1,
	const RawMatrix<_R, _C, _T>& m2)
{
	auto ret = m1;
	auto& l1 = *rmat::ViewAsArray(ret);
	auto& l2 = *rmat::ViewAsArray(m2);
	for (std::size_t i = 0; i < _R*_C; ++i)
		l1[i] += l2[i];
	return ret;
}
/**Add two matrices together. Must be the same size.
\param m1 The first matrix.
\param m2 The second matrix.
\return The result matrix.*/
template<std::size_t _C, typename _T>
RowMatrix<_C, _T> Add(const RowMatrix<_C, _T>& m1,
	const RowMatrix<_C, _T>& m2)
{
	auto ret = m1;
	for (std::size_t i = 0; i < _C; ++i)
		ret[i] += m2[i];
	return ret;
}
/**Add two matrices together. Must be the same size.
\param m1 The first matrix.
\param m2 The second matrix.
\return The result matrix.*/
template<std::size_t _C, typename _T>
RowMatrix<_C, _T> Sub(const RowMatrix<_C, _T>& m1,
	const RowMatrix<_C, _T>& m2)
{
	auto ret = m1;
	for (std::size_t i = 0; i < _C; ++i)
		ret[i] -= m2[i];
	return ret;
}
/**Add two matrices together. Must be the same size.
\param m1 The first matrix.
\param m2 The second matrix.
\return The result matrix.*/
template<std::size_t _R, typename _T>
ColumnMatrix<_R, _T> Add(const ColumnMatrix<_R, _T>& m1,
	const ColumnMatrix<_R, _T>& m2)
{
	auto ret = m1;
	auto& retA = *rmat::ViewAsArray(ret);
	auto& a2 = *rmat::ViewAsArray(m2);
	for (std::size_t i = 0; i < _R; ++i)
		retA[i] += a2[i];
	return ret;
}
/**Subtraction of matrix. Must be the same size. (m1 - m2)
\param m1 The first matrix.
\param m2 The second matrix.
\return The result matrix.*/
template<std::size_t _R, std::size_t _C, typename _T>
RawMatrix<_R, _C, _T> Sub(const RawMatrix<_R, _C, _T>& m1,
	const RawMatrix<_R, _C, _T>& m2)
{
	auto ret = m1;
	auto& l1 = *rmat::ViewAsArray(ret);
	auto& l2 = *rmat::ViewAsArray(m2);
	for (std::size_t i = 0; i < _R*_C; ++i)
		l1[i] -= l2[i];
	return ret;
}
//////////////////////////////////////////////////////////////////////////////////RAW MATRIX///////
///////////////////////////////////////////////////////////////////////////////Row/Col Ops/////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/**Do a complex row operation. [s1*r1 + s2*r2 -> r2]
\param r1 The number of the first row.
\param s1 The scalar for the first row.
\param r2 The second row number
\param s2 The scalar for the second row.
\param m The RawMatrix.  Will be modified.*/
template<std::size_t _R, std::size_t _C, typename _T>
void RowOp(RawMatrix<_R, _C, _T>& m, std::size_t r1, const _T& s1,
	std::size_t r2, const _T& s2)
{
	auto row1 = *((RowMatrix<_C, _T>*)&m[r1]);
	auto& row2 = *((RowMatrix<_C, _T>*)&m[r2]);
	Scale(row1, s1);
	Scale(row2, s2);
	row2 = Add(row1, row2);
}
/**Do a complex column operation. [s1*c1 + s2*c2 -> c2]
\param c1 The number of the first col.
\param s1 The scalar for the first col.
\param c2 The second col number
\param s2 The scalar for the second col.
\param m The RawMatrix.  Will be modified.*/
template<std::size_t _R, std::size_t _C, typename _T>
void ColumnOp(RawMatrix<_R, _C, _T>& m, std::size_t c1, const _T& s1,
	std::size_t c2, const _T& s2)
{
	auto col1 = ColumnOf(c1, m);
	auto col2 = ColumnOf(c2, m);
	Scale(col1, s1);
	Scale(col2, s2);
	SetColumn(c2, Add(col1, col2), m);
}

/**Swap tow rows.
\param r1 The first row.
\param r2 The second row.
\param m The matrix for which to swap rows.*/
template<std::size_t _R, std::size_t _C, typename _T>
void RowSwap(RawMatrix<_R, _C, _T>& m, std::size_t r1, std::size_t r2)
{
	RowMatrix<_C, _T> t = m[r1];
	m[r1] = m[r2];
	m[r2] = t;
}

/**Swap two columns.
\param c1 The first column.
\param c2 The second column.
\param m The matrix for which to swap columns.*/
template<std::size_t _R, std::size_t _C, typename _T>
void ColumnSwap(RawMatrix<_R, _C, _T>& m, std::size_t c1, std::size_t c2)
{
	auto col1 = ColumnOf(c1, m);
	SetColumn(c1, ColumnOf(c2, m), m);
	SetColumn(c2, col1, m);
}

/**Transpose a matrix.
\param m The matrix to transpose.
\return The transposed matrix.*/
template<std::size_t _R, std::size_t _C, typename _T>
RawMatrix<_C, _R, _T> Transpose(const RawMatrix<_R, _C, _T>& m)
{
	RawMatrix<_C, _R, _T> ret;
	for (std::size_t i = 0; i < _R; ++i)
		for (std::size_t j = 0; j < _C; ++j)
			ret[i][j] = m[j][i];
	return ret;
}

/**Get a zero matrix.
\return A matrix will all zero.*/
template<std::size_t _R, std::size_t _C, typename _T>
RawMatrix<_C, _R, _T> ZeroMatrix()
{
	RawMatrix<_C, _R, _T> ret = {};
	return ret;
}

}
///////////////////////////////////////////////////////////////////////////////END RMAT NS/////////
///////////////////////////////////////////////////////////////////////////////END RAW MATRIX//////
///////////////////////////////////////////////////////////////////////////////END RAW MATRIX//////


///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////MATRIX/////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/**Matrix class.*/
template<std::size_t Rows, std::size_t Columns, typename DataType>
class Matrix
{
public:
	/**The data type to use.*/
	using T = DataType;
	/**The amount of rows.*/
	const static std::size_t _R = Rows;
	/**The amount of columnes.*/
	const static std::size_t _C = Columns;
	/**The Self type.*/
	using SelfT = Matrix<_R, _C, T>;
	/**The transpose type.*/
	using TransSelfT = Matrix<_C, _R, T>;
	/************************************************************Constructors*/

	/***************************************************************Utilities*/

	/**Get a pointer to the data as type T* .
	\return A pointer to the data.*/
	T* Begin()
	{
		return (T*)&m_data;
	}
	/**Get a pointer to the data as type T* .
	\return A pointer to the data.*/
	const T* Begin()const
	{
		return (const T*)&m_data;
	}
	/**Get one past the end pointer to data.
	\param A one-past-the-end pointer.*/
	T* End()
	{
		return ((T*)&m_data) + _R*_C;
	}
	/**Get one past the end pointer to data.
	\param A one-past-the-end pointer.*/
	const T* End()const
	{
		return ((const T*)&m_data) + _R*_C;
	}

	/**Augment (attach) a column to the matrix.
	\param c The colum to attach (the the last position).
	\return The new matrix.*/
	Matrix<_R, _C + 1, T> AugmentColumn(const Matrix<_R, 1, T>& c) const
	{
		Matrix<_R, _C + 1, T> ret;
		for (std::size_t i = 0; i < _C; ++i)
			ret.SetColumn(i, Column(i));
		ret.SetColumn(_C, c);
		return ret;
	}
	/**Augment (attach) a row to the matrix.
	\param r The row to attach (the the last position).
	\return The new matrix.*/
	Matrix<_R + 1, _C, T> AugmentRow(const Matrix<1, _C, T>& r) const
	{
		Matrix<_R + 1, _C, T> ret;
		for (std::size_t i = 0; i < _R; ++i)
			ret.SetRow(i, Row(i));
		ret.SetRow(_R, r);
		return ret;
	}

	/**Flip the matrix horizontally.*/
	void FlipHV()
	{
		FlipH();
		FlipV();
	}
	/**Flip the matrix horizontally.*/
	void FlipH()
	{
		const static std::size_t amt = _C / 2;
		for (std::size_t i = 0; i < amt; ++i)
			SwapColumn(i, (_C - 1) - i);
	}
	/**Flip the matrix virtically.*/
	void FlipV()
	{
		const static std::size_t amt = _R / 2;
		for (std::size_t i = 0; i < amt; ++i)
			SwapRow(i, (_R - 1) - i);
	}

	/**Get a std::array with the copied data. It will be a single array of size
	`Rows`*`Columns` of type `DataType`.
	\return An std::array with all the data copied to it.*/
	std::array<T, _R*_C> RawCopy() const
	{
		return *rmat::ViewAsArray(m_data);
	}

	/**Get an column as an array.
	\param c The column number to get.
	\return A Matrix<Rows, 1,T> object with the column.*/
	Matrix<_R, 1, T> Column(std::size_t c) const
	{
		Matrix<_R, 1, T> col;
		for (std::size_t i = 0; i < _R; ++i)
			col.Set(i, 0, Get(i, c));
		return col;
	}

	/**Get a row.
	\param r The row to get.
	\return A Matrix<1,Columns,T> with the row.*/
	Matrix<1, _C, T> Row(std::size_t r) const
	{
		Matrix<1, _C, T> row;
		for (std::size_t i = 0; i < _C; ++i)
			row.Set(0, i, Get(r, i));
		return row;
	}

	/**Set a column to a Matrix<Rows, 1,T>
	\param c The column to set.
	\param col The single row matrix.*/
	void SetColumn(std::size_t c, const Matrix<Rows, 1, T>& col)
	{
		for (std::size_t i = 0; i < _R; ++i)
			Set(i, c, col.Get(i, 0));
	}
	/**Set a row to a Matrix<1, Columns,T>
	\param r The row to set.
	\param row The single col matrix.*/
	void SetRow(std::size_t r, const Matrix<1, Columns, T>& row)
	{
		for (std::size_t i = 0; i < _C; ++i)
			Set(r, i, row.Get(0, i));
	}

	/**Swap two rows.
	\param r1 The first row in the swap.
	\param r2 The second row in the swap.*/
	void SwapRow(std::size_t r1, std::size_t r2)
	{
		rmat::RowSwap(m_data, r1, r2);
	}

	/**Swap two columns.
	\param c1 The first in the set.
	\param c2 The second in the set.*/
	void SwapColumn(std::size_t c1, std::size_t c2)
	{
		rmat::ColumnSwap(m_data, c1, c2);
	}

	/**Scale a row.
	\param r The row.
	\param s The scalar.*/
	void ScaleRow(std::size_t r, const T& s)
	{
		rmat::Scale(m_data[r], s);
	}
	/**Scale a column.
	\param c The column.
	\param s The scalar.*/
	void ScaleColumn(std::size_t c, const T& s)
	{
		for (std::size_t i = 0; i < _R; ++i)
			Set(i, c, Get(i, c)*(s));
	}
	/**Do a complex row operation. [s1*r1 + s2*r2 --> r2]
	\param s1 The value for which to scale the argument row.
	\param r1 The argument row.
	\param s2 The value for which to scale the dest row.
	\param r2 The row that is the destination.*/
	void ComplexRowOp(const T& s1, std::size_t r1, const T& s2, std::size_t r2)
	{
		rmat::RowOp(m_data, r1, s1, r2, s2);
	}
	/**Do a complex col operation. [s1*c1 + s2*c2 --> c2]
	\param s1 The value for which to scale the argument col.
	\param c1 The argument col.
	\param s2 The value for which to scale the dest col.
	\param c2 The col that is the destination.*/
	void ComplexColumnOp(const T& s1, std::size_t c1, const T& s2,
		std::size_t c2)
	{
		rmat::ColumnOp(m_data, c1, s1, c2, s2);
	}

	/**Print a matrix to a string.
	\return The matrix as a string.*/
	std::string ToString() const
	{
		return rmat::ToString(m_data);
	}

	/**Fill a row with a value.  Will fill with assignment (copy or move)
	operator.
	\param v The value to fill.
	\param r The row to fill.*/
	void PushFillRow(std::size_t r, const T& v)
	{
		for (std::size_t i = 0; i < _C; ++i)
		{
			Set(r, i, (v));
		}
	}
	/**Fill a row with a value.  Will fill with assignment (copy or move)
	operator.
	\param args The constructors to do inplace construction apon.
	\param r The row to fill.*/
	template<typename...Args>
	void EmplaceFillRow(std::size_t r, Args&&...args)
	{
		for (std::size_t i = 0; i < _C; ++i)
			Set(r, i, std::forward<Args>(args)...);
	}
	/**Fill a column with a value.  Will fill with assignment (copy or move)
	operator.
	\param v The value to fill.
	\param c The column to fill.*/
	void PushFillColumn(std::size_t c, const T& v)
	{
		for (std::size_t i = 0; i < _R; ++i)
		{
			Set(i, c, (v));
		}
	}
	/**Fill a column with a value.  Will fill with assignment (copy or move)
	operator.
	\param args The constructors to do inplace construction apon.
	\param c The column to fill.*/
	template<typename...Args>
	void EmplaceFillColumn(std::size_t c, Args&&...args)
	{
		for (std::size_t i = 0; i < _R; ++i)
		{
			Set(i, c, std::forward<Args>(args)...);
		}
	}
	/**Fill all with a value.  Will fill with assignment (copy or move)
	operator.
	\param v The value to fill.*/
	void PushFill(const T& v)
	{

		for (std::size_t i = 0; i < _R; ++i)
			for (std::size_t j = 0; j < _C; ++j)
				Set(i, j, (v));
	}
	/**Fill all with a value.  Will fill with assignment (copy or move)
	operator.
	\param args The constructors to do inplace construction apon.*/
	template<typename...Args>
	void EmplaceFill(Args&&...args)
	{
		for (std::size_t i = 0; i < _R; ++i)
			for (std::size_t j = 0; j < _C; ++j)
				Set(i, j, std::forward<Args>(args)...);
	}

	/**Assign a value.
	\param r The row number.
	\param c The column number.
	\param args The args for which to construct the value in place.
	\return a reference to the cell.*/
	template<typename...Args>
	T& Set(std::size_t r, std::size_t c, Args&&...args)
	{
		if (r >= _R || c >= _C)
			throw std::out_of_range("One or more index out of range.");
		T& x = m_data[r][c];
		new (&x) T(std::forward<Args>(args)...);
		return x;
	}
	/**Assign a value.
	\param r The row number.
	\param c The column number.
	\param t The value for which to assign.
	\return a reference to the cell.*/
	T& Set(std::size_t r, std::size_t c, const T& t)
	{
		if (r >= _R || c >= _C)
			throw std::out_of_range("One or more index out of range.");
		T& x = m_data[r][c];
		x = (t);
		return x;
	}
	/**Get a cell.
	\param r The row.
	\param c The column.
	\return A reference to the cell.*/
	T& Get(std::size_t r, std::size_t c)
	{
		return m_data[r][c];
	}
	/**Get a cell.
	\param r The row.
	\param c The column.
	\return A reference to the cell.*/
	const T& Get(std::size_t r, std::size_t c) const
	{
		return m_data[r][c];
	}

	/**Get the transpose of the matrix.
	\return A transpose copy of the matrix.*/
	TransSelfT Transpose() const
	{
		TransSelfT ret;
		ret.m_data = rmat::Transpose(m_data);
		return ret;
	}

	/**********************************************************Static helpers*/

	/**Get A string of a fundamental object.
	\param t The fundamental data to stringify.
	\return The string version.*/
	template<typename T>
	static std::enable_if_t<
		std::is_fundamental<std::decay_t<T>>::value, std::string>
		StringHelper(const T& t)
	{
		return std::to_string((t));
	}
	/**Get A string of a non-fundamental object.  The object type T must have
	`operator std::string()` (const optional) implemented.
	\param t The fundamental data to stringify.
	\return The string version.*/
	template<typename T>
	static std::enable_if_t<
		!std::is_fundamental<std::decay_t<T>>::value, std::string>
		StringHelper(const T& t)
	{
		return (std::string) t;
	}
	/**Create an identity matrix.
	\param s The number for the identity values.
	\tparam Size The size of the matrix (will be square).*/
	template<std::size_t Size, typename T>
	static Matrix<Size, Size, T> IdentityMatrix(const T& s)
	{
		Matrix<Size, Size, T> ret;
		ret.EmplaceFill(T(0));
		for (std::size_t i = 0; i < Size; ++i)
			ret.Set(i, i, (s));
		return ret;
	}
	/**Create a matrix. Example: cg::Matrix<4,3,int>::Create(
		{{1,2,3},{4,5,6},{7,8,9},{10,11,12}});
	\param rows The rows.*/
	static Matrix<_R, _C, T> Create(
		std::initializer_list<std::initializer_list<T>> rows)
	{
		Matrix<_R, _C, T> ret;
		auto s = rows.size();
		if (s != _R)
			throw std::length_error("Incorrect amount of rows.");
		for (std::size_t i = 0; i < s; ++i)
			for (std::size_t j = 0; j < rows.begin()->size(); ++j)
				ret.Set(i, j, *((rows.begin() + i)->begin() + j));
		return ret;
	}
	/***************************************************************Operators*/
	/**Compare for equality.
	\param other The thing to compare.
	\return True if this object has exactly the same values as other.*/
	template<std::size_t _Rw2, std::size_t _Cl2>
	bool operator!=(const Matrix<_Rw2, _Cl2, T>& other)
	{
		static_assert(_Rw2 == _R && _Cl2 == _C, "Cannot compare different "
			"size matrices.");
		for (std::size_t i = 0; i < _R; ++i)
			for (std::size_t j = 0; j < _C; ++j)
				if (Get(i, j) != other.Get(i, j))
					return true;
		return false;
	}
	/**Compare for equality.
	\param other The thing to compare.
	\return True if this object has exactly the same values as other.*/
	template<std::size_t _Rw2, std::size_t _Cl2>
	bool operator==(const Matrix<_Rw2, _Cl2, T>& other)
	{
		return !(*this != other);
	}

	/**Do some matrix math.
	\param other The thing to work with. Must be the same exact size
	for + and -
	\return A reference to this.*/
	SelfT& operator+=(const SelfT& other)
	{
		for (std::size_t i = 0; i < _R; ++i)
			m_data[i] = rmat::Add(m_data[i], other.m_data[i]);
		return *this;
	}
	/**Do some matrix math.
	\param other The thing to work with. Must be the same exact size
	for + and -
	\return A reference to this.*/
	SelfT& operator-=(const SelfT& other)
	{
		for (std::size_t i = 0; i < _R; ++i)
			m_data[i] = rmat::Sub<Columns,DataType>
			(m_data[i], other.m_data[i]);
		return *this;
	}
	/**Do some matrix math.
	\param other The thing to work with. Must be the same exact size
	for + and -
	\return A reference to this.*/
	SelfT& operator*=(const T& other)
	{
		for (std::size_t i = 0; i < _R; ++i)
			rmat::Scale(m_data[i], other);
		return *this;
	}
	/**Do some matrix math.
	\param other The thing to work with. Must be the same exact size
	for + and -
	\return A reference to this.*/
	SelfT& operator/=(const T& other)
	{
		for (std::size_t i = 0; i < _R; ++i)
			cg::rmat::Scale(m_data[i], 1 / other);
		return *this;
	}
	/**Do some matrix math.
	\param other The thing to work with. Must be the same exact size
	for + and -
	\return A reference to this.*/
	SelfT operator*(const T& other) const
	{
		auto copy = *this;
		return copy *= other;
	}
	/**Do some matrix math.
	\param other The thing to work with. Must be the same exact size
	for + and -
	\return A reference to this.*/
	SelfT operator/(const T& other) const
	{
		auto copy = *this;
		return copy /= other;
	}

	/**Get an element.
	\param r The row to get.
	\return A reference to the Row.*/
	cg::rmat::RowMatrix<Columns, DataType>&
		operator[](std::size_t r)
	{
		return (cg::rmat::RowMatrix<Columns, DataType>&) m_data[r];
	}
	/**Get an element.
	\param r The row to get.
	\return A reference to the Row.*/
	const cg::rmat::RowMatrix<Columns,DataType>& 
		operator[](std::size_t r) const
	{
		return (const cg::rmat::RowMatrix<Columns, DataType>&) m_data[r];
	}
	/**Get the negative.
	\return A copy of the matrix as an additive inverse of this matrix.*/
	Matrix<_R, _C, T> operator-() const
	{
		SelfT ret;
		ret.m_data = m_data;
		rmat::Scale<_R,_C,T>(ret.m_data, -1);
		return ret;
	}

protected:
	template<std::size_t a, std::size_t b, typename T>
	friend class Matrix;
	/**The arrays.*/
	cg::rmat::RawMatrix<Rows, Columns, DataType> m_data;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////VECTOR 2///////
///////////////////////////////////////////////////////////////////////////////////////////////////

/**2d Vector class.*/
template<typename DataType>
class Vector2 : public Matrix<2, 1, DataType>
{
public:
	/**The data type to use.*/
	using T = DataType;
	/**The Self type.*/
	using SelfT = Vector2<T>;
	/**The transpose type.*/
	using TransSelfT = typename Matrix<2, 1, T>::TransSelfT;


	/************************************************************Constructors*/

	/**Create the vector with an X and Y value.
	\param x The x value.
	\param y The y value.
	\return The new vector.*/
	static Vector2 Create(const T& x, const T& y)
	{
		Vector2<T> v;
		v.X(x);
		v.Y(y);
		return v;
	}
	/**Create a vector from a 2x1 matrix.
	\param mat The matrix.
	\return The new vector.*/
	static Vector2 Create(const cg::Matrix<2, 1, T>& mat)
	{
		Vector2<T> v;
		std::memcpy(&v.m_data, mat.Begin(), 2 * sizeof(T));
		return v;
	}

	/***************************************************************Operators*/

	/**Copy a 3x1 matrix row to this vector.
	\param mat The matrix.*/
	void operator=(const cg::Matrix<3, 1, T>& mat)
	{
		std::memcpy(&this->m_data, 
            mat.Begin(), 3 * sizeof(T));
	}
	/***************************************************************Utilities*/

	/**Set the X value.
	\param args The args for which to construct X.*/
	void X(const T& args)
	{
		Matrix<2, 1, T>::Set(0, 0, (args));
	}
	/**Set the Y value.
	\param args The args for which to construct Y.*/
	void Y(const T& args)
	{
		Matrix<2, 1, T>::Set(1, 0, (args));
	}
	/**Get the X value.
	\return The X value as a const reference.*/
	const T& X() const
	{
		return Matrix<2, 1, DataType>::Get(0, 0);
	}
	/**Get the X value.
	\return The X value as a const reference.*/
	const T& Y() const
	{
		return Matrix<2, 1, DataType>::Get(1, 0);
	}
	/**Set the values.
	\param x The x value.
	\param y The y value.*/
	void Set(const T& x, const T& y)
	{
		X((x));
		Y((y));
	}
private:

};

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////VECTOR 3///////
///////////////////////////////////////////////////////////////////////////////////////////////////

/**3d Vector class.*/
template<typename DataType>
class Vector3 : public Matrix<3, 1, DataType>
{
public:
	/**The data type to use.*/
	using T = DataType;
	/**The Self type.*/
	using SelfT = Vector3<T>;
	/**The transpose type.*/
	using TransSelfT = typename Matrix<3, 1, T>::TransSelfT;


	/************************************************************Constructors*/
	/**Create the vector with an X and Y value.
	\param x The x value.
	\param y The y value.
	\param z The z value.
	\return The new vector.*/
	static Vector3 Create(const T& x, const T& y, const T& z)
	{
		Vector3<T> v;
		v.X((x));
		v.Y((y));
		v.Z((z));
		return v;
	}
	/**Create a vector from a 3x1 matrix.
	\param mat The matrix.
	\return The new vector.*/
	static Vector3 Create(const cg::Matrix<3, 1, T>& mat)
	{
		Vector3<T> v;
		std::memcpy(&v.m_data, mat.Begin(), 3 * sizeof(T));
		return v;
	}

	/***************************************************************Operators*/

	/**Copy a 3x1 matrix row to this vector.
	\param mat The matrix.*/
	void operator=(const cg::Matrix<3, 1, T>& mat)
	{
		std::memcpy(&this->m_data, mat.Begin(), 3 * sizeof(T));
	}

	/***************************************************************Utilities*/

	/**Set the X value.
	\param args The args for which to construct X.*/
	void X(const T& args)
	{
		Matrix<3, 1, T>::Set(0, 0, (args));
	}
	/**Set the Y value.
	\param args The args for which to construct Y.*/
	void Y(const T& args)
	{
		Matrix<3, 1, T>::Set(1, 0, (args));
	}
	/**Set the Z value.
	\param args The args for which to construct Z.*/
	void Z(const T& args)
	{
		Matrix<3, 1, T>::Set(2, 0, (args));
	}
	/**Get the X value.
	\return The X value as a const reference.*/
	const T& X() const
	{
		return Matrix<3, 1, DataType>::Get(0, 0);
	}
	/**Get the X value.
	\return The X value as a const reference.*/
	const T& Y() const
	{
		return Matrix<3, 1, DataType>::Get(1, 0);
	}
	/**Get the Z value.
	\return The Z value as a const reference.*/
	const T& Z() const
	{
		return Matrix<3, 1, DataType>::Get(2, 0);
	}
	/**Set the values.
	\param x The x value.
	\param y The y value.
	\param z The z value.*/
	void Set(const T& x, const T& y, const T& z)
	{
		X((x));
		Y((y));
		Z((z));
	}
private:

};

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////VECTOR 4///////
///////////////////////////////////////////////////////////////////////////////////////////////////

/**4d Vector class.*/
template<typename DataType>
class Vector4 : public Matrix<4, 1, DataType>
{
public:
	/**The data type to use.*/
	using T = DataType;
	/**The Self type.*/
	using SelfT = Vector4<T>;
	/**The transpose type.*/
	using TransSelfT = typename Matrix<4, 1, T>::TransSelfT;


	/************************************************************Constructors*/

	/**Create the vector with an X and Y value.
	\param x The x value.
	\param y The y value.
	\param z The z value.
	\param w The w value
	\return The new vector.*/
	static Vector4 Create(const T& x, const T& y, const T& z, const T& w)
	{
		Vector4<T> v;
		v.X((x));
		v.Y((y));
		v.Z((z));
		v.W((w));
		return v;
	}
	/**Create a vector from a 4x1 matrix.
	\param mat The matrix.
	\return The new vector.*/
	static Vector4 Create(const cg::Matrix<4, 1, T>& mat)
	{
		Vector4<T> v;
		std::memcpy(&v.m_data, mat.Begin(), 4 * sizeof(T));
		return v;
	}

	/**********************************************************Static Helpers*/

	/**Get the X-AXIS/Roll
	\param amt The amount in radians to rotate.
	\return A prebuild matrix that when multiplied to, will produce a rotation
	around the X axis.*/
	static Matrix<4, 4, T> RotXMatrixR(const T& amt)
	{
		auto rad = amt;
		auto mat = Matrix<4, 4, T>::Create({
			{ 1, 0, 0,0 },
			{ 0,(T)std::cos(rad),(T)-std::sin(rad),0 },
			{ 0,(T)std::sin(rad),(T)std::cos(rad),0 },
			{ 0,0,0,1 }
		});
		return mat;
	}
	/**Get the X-AXIS/Roll
	\param amt The amount in degrees to rotate.
	\return A prebuild matrix that when multiplied to, will produce a rotation
	around the X axis.*/
	static Matrix<4, 4, T> RotXMatrixD(const T& amt)
	{
		return RotXMatrixR(amt*(T)DegreeToRadian);
	}
	/**Get the Y-AXIS/Pitch
	\param amt The amount in radians to rotate.
	\return A prebuild matrix that when multiplied to, will produce a rotation
	around the Y axis.*/
	static Matrix<4, 4, T> RotYMatrixR(const T& amt)
	{
		auto rad = amt;
		auto mat = Matrix<4, 4, T>::Create({
			{ (T)std::cos(rad),0,(T)std::sin(rad),0 },
			{ 0,1,0,0 },
			{ (T)-std::sin(rad),0,(T)std::cos(rad),0 },
			{ 0,0,0,1 }
		});
		return mat;
	}
	/**Get the Y-AXIS/Pitch
	\param amt The amount in degrees to rotate.
	\return A prebuild matrix that when multiplied to, will produce a rotation
	around the Y axis.*/
	static Matrix<4, 4, T> RotYMatrixD(const T& amt)
	{
		return RotYMatrixR(amt* (T)DegreeToRadian);
	}
	/**Get the Z-AXIS/Yaw
	\param amt The amount in radians to rotate.
	\return A prebuild matrix that when multiplied to, will produce a rotation
	around the Z axis.*/
	static Matrix<4, 4, T> RotZMatrixR(const T& amt)
	{
		auto rad = amt;
		auto mat = Matrix<4, 4, T>::Create({
			{ (T)std::cos(rad),(T)-std::sin(rad),0,0 },
			{ (T)std::sin(rad),(T)std::cos(rad),0,0 },
			{ 0,0, 1,0 },
			{ 0,0,0,1 }
		});
		return mat;
	}
	/**Get the Z-AXIS/Yaw
	\param amt The amount in degrees to rotate.
	\return A prebuild matrix that when multiplied to, will produce a rotation
	around the Z axis.*/
	static Matrix<4, 4, T> RotZMatrixD(const T& amt)
	{
		return RotZMatrixR(amt* (T)DegreeToRadian);
	}

	/**Create a translation matrix.
	\param x The X amount to translate.
	\param y The Y amount to translate.
	\param z The Z amount to translate.
	\return A matrix that when multiplied to the vector will translate it.*/
	static Matrix<4, 4, T> TranslationMatrix(const T& x, const T& y,
		const T& z)
	{
		auto mat = Matrix<4, 4, T>::Create({
			{ 1,0,0,x },
			{ 0,1,0,y },
			{ 0,0,1,z },
			{ 0,0,0,1 }
		});
		return mat;
	}
	/**Create a scaling matrix.
	\param x The X amount to scale.
	\param y The Y amount to scale.
	\param z The Z amount to scale.
	\return A matrix that when multiplied to the vector will scale it.*/
	static Matrix<4, 4, T> ScalarMatrix(const T& x, const T& y,
		const T& z)
	{
		auto mat = Matrix<4, 4, T>::Create({
			{ x,0,0,0 },
			{ 0,y,0,0 },
			{ 0,0,z,0 },
			{ 0,0,0,1 }
		});
		return mat;
	}
	/**Create a Rotation matrix in degrees.
	\param x The rotation about X.
	\param y The rotation about Y.
	\param z The rotation about Z.
	\return A matrix that when multiplied to the vector will rotate it.*/
	static Matrix<4, 4, T> RotationMatrixD(const T& x, const T& y,
		const T& z)
	{
		return RotXMatrixD(x) * RotYMatrixD(y) * RotZMatrixD(z);
	}
	/**Create a Rotation matrix in radians.
	\param x The rotation about X.
	\param y The rotation about Y.
	\param z The rotation about Z.
	\return A matrix that when multiplied to the vector will rotate it.*/
	static Matrix<4, 4, T> RotationMatrixR(const T& x, const T& y,
		const T& z)
	{
		return RotXMatrixR(x) * RotYMatrixR(y) * RotZMatrixR(z);
	}

	/***************************************************************Operators*/

	/**Copy a 4x1 matrix row to this vector.
	\param mat The matrix.*/
	void operator=(const cg::Matrix<4, 1, T>& mat)
	{
		std::memcpy(&this->m_data, mat.Begin(), 4 * sizeof(T));
	}

	/***************************************************************Utilities*/

	/**Set the X value.
	\param args The args for which to construct X.*/
	void X(const T& args)
	{
		Matrix<4, 1, T>::Set(0, 0, (args));
	}
	/**Set the Y value.
	\param args The args for which to construct Y.*/
	void Y(const T& args)
	{
		Matrix<4, 1, T>::Set(1, 0, (args));
	}
	/**Set the Z value.
	\param args The args for which to construct Z.*/
	void Z(const T& args)
	{
		Matrix<4, 1, T>::Set(2, 0, (args));
	}
	/**Set the W value.
	\param args The args for which to construct W.*/
	void W(const T& args)
	{
		Matrix<4, 1, T>::Set(3, 0, (args));
	}
	/**Get the X value.
	\return The X value as a const reference.*/
	const T& X() const
	{
		return Matrix<4, 1, DataType>::Get(0, 0);
	}
	/**Get the X value.
	\return The X value as a const reference.*/
	const T& Y() const
	{
		return Matrix<4, 1, DataType>::Get(1, 0);
	}
	/**Get the Z value.
	\return The Z value as a const reference.*/
	const T& Z() const
	{
		return Matrix<4, 1, DataType>::Get(2, 0);
	}
	/**Get the W value.
	\return The W value as a const reference.*/
	const T& W() const
	{
		return Matrix<4, 1, DataType>::Get(3, 0);
	}
	/**Set the values.
	\param x The x value.
	\param y The y value.
	\param z The z value.
	\param w The w value.*/
	void Set(const T& x, const T& y, const T& z, const T& w)
	{
		X((x));
		Y((y));
		Z((z));
		W((w));
	}
	/**Translate this vector assuming its a set of Homogeneous cordinates.
	\param x The x translation.
	\param y The y translation.
	\param z The Z translation.
	\return A reference to this.*/
	SelfT& HTranslate(const T& x, const T& y, const T& z)
	{
		X(X() + (x));
		Y(Y() + (y));
		Z(Z() + (z));
		return *this;
	}
	/**Scale this vector assuming its a set of Homogeneous cordinates.
	\param x The x translation.
	\param y The y translation.
	\param z The Z translation.
	\return A reference to this.*/
	SelfT& HScale(const T& x, const T& y, const T& z)
	{
		X(X() * (x));
		Y(Y() * (y));
		Z(Z() * (z));
		return *this;
	}
	/**Rotate this around the X axis.
	\param amt The amount to rotate in degrees.
	\return A reference to this.*/
	SelfT& RotAboutXD(const T& amt)
	{
		*this = RotXMatrixD(amt) * *this;
		return *this;
	}
	/**Rotate this around the X axis.
	\param amt The amount to rotate in radians.
	\return A reference to this.*/
	SelfT& RotAboutXR(const T& amt)
	{
		*this = RotXMatrixR(amt) * *this;
		return *this;
	}
	/**Rotate this around the X axis.
	\param amt The amount to rotate in degrees.
	\return A reference to this.*/
	SelfT& RollD(const T& amt) {
		return RotAboutXD(amt);
	}
	/**Rotate this around the X axis.
	\param amt The amount to rotate in radians.
	\return A reference to this.*/
	SelfT& RollR(const T& amt) {
		return RotAboutXR(amt);
	}
	/**Rotate this around the Y axis.
	\param amt The amount to rotate in degrees.
	\return A reference to this.*/
	SelfT& RotAboutYD(const T& amt)
	{
		*this = RotYMatrixD(amt) * *this;
		return *this;
	}
	/**Rotate this around the Y axis.
	\param amt The amount to rotate in radians.
	\return A reference to this.*/
	SelfT& RotAboutYR(const T& amt)
	{
		*this = RotYMatrixR(amt) * *this;
		return *this;
	}
	/**Rotate this around the Y axis.
	\param amt The amount to rotate in degrees.
	\return A reference to this.*/
	SelfT& PitchD(const T& amt) {
		return RotAboutYD(amt);
	}
	/**Rotate this around the Y axis.
	\param amt The amount to rotate in radians.
	\return A reference to this.*/
	SelfT& PitchR(const T& amt) {
		return RotAboutYR(amt);
	}
	/**Rotate this around the Z axis.
	\param amt The amount to rotate in degrees.
	\return A reference to this.*/
	SelfT& RotAboutZD(const T& amt)
	{
		*this = RotZMatrixD(amt) * *this;
		return *this;
	}
	/**Rotate this around the Z axis.
	\param amt The amount to rotate in radians.
	\return A reference to this.*/
	SelfT& RotAboutZR(const T& amt)
	{
		*this = RotZMatrixR(amt) * *this;
		return *this;
	}
	/**Rotate this around the Z axis.
	\param amt The amount to rotate in degrees.
	\return A reference to this.*/
	SelfT& YawD(const T& amt) {
		return RotAboutZD(amt);
	}
	/**Rotate this around the Z axis.
	\param amt The amount to rotate in radians.
	\return A reference to this.*/
	SelfT& YawR(const T& amt) {
		return RotAboutZR(amt);
	}
	/**Get this vector crossed with another. W will stay as this vectors W val.
	\param v The other vector.
	\return The result of this CROSS v*/
	SelfT Cross(const SelfT& v) const
	{
		return SelfT::Create(
			/*X*/Y()*v.Z() - Z()*v.Y(),
			/*Y*/Z()*v.X() - X()*v.Z(),
			/*Z*/X()*v.Y() - Y()*v.X(),
			/*W*/W()
		);
	}
private:

};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////Spherical Vec///
///////////////////////////////////////////////////////////////////////////////////////////////////

/**A class for a homogenous spherical vector.*/
template<typename DataType>
class SphericalVector : public Matrix<4, 1, DataType>
{
public:

	/**The data type to use.*/
	using T = DataType;
	/**The Self type.*/
	using SelfT = SphericalVector<T>;
	/**The transpose type.*/
	using TransSelfT = typename Matrix<4, 1, T>::TransSelfT;


	/************************************************************Constructors*/

	/**Create the vector with pre set values.
	\param r The RHO value.
	\param t The THETA value.
	\param p The PHI value.
	\param w The w value (projective element)*/
	static SphericalVector Create(const T& r, const T& t, const T& p,
		const T& w)
	{
		SphericalVector<T> v;
		v.R((Matrix<2, 1, DataType>::x));
		v.Th((Matrix<2, 1, DataType>::y));
		v.P((Matrix<2, 1, DataType>::z));
		v.W((Matrix<2, 1, DataType>::w));
		return v;
	}


	/***************************************************************Utilities*/

	/**Set the RHO value.
	\param v The value to set.*/
	void R(const T& v)
	{
		Matrix<2, 1, DataType>::Set(0, 0, v);
	}
	/**Get the RHO value.
	\return The RHO value as a const reference.*/
	const T& R() const
	{
		return Matrix<2, 1, DataType>::Get(0, 0);
	}
	/**Set the THETA value.
	\param v The value to set.*/
	void Th(const T& v)
	{
		Matrix<2, 1, DataType>::Set(1, 0, v);
	}
	/**Get the THETA value.
	\return The THETA value as a const reference.*/
	const T& Th() const
	{
		return Matrix<2, 1, DataType>::Get(1, 0);
	}
	/**Set the PHI value.
	\param v The value to set.*/
	void P(const T& v)
	{
		Matrix<2, 1, DataType>::Set(2, 0, v);
	}
	/**Get the PHI value.
	\return The PHI value as a const reference.*/
	const T& P() const
	{
		return Matrix<2, 1, DataType>::Get(2, 0);
	}
	/**Set the W value.
	\param args The args for which to construct W.*/
	void W(const T& args)
	{
		Matrix<4, 1, T>::Set(3, 0, (args));
	}
	/**Get the W value.
	\return The W value as a const reference.*/
	const T& W() const
	{
		return Matrix<2, 1, DataType>::Get(3, 0);
	}
private:

};

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////GLOBAL OPERATORS///
///////////////////////////////////////////////////////////////////////////////////////////////////

/**Do math with scalar first.
\param scalar The scalar.
\param mat The matrix.*/
template<std::size_t Rows, std::size_t Columns, typename DataType>
Matrix<Rows, Columns, DataType> operator*(const DataType& scalar,
	const Matrix<Rows, Columns, DataType>& mat)
{
	return mat * scalar;
}

/**Multiply a single column by a single row.
\param row The row.
\param col The column.*/
template<std::size_t S, typename T>
T operator*(const Matrix<1, S, T>& row, const Matrix<S, 1, T>& col)
{
	T ret = 0;
	for (std::size_t i = 0; i < S; ++i)
		ret += row[0][i] * col[i][0];
	return ret;
}
/**Multiply two matrices.  See matrix multiplication rules. arguments must be
is the proper order.
\param lhs The first.
\param rhs The second.*/
template<std::size_t S1, std::size_t S2, std::size_t S3, typename T>
Matrix<S1, S3, T> operator*(const Matrix<S1, S2, T>& lhs,
	const Matrix<S2, S3, T>& rhs)
{
	Matrix<S1, S3, T> ret;
	const static std::size_t R = S1;
	const static std::size_t C = S3;
	const static std::size_t multAmt = S1 * S3;
	for (std::size_t i = 0; i < R; ++i)
		for (std::size_t j = 0; j < C; ++j)
		{
			ret[i][j] = lhs.Row(i) * rhs.Column(j);
		}
	return ret;
}

/**Do dot product on two same size vectors.
\param v1 The first vector.
\param v2 The second vector.
\return The DOT product as T.*/
template<typename T>
T operator*(const cg::Vector2<T>& v1, const cg::Vector2<T>& v2)
{
	return v2.Transpose()  * v1;
}
/**Do dot product on two same size vectors.
\param v1 The first vector.
\param v2 The second vector.
\return The DOT product as T.*/
template<typename T>
T operator*(const cg::Vector3<T>& v1, const cg::Vector3<T>& v2)
{
	return v2.Transpose() *v1;
}
/**Do dot product on two same size vectors.
\param v1 The first vector.
\param v2 The second vector.
\return The DOT product as T.*/
template<typename T>
T operator*(const cg::Vector4<T>& v1, const cg::Vector4<T>& v2)
{
	return v2.Transpose()* v1;
}




template<typename T>
using Matrix3 = Matrix<3, 3, T>;
template<typename T>
using Matrix2 = Matrix<2, 2, T>;

template class Vector2<float>;
template class Vector3<float>;
template class Vector4<float>;
template class Matrix<10, 10, int>;
template class Matrix<10, 10, float>;



}