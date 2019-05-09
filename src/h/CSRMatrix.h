#pragma once
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>
#include "SparseMatrix.h"


template <typename T> class CSRMatrix : public SparseMatrix<T>
{
private:
	std::vector<int>* _tempColumns;
public:
	int size;
	int elementCount;
	T* values;
	int* columns;
	int* rowIndexs;

	CSRMatrix(int m) : SparseMatrix<T>(m)
	{
		size = m;
		values = nullptr;
		columns = nullptr;
		rowIndexs = new int[size + 1];
		_tempColumns = nullptr;
	}

	// call beginPositionMark() will allocate _tempColumns,
	// so that markPostion() can be called.
	void beginPostionMark()
	{
		_tempColumns = new std::vector<int>[size];

	}

	void markPosition(int row, int column)
	{
		// insert column
		_tempColumns[row - 1].push_back(column);
	}

	void Allocate()
	{

		for (int row = 0; row < size; ++row)
		{
			std::vector<int>& v = _tempColumns[row];
			std::sort(v.begin(), v.end());
			v.erase(std::unique(v.begin(), v.end()), v.end());
		}

		// first, calculate row indexs. Notice row index starts from 1
		rowIndexs[0] = 1;
		for (int row = 0; row < size; ++row)
		{
			rowIndexs[row + 1] = rowIndexs[row] + int(_tempColumns[row].size());
		}

		// allocate columns
		elementCount = rowIndexs[size] - rowIndexs[0];
		columns = new int[elementCount];

		// write columns
		int count = 0;
		for (int row = 0; row < size; ++row)
		{
			// column already sorted in _tempColumns
			for (const auto& column : _tempColumns[row])
			{
				columns[count] = column;
				count++;
			}
		}

		// free _tempColumns
		delete[] _tempColumns;
		_tempColumns = nullptr;

		// allocate values last to save memory
		values = new T[elementCount];
		for (int i = 0; i < elementCount; ++i)
			values[i] = T(0);
	}

	// get item at (row, column)
	T& operator()(unsigned row, unsigned column)
	{
		if (row > column)
		{
			return this->operator()(column, row);
		}
		int offset1 = rowIndexs[row - 1] - 1;
		int offset2 = rowIndexs[row] - 1;
		// index lies in [offset1, offset2)
		// find index by bisection

		while (offset2 != (offset1 + 1))
		{
			int offset = (offset1 + offset2) / 2;
			if (columns[offset] > int(column))
			{
				offset2 = offset;
			}
			else
			{
				offset1 = offset;
			}
		}

		return values[offset1];
	}

	int dim() const { return size; }

	void Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND)
	{
		//  Assemble global stiffness matrix
		for (unsigned int j = 0; j < ND; j++)
		{
			unsigned int Lj = LocationMatrix[j];    // Global equation number corresponding to jth DOF of the element
			if (!Lj)
				continue;

			//      Address of diagonal element of column j in the one dimensional element stiffness matrix
			unsigned int DiagjElement = (j + 1)*j / 2 + 1;

			for (unsigned int i = 0; i <= j; i++)
			{
				unsigned int Li = LocationMatrix[i];    // Global equation number corresponding to ith DOF of the element

				if (!Li)
					continue;

				(*this)(Li, Lj) += Matrix[DiagjElement + j - i - 1];
			}
		}

		return;
	};
	~CSRMatrix()
	{
		delete[] rowIndexs;
		delete[] values;
		delete[] columns;
	}

};

template <typename T> std::ostream& operator<<(std::ostream& out, const CSRMatrix<T>& mat)
{
	out << "CSR Matrix, size = " << mat.size << std::endl;
	out << "values = " << std::endl << "(";
	for (int i = 0; i < mat.elementCount; ++i)
	{
		out << std::setw(14) << mat.values[i];
	}
	out << ")\ncolumns = \n(";
	for (std::size_t i = 0; i < mat.elementCount; ++i)
	{
		out << std::setw(14) << mat.columns[i];
	}
	out << ")\nrowIndexs = \n(";
	for (std::size_t i = 0; i <= mat.size; ++i)
	{
		out << std::setw(14) << mat.rowIndexs[i];
	}
	out << ")";
	return out;
}
