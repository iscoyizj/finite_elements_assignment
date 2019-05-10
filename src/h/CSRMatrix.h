#pragma once
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>



template <typename T> class CSparseMatrix
{
private:
	std::vector<int>* col_foralloc;
	
public:
	int size;
	int num_ele;
	T* values;
	int* col_index;
	int* row_index;

	CSparseMatrix(int n)
	{
		size = n;
		values = nullptr;
		col_index = nullptr;
		row_index = new int[size + 1];
		col_foralloc = new std::vector<int>[size];
	}

	void markPosition(int row, int column)
	{
		col_foralloc[row - 1].push_back(column);
	}

	void Allocate()
	{
		for (int row = 0; row < size; ++row)
		{
			std::vector<int>& v = col_foralloc[row];
			std::sort(v.begin(), v.end());
			v.erase(std::unique(v.begin(), v.end()), v.end());
		}
		row_index[0] = 1;
		for (int row = 0; row < size; ++row)
		{
			row_index[row + 1] = row_index[row] + int(col_foralloc[row].size());
		}
		// allocate columns
		num_ele = row_index[size] - row_index[0];
		col_index = new int[num_ele];

		// write columns
		int count = 0;
		for (int row = 0; row < size; ++row)
		{
			for (int column : col_foralloc[row])
			{
				col_index[count] = column;
				count++;
			}
		}
		delete[] col_foralloc;
		col_foralloc = nullptr;

		// allocate values last to save memory
		values = new T[num_ele];
		for (int i = 0; i < num_ele; ++i)
			values[i] = T(0);
	}

	T& operator()(int row, int column)
	{
		if (row > column)
		{
			return this->operator()(column, row);
		}
		int offset1 = row_index[row - 1] - 1;
		int offset2 = row_index[row] - 1;
	

		while (offset2 != (offset1 + 1))
		{
			int offset = (offset1 + offset2) / 2;
			if (col_index[offset] > column)
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

//	int dim() const { return size; }

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
	~CSparseMatrix()
	{
		delete[] row_index;
		delete[] values;
		delete[] col_index;
	}

};

