/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Solver.h"
#include <cmath>
#include <iostream>
#ifdef _PARDISO_
#include<mkl.h>
#include <algorithm>
#include <cfloat>

#endif // _PARDISO_

using namespace std;

// LDLT facterization
void CLDLTSolver::LDLT()
{
	unsigned int N = K->dim();
    unsigned int* ColumnHeights = K->GetColumnHeights();   // Column Hights

	for (unsigned int j = 2; j <= N; j++)      // Loop for column 2:n (Numbering starting from 1)
	{
        // Row number of the first non-zero element in column j (Numbering starting from 1)
		unsigned int mj = j - ColumnHeights[j-1];
        
		for (unsigned int i = mj+1; i <= j-1; i++)	// Loop for mj+1:j-1 (Numbering starting from 1)
		{
            // Row number of the first nonzero element in column i (Numbering starting from 1)
			unsigned int mi = i - ColumnHeights[i-1];

			double C = 0.0;
			for (unsigned int r = max(mi, mj); r <= i-1; r++)
				C += (*K)(r,i) * (*K)(r,j);		// C += L_ri * U_rj

			(*K)(i,j) -= C;	// U_ij = K_ij - C
		}

		for (unsigned int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)
		{
			double Lrj = (*K)(r,j) / (*K)(r,r);	// L_rj = U_rj / D_rr
			(*K)(j,j) -= Lrj * (*K)(r,j);	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			(*K)(r,j) = Lrj;
		}

        if (fabs((*K)(j,j)) <= FLT_MIN)
        {
            cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
            	 << "    Euqation no = " << j << endl
            	 << "    Pivot = " << (*K)(j,j) << endl;
            
            exit(4);
        }
    }
};

// Solve displacement by back substitution
void CLDLTSolver::BackSubstitution(double* Force)
{
	unsigned int N = K->dim();
    unsigned int* ColumnHeights = K->GetColumnHeights();   // Column Hights

//	Reduce right-hand-side load vector (LV = R)
	for (unsigned int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
	{
        unsigned int mi = i - ColumnHeights[i-1];

		for (unsigned int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1
			Force[i-1] -= (*K)(j,i) * Force[j-1];	// V_i = R_i - sum_j (L_ji V_j)
	}

//	Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
	for (unsigned int i = 1; i <= N; i++)	// Loop for i=1:N
		Force[i-1] /= (*K)(i,i);	// Vbar = D^(-1) V

	for (unsigned int j = N; j >= 2; j--)	// Loop for j=N:2
	{
        unsigned int mj = j - ColumnHeights[j-1];

		for (unsigned int i = mj; i <= j-1; i++)	// Loop for i=mj:j-1
			Force[i-1] -= (*K)(i,j) * Force[j-1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)
	}
};

#ifdef _PARDISO_
void CSRSolver::solve(double* Force, unsigned NLCase)
{
	void* pt[64];
	for (unsigned _ = 0; _ < 64; _++) pt[_] = 0;
	int num_procs;
	const int mtype = 2;
	int iparm[64] = { 0 };
	pardisoinit(pt, &mtype, iparm);
	
	iparm[1] = 2; // The parallel (OpenMP) version of the nested dissection algorithm.
	// Try 3 later
	iparm[5] = 1; // write back to Force
	char    *var;
	num_procs = 4;
	iparm[2] = num_procs;
	
	const int maxfct(1),mnum(1);
	const int size = K->size;
	double* values = K->values;
	int* column_index = K->col_index;
	int* row_index = K->row_index;

	const int nrhs = NLCase;
	double* rhs = Force;

	int phase = 13;
	double* res = new double[nrhs*size];
	for (std::size_t _ = 0; _ < nrhs*size; _++) res[_] = 0;

	int msglvl = 0; // print info

	int* perm = new int[size];
	int error;

	pardiso(
		pt, // handle to some shit
		&maxfct, // maxfct
		&mnum, // mnum
		&mtype, // sym pos matrix
		&phase, // go through all
		&size,  // size of matrix
		values,
		row_index,
		column_index,
		perm, 
		&nrhs,
		iparm, 
		&msglvl, // print info or not
		rhs,
		res,
		&error // see if any error
	);
	if (error)
	{
		std::cerr << "ERROR IN PARDISO SOLVER: " << error << std::endl;
		exit(8);
	}
	delete[] perm;
	delete[] res;
}
#endif // _PARDISO_



