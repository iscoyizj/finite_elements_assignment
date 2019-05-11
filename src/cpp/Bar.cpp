/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Bar.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CBar::CBar()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 6;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CBar::~CBar()
{
}

//	Read element data from stream Input
bool CBar::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl 
			 << "    Expected element : " << Ele + 1 << endl
			 << "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
    ElementMaterial_ = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CBar::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CBar::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int CBar::SizeOfStiffnessMatrix() { return 21; }

//Caculate Gravity of Elements
void CBar::GravityCalculation(double* ptr_force)
{
	double g = 9.8;
	double DX[3];
	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
	}
	double leng = 0;
	for (unsigned int i = 0; i < 3; i++)
	{
		leng += DX[i] * DX[i];
	}
	leng = sqrt(leng);
	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element
	weight = g * leng * material_->Area * material_->density;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBar::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate bar length
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[6];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];
	DX2[3] = DX[0] * DX[1];
	DX2[4] = DX[1] * DX[2];
	DX2[5] = DX[0] * DX[2];

	double L2 = DX2[0] + DX2[1] + DX2[2];
	double L = sqrt(L2);

//	Calculate element stiffness matrix

	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double k = material_->E * material_->Area / L / L2;

	Matrix[0] = k*DX2[0];
	Matrix[1] = k*DX2[1];
	Matrix[2] = k*DX2[3];
	Matrix[3] = k*DX2[2];
	Matrix[4] = k*DX2[4];
	Matrix[5] = k*DX2[5];
	Matrix[6] = k*DX2[0];
	Matrix[7] = -k*DX2[5];
	Matrix[8] = -k*DX2[3];
	Matrix[9] = -k*DX2[0];
	Matrix[10] = k*DX2[1];
	Matrix[11] = k*DX2[3];
	Matrix[12] = -k*DX2[4];
	Matrix[13] = -k*DX2[1];
	Matrix[14] = -k*DX2[3];
	Matrix[15] = k*DX2[2];
	Matrix[16] = k*DX2[4];
	Matrix[17] = k*DX2[5];
	Matrix[18] = -k*DX2[2];
	Matrix[19] = -k*DX2[4];
	Matrix[20] = -k*DX2[5];
}


void CBar::ElementMass(double* Mass)
{
	clear(Mass, SizeOfStiffnessMatrix());
	//	Calculate bar length
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[6];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];
	DX2[3] = DX[0] * DX[1];
	DX2[4] = DX[1] * DX[2];
	DX2[5] = DX[0] * DX[2];

	double L2 = DX2[0] + DX2[1] + DX2[2];
	double L = sqrt(L2);

	DX2[0] = DX2[0] / L2;
	DX2[1] = DX2[1] / L2;
	DX2[2] = DX2[2] / L2;
	DX2[3] = DX2[3] / L2;
	DX2[4] = DX2[4] / L2;
	DX2[5] = DX2[5] / L2;

	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double area = material_->Area;
	double density = material_->density;
	double m = density * area*L;

	Mass[0] = (1 - 2 * DX2[0] / 3)*m;
	Mass[1] = -2 * DX2[1] / 3 * m;
	Mass[2] = -2 * DX2[3] / 3 * m;
	Mass[3] = (1 - 2 * DX2[2] / 3)*m;
	Mass[4] = -2 * DX2[4] / 3 * m;
	Mass[5] = -2 * DX2[5] / 3 * m;
	Mass[6] = (1 - 2 * DX2[0] / 3)*m;
	Mass[7] = -5 * DX2[5] / 6 * m;
	Mass[8] = -5 * DX2[3] / 6 * m;
	Mass[9] = (1 - 5 * DX2[0] / 6)*m;
	Mass[10] = -2 * DX2[1] / 3 * m;
	Mass[11] = -2 * DX2[3] / 3 * m;
	Mass[12] = -5 * DX2[4] / 6 * m;
	Mass[13] = (1 - 5 * DX2[1] / 6)*m;
	Mass[14] = -5 * DX2[3] / 6 * m;
	Mass[15] = -2 * DX2[4] / 3 * m;
	Mass[16] = -2 * DX2[5] / 3 * m;
	Mass[17] = (1 - 2 * DX2[0] / 3)*m;
	Mass[18] = -5 * DX2[4] / 6 * m;
	Mass[19] = (1 - 5 * DX2[2] / 6)*m;
	Mass[20] = -5 * DX2[5] / 6 * m;
}

//	Calculate element stress 
void CBar::ElementStress(double* stress, double* Displacement)
{
	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i+3] = -S[i];
	}
	
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += S[i] * Displacement[LocationMatrix_[i]-1];
	}
}

void CBar::RecoverElementStress(double* Displacement, double* A)
{

	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i] * DX[i];
	}
	int N1;
	int N2;


	N1 = nodes_[0]->NodeNumber - 1;
	N2 = nodes_[1]->NodeNumber - 1;

	A[73 * N1 + 72] += 1;
	A[73 * N2 + 72] += 2;

	for (unsigned int l = 0; l < 3; l++)
	{
		A[73 * N1 + l] = 0;
		A[73 * N2 + l] = 0;

	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i + 3] = -S[i];
	}

	A[72 + 73 * N1] += 1;
	A[72 + 73 * N2] += 1;

	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			A[73 * N1 + 3] += S[i] * Displacement[LocationMatrix_[i] - 1];
		A[73 * N2 + 3] += S[i] * Displacement[LocationMatrix_[i] - 1];
	}

	nodes_[0]->stress_node[0] = A[73 * N1 + 3];
	nodes_[1]->stress_node[0] = A[73 * N2 + 3];

}
//	Calculate element stress for plot
void CBar::ElementStressplot1(double* newlocation, double* Displacement)
{
	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element
	for (unsigned int j = 0; j < NEN_; j++)
	{
		for (unsigned int i = 0; i < 3; i++)
		{
			if (LocationMatrix_[i + 3 * j])
				newlocation[i + 9 * j] = nodes_[j]->XYZ[i] + Displacement[LocationMatrix_[i + 3 * j] - 1];
			else
				newlocation[i + 9 * j] = nodes_[j]->XYZ[i];

		}

		for (unsigned int i = 3; i < 9; i++)
		{
			newlocation[i + 9 * j] = nodes_[j]->stress_node[i - 3];
		}
	}
}