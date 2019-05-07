/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Q4.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CQ4::CQ4()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 8;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CQ4::~CQ4()
{
}

//	Read element data from stream Input
bool CQ4::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2, N3, N4;	// four node numbers

	Input >> N1 >> N2 >> N3 >> N4>> MSet;
    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[0]->bcode[2]=1;
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[1]->bcode[2]=1;
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[2]->bcode[2]=1;
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[3]->bcode[2]=1;

	return true;
}

//	Write element data to stream
void CQ4::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber 
		   << setw(9) << nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CQ4::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 2; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
	
/*	printf("LM\n");
	for (int ii=0;ii<8;ii++)
		printf("%d\n",LocationMatrix_[ii]);
*/
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 4Q element, element stiffness is a 8x8 matrix, whose upper triangular part
//	has 36 elements
unsigned int CQ4::SizeOfStiffnessMatrix() { return 36; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element

void CQ4::CalculateBe (double* Be, double eta, double psi)
{
	
//GN matrix
    double GN[4];
	
	GN[0] = (eta-1)/4;
	GN[1] = (eta+1)/4;
	GN[2] = (psi-1)/4;
	GN[3] = (psi+1)/4;
	
//Jaccobi
// [0  2]
// [1  4]
	double Ja[4];
	Ja[0] = GN[0]*nodes_[0]->XYZ[0] - GN[0]*nodes_[1]->XYZ[0] + GN[1]*nodes_[2]->XYZ[0] - GN[1]*nodes_[3]->XYZ[0];
	Ja[1] = GN[2]*nodes_[0]->XYZ[0] - GN[3]*nodes_[1]->XYZ[0] + GN[3]*nodes_[2]->XYZ[0] - GN[2]*nodes_[3]->XYZ[0];
	Ja[2] = GN[0]*nodes_[0]->XYZ[1] - GN[0]*nodes_[1]->XYZ[1] + GN[1]*nodes_[2]->XYZ[1] - GN[1]*nodes_[3]->XYZ[1];
	Ja[3] = GN[2]*nodes_[0]->XYZ[1] - GN[3]*nodes_[1]->XYZ[1] + GN[3]*nodes_[2]->XYZ[1] - GN[2]*nodes_[3]->XYZ[1];
	
	double JaDet;
	JaDet = Ja[0]*Ja[3]-Ja[2]*Ja[1];
//inversion
	double JaInv[4];
	JaInv[0] = Ja[3]/JaDet;
	JaInv[1] = -Ja[1]/JaDet;
	JaInv[2] = -Ja[2]/JaDet;
	JaInv[3] = Ja[0]/JaDet;

//Be [1x 1y 2x 2y 3x 3y 4x 4y JaDet]
	Be[0] = JaInv[0]*GN[0] + JaInv[2]*GN[2];
	Be[1] = JaInv[1]*GN[0] + JaInv[3]*GN[2];
	Be[2] = -JaInv[0]*GN[0] - JaInv[2]*GN[3];
	Be[3] = -JaInv[1]*GN[0] - JaInv[3]*GN[3];
	Be[4] = JaInv[0]*GN[1] + JaInv[2]*GN[3];
	Be[5] = JaInv[1]*GN[1] + JaInv[3]*GN[3];
	Be[6] = -JaInv[0]*GN[1] - JaInv[2]*GN[2];
	Be[7] = -JaInv[1]*GN[1] - JaInv[3]*GN[2];
	Be[8] = JaDet;

}


void CQ4::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	
//constitutive matrix
	double CM[3];
	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);
	CM[0] = material_->E / (1- material_->Nu * material_->Nu);
	CM[1] = CM[0] * material_->Nu;
	CM[2] = (CM[0] - CM[1])/2;
	
	double Be[9];
	double Ga=0.57735;
	
	CalculateBe(Be,-Ga,-Ga);
	Matrix[0] = Matrix[0] + ( Be[0]*CM[0]*Be[0] + Be[1]*CM[2]*Be[1] )*Be[8]; //b0*c0*conj(b0) + b1*c2*conj(b1)
	Matrix[1] = Matrix[1] + ( Be[0]*CM[2]*Be[0] + Be[1]*CM[0]*Be[1] )*Be[8]; //b0*c2*conj(b0) + b1*c0*conj(b1)
	Matrix[2] = Matrix[2] + ( Be[1]*CM[1]*Be[0] + Be[0]*CM[2]*Be[1] )*Be[8]; //b1*c1*conj(b0) + b0*c2*conj(b1)
	Matrix[3] = Matrix[3] + ( Be[2]*CM[0]*Be[2] + Be[3]*CM[2]*Be[3] )*Be[8]; //b2*c0*conj(b2) + b3*c2*conj(b3)
	Matrix[4] = Matrix[4] + ( Be[2]*CM[1]*Be[1] + Be[3]*CM[2]*Be[0] )*Be[8]; //b2*c1*conj(b1) + b3*c2*conj(b0)
	Matrix[5] = Matrix[5] + ( Be[2]*CM[0]*Be[0] + Be[3]*CM[2]*Be[1] )*Be[8]; //b2*c0*conj(b0) + b3*c2*conj(b1)
	Matrix[6] = Matrix[6] + ( Be[2]*CM[2]*Be[2] + Be[3]*CM[0]*Be[3] )*Be[8]; //b2*c2*conj(b2) + b3*c0*conj(b3)
	Matrix[7] = Matrix[7] + ( Be[3]*CM[1]*Be[2] + Be[2]*CM[2]*Be[3] )*Be[8]; //b3*c1*conj(b2) + b2*c2*conj(b3)
	Matrix[8] = Matrix[8] + ( Be[2]*CM[2]*Be[0] + Be[3]*CM[0]*Be[1] )*Be[8]; //b2*c2*conj(b0) + b3*c0*conj(b1)
	Matrix[9] = Matrix[9] + ( Be[3]*CM[1]*Be[0] + Be[2]*CM[2]*Be[1] )*Be[8]; //b3*c1*conj(b0) + b2*c2*conj(b1)
	Matrix[10] = Matrix[10] + ( Be[4]*CM[0]*Be[4] + Be[5]*CM[2]*Be[5] )*Be[8]; //b4*c0*conj(b4) + b5*c2*conj(b5)
	Matrix[11] = Matrix[11] + ( Be[4]*CM[1]*Be[3] + Be[5]*CM[2]*Be[2] )*Be[8]; //b4*c1*conj(b3) + b5*c2*conj(b2)
	Matrix[12] = Matrix[12] + ( Be[4]*CM[0]*Be[2] + Be[5]*CM[2]*Be[3] )*Be[8]; //b4*c0*conj(b2) + b5*c2*conj(b3)
	Matrix[13] = Matrix[13] + ( Be[4]*CM[1]*Be[1] + Be[5]*CM[2]*Be[0] )*Be[8]; //b4*c1*conj(b1) + b5*c2*conj(b0)
	Matrix[14] = Matrix[14] + ( Be[4]*CM[0]*Be[0] + Be[5]*CM[2]*Be[1] )*Be[8]; //b4*c0*conj(b0) + b5*c2*conj(b1)
	Matrix[15] = Matrix[15] + ( Be[4]*CM[2]*Be[4] + Be[5]*CM[0]*Be[5] )*Be[8]; //b4*c2*conj(b4) + b5*c0*conj(b5)
	Matrix[16] = Matrix[16] + ( Be[5]*CM[1]*Be[4] + Be[4]*CM[2]*Be[5] )*Be[8]; //b5*c1*conj(b4) + b4*c2*conj(b5)
	Matrix[17] = Matrix[17] + ( Be[4]*CM[2]*Be[2] + Be[5]*CM[0]*Be[3] )*Be[8]; //b4*c2*conj(b2) + b5*c0*conj(b3)
	Matrix[18] = Matrix[18] + ( Be[5]*CM[1]*Be[2] + Be[4]*CM[2]*Be[3] )*Be[8]; //b5*c1*conj(b2) + b4*c2*conj(b3)
	Matrix[19] = Matrix[19] + ( Be[4]*CM[2]*Be[0] + Be[5]*CM[0]*Be[1] )*Be[8]; //b4*c2*conj(b0) + b5*c0*conj(b1)
	Matrix[20] = Matrix[20] + ( Be[5]*CM[1]*Be[0] + Be[4]*CM[2]*Be[1] )*Be[8]; //b5*c1*conj(b0) + b4*c2*conj(b1)
	Matrix[21] = Matrix[21] + ( Be[6]*CM[0]*Be[6] + Be[7]*CM[2]*Be[7] )*Be[8]; //b6*c0*conj(b6) + b7*c2*conj(b7)
	Matrix[22] = Matrix[22] + ( Be[6]*CM[1]*Be[5] + Be[7]*CM[2]*Be[4] )*Be[8]; //b6*c1*conj(b5) + b7*c2*conj(b4)
	Matrix[23] = Matrix[23] + ( Be[6]*CM[0]*Be[4] + Be[7]*CM[2]*Be[5] )*Be[8]; //b6*c0*conj(b4) + b7*c2*conj(b5)
	Matrix[24] = Matrix[24] + ( Be[6]*CM[1]*Be[3] + Be[7]*CM[2]*Be[2] )*Be[8]; //b6*c1*conj(b3) + b7*c2*conj(b2)
	Matrix[25] = Matrix[25] + ( Be[6]*CM[0]*Be[2] + Be[7]*CM[2]*Be[3] )*Be[8]; //b6*c0*conj(b2) + b7*c2*conj(b3)
	Matrix[26] = Matrix[26] + ( Be[6]*CM[1]*Be[1] + Be[7]*CM[2]*Be[0] )*Be[8]; //b6*c1*conj(b1) + b7*c2*conj(b0)
	Matrix[27] = Matrix[27] + ( Be[6]*CM[0]*Be[0] + Be[7]*CM[2]*Be[1] )*Be[8]; //b6*c0*conj(b0) + b7*c2*conj(b1)
	Matrix[28] = Matrix[28] + ( Be[6]*CM[2]*Be[6] + Be[7]*CM[0]*Be[7] )*Be[8]; //b6*c2*conj(b6) + b7*c0*conj(b7)
	Matrix[29] = Matrix[29] + ( Be[7]*CM[1]*Be[6] + Be[6]*CM[2]*Be[7] )*Be[8]; //b7*c1*conj(b6) + b6*c2*conj(b7)
	Matrix[30] = Matrix[30] + ( Be[6]*CM[2]*Be[4] + Be[7]*CM[0]*Be[5] )*Be[8]; //b6*c2*conj(b4) + b7*c0*conj(b5)
	Matrix[31] = Matrix[31] + ( Be[7]*CM[1]*Be[4] + Be[6]*CM[2]*Be[5] )*Be[8]; //b7*c1*conj(b4) + b6*c2*conj(b5)
	Matrix[32] = Matrix[32] + ( Be[6]*CM[2]*Be[2] + Be[7]*CM[0]*Be[3] )*Be[8]; //b6*c2*conj(b2) + b7*c0*conj(b3)
	Matrix[33] = Matrix[33] + ( Be[7]*CM[1]*Be[2] + Be[6]*CM[2]*Be[3] )*Be[8]; //b7*c1*conj(b2) + b6*c2*conj(b3)
	Matrix[34] = Matrix[34] + ( Be[6]*CM[2]*Be[0] + Be[7]*CM[0]*Be[1] )*Be[8]; //b6*c2*conj(b0) + b7*c0*conj(b1)
	Matrix[35] = Matrix[35] + ( Be[7]*CM[1]*Be[0] + Be[6]*CM[2]*Be[1] )*Be[8]; //b7*c1*conj(b0) + b6*c2*conj(b1)
	
	CalculateBe(Be,Ga,-Ga);
	Matrix[0] = Matrix[0] + ( Be[0]*CM[0]*Be[0] + Be[1]*CM[2]*Be[1] )*Be[8]; //b0*c0*conj(b0) + b1*c2*conj(b1)
	Matrix[1] = Matrix[1] + ( Be[0]*CM[2]*Be[0] + Be[1]*CM[0]*Be[1] )*Be[8]; //b0*c2*conj(b0) + b1*c0*conj(b1)
	Matrix[2] = Matrix[2] + ( Be[1]*CM[1]*Be[0] + Be[0]*CM[2]*Be[1] )*Be[8]; //b1*c1*conj(b0) + b0*c2*conj(b1)
	Matrix[3] = Matrix[3] + ( Be[2]*CM[0]*Be[2] + Be[3]*CM[2]*Be[3] )*Be[8]; //b2*c0*conj(b2) + b3*c2*conj(b3)
	Matrix[4] = Matrix[4] + ( Be[2]*CM[1]*Be[1] + Be[3]*CM[2]*Be[0] )*Be[8]; //b2*c1*conj(b1) + b3*c2*conj(b0)
	Matrix[5] = Matrix[5] + ( Be[2]*CM[0]*Be[0] + Be[3]*CM[2]*Be[1] )*Be[8]; //b2*c0*conj(b0) + b3*c2*conj(b1)
	Matrix[6] = Matrix[6] + ( Be[2]*CM[2]*Be[2] + Be[3]*CM[0]*Be[3] )*Be[8]; //b2*c2*conj(b2) + b3*c0*conj(b3)
	Matrix[7] = Matrix[7] + ( Be[3]*CM[1]*Be[2] + Be[2]*CM[2]*Be[3] )*Be[8]; //b3*c1*conj(b2) + b2*c2*conj(b3)
	Matrix[8] = Matrix[8] + ( Be[2]*CM[2]*Be[0] + Be[3]*CM[0]*Be[1] )*Be[8]; //b2*c2*conj(b0) + b3*c0*conj(b1)
	Matrix[9] = Matrix[9] + ( Be[3]*CM[1]*Be[0] + Be[2]*CM[2]*Be[1] )*Be[8]; //b3*c1*conj(b0) + b2*c2*conj(b1)
	Matrix[10] = Matrix[10] + ( Be[4]*CM[0]*Be[4] + Be[5]*CM[2]*Be[5] )*Be[8]; //b4*c0*conj(b4) + b5*c2*conj(b5)
	Matrix[11] = Matrix[11] + ( Be[4]*CM[1]*Be[3] + Be[5]*CM[2]*Be[2] )*Be[8]; //b4*c1*conj(b3) + b5*c2*conj(b2)
	Matrix[12] = Matrix[12] + ( Be[4]*CM[0]*Be[2] + Be[5]*CM[2]*Be[3] )*Be[8]; //b4*c0*conj(b2) + b5*c2*conj(b3)
	Matrix[13] = Matrix[13] + ( Be[4]*CM[1]*Be[1] + Be[5]*CM[2]*Be[0] )*Be[8]; //b4*c1*conj(b1) + b5*c2*conj(b0)
	Matrix[14] = Matrix[14] + ( Be[4]*CM[0]*Be[0] + Be[5]*CM[2]*Be[1] )*Be[8]; //b4*c0*conj(b0) + b5*c2*conj(b1)
	Matrix[15] = Matrix[15] + ( Be[4]*CM[2]*Be[4] + Be[5]*CM[0]*Be[5] )*Be[8]; //b4*c2*conj(b4) + b5*c0*conj(b5)
	Matrix[16] = Matrix[16] + ( Be[5]*CM[1]*Be[4] + Be[4]*CM[2]*Be[5] )*Be[8]; //b5*c1*conj(b4) + b4*c2*conj(b5)
	Matrix[17] = Matrix[17] + ( Be[4]*CM[2]*Be[2] + Be[5]*CM[0]*Be[3] )*Be[8]; //b4*c2*conj(b2) + b5*c0*conj(b3)
	Matrix[18] = Matrix[18] + ( Be[5]*CM[1]*Be[2] + Be[4]*CM[2]*Be[3] )*Be[8]; //b5*c1*conj(b2) + b4*c2*conj(b3)
	Matrix[19] = Matrix[19] + ( Be[4]*CM[2]*Be[0] + Be[5]*CM[0]*Be[1] )*Be[8]; //b4*c2*conj(b0) + b5*c0*conj(b1)
	Matrix[20] = Matrix[20] + ( Be[5]*CM[1]*Be[0] + Be[4]*CM[2]*Be[1] )*Be[8]; //b5*c1*conj(b0) + b4*c2*conj(b1)
	Matrix[21] = Matrix[21] + ( Be[6]*CM[0]*Be[6] + Be[7]*CM[2]*Be[7] )*Be[8]; //b6*c0*conj(b6) + b7*c2*conj(b7)
	Matrix[22] = Matrix[22] + ( Be[6]*CM[1]*Be[5] + Be[7]*CM[2]*Be[4] )*Be[8]; //b6*c1*conj(b5) + b7*c2*conj(b4)
	Matrix[23] = Matrix[23] + ( Be[6]*CM[0]*Be[4] + Be[7]*CM[2]*Be[5] )*Be[8]; //b6*c0*conj(b4) + b7*c2*conj(b5)
	Matrix[24] = Matrix[24] + ( Be[6]*CM[1]*Be[3] + Be[7]*CM[2]*Be[2] )*Be[8]; //b6*c1*conj(b3) + b7*c2*conj(b2)
	Matrix[25] = Matrix[25] + ( Be[6]*CM[0]*Be[2] + Be[7]*CM[2]*Be[3] )*Be[8]; //b6*c0*conj(b2) + b7*c2*conj(b3)
	Matrix[26] = Matrix[26] + ( Be[6]*CM[1]*Be[1] + Be[7]*CM[2]*Be[0] )*Be[8]; //b6*c1*conj(b1) + b7*c2*conj(b0)
	Matrix[27] = Matrix[27] + ( Be[6]*CM[0]*Be[0] + Be[7]*CM[2]*Be[1] )*Be[8]; //b6*c0*conj(b0) + b7*c2*conj(b1)
	Matrix[28] = Matrix[28] + ( Be[6]*CM[2]*Be[6] + Be[7]*CM[0]*Be[7] )*Be[8]; //b6*c2*conj(b6) + b7*c0*conj(b7)
	Matrix[29] = Matrix[29] + ( Be[7]*CM[1]*Be[6] + Be[6]*CM[2]*Be[7] )*Be[8]; //b7*c1*conj(b6) + b6*c2*conj(b7)
	Matrix[30] = Matrix[30] + ( Be[6]*CM[2]*Be[4] + Be[7]*CM[0]*Be[5] )*Be[8]; //b6*c2*conj(b4) + b7*c0*conj(b5)
	Matrix[31] = Matrix[31] + ( Be[7]*CM[1]*Be[4] + Be[6]*CM[2]*Be[5] )*Be[8]; //b7*c1*conj(b4) + b6*c2*conj(b5)
	Matrix[32] = Matrix[32] + ( Be[6]*CM[2]*Be[2] + Be[7]*CM[0]*Be[3] )*Be[8]; //b6*c2*conj(b2) + b7*c0*conj(b3)
	Matrix[33] = Matrix[33] + ( Be[7]*CM[1]*Be[2] + Be[6]*CM[2]*Be[3] )*Be[8]; //b7*c1*conj(b2) + b6*c2*conj(b3)
	Matrix[34] = Matrix[34] + ( Be[6]*CM[2]*Be[0] + Be[7]*CM[0]*Be[1] )*Be[8]; //b6*c2*conj(b0) + b7*c0*conj(b1)
	Matrix[35] = Matrix[35] + ( Be[7]*CM[1]*Be[0] + Be[6]*CM[2]*Be[1] )*Be[8]; //b7*c1*conj(b0) + b6*c2*conj(b1)

	
	CalculateBe(Be,-Ga,Ga);
	Matrix[0] = Matrix[0] + ( Be[0]*CM[0]*Be[0] + Be[1]*CM[2]*Be[1] )*Be[8]; //b0*c0*conj(b0) + b1*c2*conj(b1)
	Matrix[1] = Matrix[1] + ( Be[0]*CM[2]*Be[0] + Be[1]*CM[0]*Be[1] )*Be[8]; //b0*c2*conj(b0) + b1*c0*conj(b1)
	Matrix[2] = Matrix[2] + ( Be[1]*CM[1]*Be[0] + Be[0]*CM[2]*Be[1] )*Be[8]; //b1*c1*conj(b0) + b0*c2*conj(b1)
	Matrix[3] = Matrix[3] + ( Be[2]*CM[0]*Be[2] + Be[3]*CM[2]*Be[3] )*Be[8]; //b2*c0*conj(b2) + b3*c2*conj(b3)
	Matrix[4] = Matrix[4] + ( Be[2]*CM[1]*Be[1] + Be[3]*CM[2]*Be[0] )*Be[8]; //b2*c1*conj(b1) + b3*c2*conj(b0)
	Matrix[5] = Matrix[5] + ( Be[2]*CM[0]*Be[0] + Be[3]*CM[2]*Be[1] )*Be[8]; //b2*c0*conj(b0) + b3*c2*conj(b1)
	Matrix[6] = Matrix[6] + ( Be[2]*CM[2]*Be[2] + Be[3]*CM[0]*Be[3] )*Be[8]; //b2*c2*conj(b2) + b3*c0*conj(b3)
	Matrix[7] = Matrix[7] + ( Be[3]*CM[1]*Be[2] + Be[2]*CM[2]*Be[3] )*Be[8]; //b3*c1*conj(b2) + b2*c2*conj(b3)
	Matrix[8] = Matrix[8] + ( Be[2]*CM[2]*Be[0] + Be[3]*CM[0]*Be[1] )*Be[8]; //b2*c2*conj(b0) + b3*c0*conj(b1)
	Matrix[9] = Matrix[9] + ( Be[3]*CM[1]*Be[0] + Be[2]*CM[2]*Be[1] )*Be[8]; //b3*c1*conj(b0) + b2*c2*conj(b1)
	Matrix[10] = Matrix[10] + ( Be[4]*CM[0]*Be[4] + Be[5]*CM[2]*Be[5] )*Be[8]; //b4*c0*conj(b4) + b5*c2*conj(b5)
	Matrix[11] = Matrix[11] + ( Be[4]*CM[1]*Be[3] + Be[5]*CM[2]*Be[2] )*Be[8]; //b4*c1*conj(b3) + b5*c2*conj(b2)
	Matrix[12] = Matrix[12] + ( Be[4]*CM[0]*Be[2] + Be[5]*CM[2]*Be[3] )*Be[8]; //b4*c0*conj(b2) + b5*c2*conj(b3)
	Matrix[13] = Matrix[13] + ( Be[4]*CM[1]*Be[1] + Be[5]*CM[2]*Be[0] )*Be[8]; //b4*c1*conj(b1) + b5*c2*conj(b0)
	Matrix[14] = Matrix[14] + ( Be[4]*CM[0]*Be[0] + Be[5]*CM[2]*Be[1] )*Be[8]; //b4*c0*conj(b0) + b5*c2*conj(b1)
	Matrix[15] = Matrix[15] + ( Be[4]*CM[2]*Be[4] + Be[5]*CM[0]*Be[5] )*Be[8]; //b4*c2*conj(b4) + b5*c0*conj(b5)
	Matrix[16] = Matrix[16] + ( Be[5]*CM[1]*Be[4] + Be[4]*CM[2]*Be[5] )*Be[8]; //b5*c1*conj(b4) + b4*c2*conj(b5)
	Matrix[17] = Matrix[17] + ( Be[4]*CM[2]*Be[2] + Be[5]*CM[0]*Be[3] )*Be[8]; //b4*c2*conj(b2) + b5*c0*conj(b3)
	Matrix[18] = Matrix[18] + ( Be[5]*CM[1]*Be[2] + Be[4]*CM[2]*Be[3] )*Be[8]; //b5*c1*conj(b2) + b4*c2*conj(b3)
	Matrix[19] = Matrix[19] + ( Be[4]*CM[2]*Be[0] + Be[5]*CM[0]*Be[1] )*Be[8]; //b4*c2*conj(b0) + b5*c0*conj(b1)
	Matrix[20] = Matrix[20] + ( Be[5]*CM[1]*Be[0] + Be[4]*CM[2]*Be[1] )*Be[8]; //b5*c1*conj(b0) + b4*c2*conj(b1)
	Matrix[21] = Matrix[21] + ( Be[6]*CM[0]*Be[6] + Be[7]*CM[2]*Be[7] )*Be[8]; //b6*c0*conj(b6) + b7*c2*conj(b7)
	Matrix[22] = Matrix[22] + ( Be[6]*CM[1]*Be[5] + Be[7]*CM[2]*Be[4] )*Be[8]; //b6*c1*conj(b5) + b7*c2*conj(b4)
	Matrix[23] = Matrix[23] + ( Be[6]*CM[0]*Be[4] + Be[7]*CM[2]*Be[5] )*Be[8]; //b6*c0*conj(b4) + b7*c2*conj(b5)
	Matrix[24] = Matrix[24] + ( Be[6]*CM[1]*Be[3] + Be[7]*CM[2]*Be[2] )*Be[8]; //b6*c1*conj(b3) + b7*c2*conj(b2)
	Matrix[25] = Matrix[25] + ( Be[6]*CM[0]*Be[2] + Be[7]*CM[2]*Be[3] )*Be[8]; //b6*c0*conj(b2) + b7*c2*conj(b3)
	Matrix[26] = Matrix[26] + ( Be[6]*CM[1]*Be[1] + Be[7]*CM[2]*Be[0] )*Be[8]; //b6*c1*conj(b1) + b7*c2*conj(b0)
	Matrix[27] = Matrix[27] + ( Be[6]*CM[0]*Be[0] + Be[7]*CM[2]*Be[1] )*Be[8]; //b6*c0*conj(b0) + b7*c2*conj(b1)
	Matrix[28] = Matrix[28] + ( Be[6]*CM[2]*Be[6] + Be[7]*CM[0]*Be[7] )*Be[8]; //b6*c2*conj(b6) + b7*c0*conj(b7)
	Matrix[29] = Matrix[29] + ( Be[7]*CM[1]*Be[6] + Be[6]*CM[2]*Be[7] )*Be[8]; //b7*c1*conj(b6) + b6*c2*conj(b7)
	Matrix[30] = Matrix[30] + ( Be[6]*CM[2]*Be[4] + Be[7]*CM[0]*Be[5] )*Be[8]; //b6*c2*conj(b4) + b7*c0*conj(b5)
	Matrix[31] = Matrix[31] + ( Be[7]*CM[1]*Be[4] + Be[6]*CM[2]*Be[5] )*Be[8]; //b7*c1*conj(b4) + b6*c2*conj(b5)
	Matrix[32] = Matrix[32] + ( Be[6]*CM[2]*Be[2] + Be[7]*CM[0]*Be[3] )*Be[8]; //b6*c2*conj(b2) + b7*c0*conj(b3)
	Matrix[33] = Matrix[33] + ( Be[7]*CM[1]*Be[2] + Be[6]*CM[2]*Be[3] )*Be[8]; //b7*c1*conj(b2) + b6*c2*conj(b3)
	Matrix[34] = Matrix[34] + ( Be[6]*CM[2]*Be[0] + Be[7]*CM[0]*Be[1] )*Be[8]; //b6*c2*conj(b0) + b7*c0*conj(b1)
	Matrix[35] = Matrix[35] + ( Be[7]*CM[1]*Be[0] + Be[6]*CM[2]*Be[1] )*Be[8]; //b7*c1*conj(b0) + b6*c2*conj(b1)
	

	CalculateBe(Be,Ga,Ga);
	Matrix[0] = Matrix[0] + ( Be[0]*CM[0]*Be[0] + Be[1]*CM[2]*Be[1] )*Be[8]; //b0*c0*conj(b0) + b1*c2*conj(b1)
	Matrix[1] = Matrix[1] + ( Be[0]*CM[2]*Be[0] + Be[1]*CM[0]*Be[1] )*Be[8]; //b0*c2*conj(b0) + b1*c0*conj(b1)
	Matrix[2] = Matrix[2] + ( Be[1]*CM[1]*Be[0] + Be[0]*CM[2]*Be[1] )*Be[8]; //b1*c1*conj(b0) + b0*c2*conj(b1)
	Matrix[3] = Matrix[3] + ( Be[2]*CM[0]*Be[2] + Be[3]*CM[2]*Be[3] )*Be[8]; //b2*c0*conj(b2) + b3*c2*conj(b3)
	Matrix[4] = Matrix[4] + ( Be[2]*CM[1]*Be[1] + Be[3]*CM[2]*Be[0] )*Be[8]; //b2*c1*conj(b1) + b3*c2*conj(b0)
	Matrix[5] = Matrix[5] + ( Be[2]*CM[0]*Be[0] + Be[3]*CM[2]*Be[1] )*Be[8]; //b2*c0*conj(b0) + b3*c2*conj(b1)
	Matrix[6] = Matrix[6] + ( Be[2]*CM[2]*Be[2] + Be[3]*CM[0]*Be[3] )*Be[8]; //b2*c2*conj(b2) + b3*c0*conj(b3)
	Matrix[7] = Matrix[7] + ( Be[3]*CM[1]*Be[2] + Be[2]*CM[2]*Be[3] )*Be[8]; //b3*c1*conj(b2) + b2*c2*conj(b3)
	Matrix[8] = Matrix[8] + ( Be[2]*CM[2]*Be[0] + Be[3]*CM[0]*Be[1] )*Be[8]; //b2*c2*conj(b0) + b3*c0*conj(b1)
	Matrix[9] = Matrix[9] + ( Be[3]*CM[1]*Be[0] + Be[2]*CM[2]*Be[1] )*Be[8]; //b3*c1*conj(b0) + b2*c2*conj(b1)
	Matrix[10] = Matrix[10] + ( Be[4]*CM[0]*Be[4] + Be[5]*CM[2]*Be[5] )*Be[8]; //b4*c0*conj(b4) + b5*c2*conj(b5)
	Matrix[11] = Matrix[11] + ( Be[4]*CM[1]*Be[3] + Be[5]*CM[2]*Be[2] )*Be[8]; //b4*c1*conj(b3) + b5*c2*conj(b2)
	Matrix[12] = Matrix[12] + ( Be[4]*CM[0]*Be[2] + Be[5]*CM[2]*Be[3] )*Be[8]; //b4*c0*conj(b2) + b5*c2*conj(b3)
	Matrix[13] = Matrix[13] + ( Be[4]*CM[1]*Be[1] + Be[5]*CM[2]*Be[0] )*Be[8]; //b4*c1*conj(b1) + b5*c2*conj(b0)
	Matrix[14] = Matrix[14] + ( Be[4]*CM[0]*Be[0] + Be[5]*CM[2]*Be[1] )*Be[8]; //b4*c0*conj(b0) + b5*c2*conj(b1)
	Matrix[15] = Matrix[15] + ( Be[4]*CM[2]*Be[4] + Be[5]*CM[0]*Be[5] )*Be[8]; //b4*c2*conj(b4) + b5*c0*conj(b5)
	Matrix[16] = Matrix[16] + ( Be[5]*CM[1]*Be[4] + Be[4]*CM[2]*Be[5] )*Be[8]; //b5*c1*conj(b4) + b4*c2*conj(b5)
	Matrix[17] = Matrix[17] + ( Be[4]*CM[2]*Be[2] + Be[5]*CM[0]*Be[3] )*Be[8]; //b4*c2*conj(b2) + b5*c0*conj(b3)
	Matrix[18] = Matrix[18] + ( Be[5]*CM[1]*Be[2] + Be[4]*CM[2]*Be[3] )*Be[8]; //b5*c1*conj(b2) + b4*c2*conj(b3)
	Matrix[19] = Matrix[19] + ( Be[4]*CM[2]*Be[0] + Be[5]*CM[0]*Be[1] )*Be[8]; //b4*c2*conj(b0) + b5*c0*conj(b1)
	Matrix[20] = Matrix[20] + ( Be[5]*CM[1]*Be[0] + Be[4]*CM[2]*Be[1] )*Be[8]; //b5*c1*conj(b0) + b4*c2*conj(b1)
	Matrix[21] = Matrix[21] + ( Be[6]*CM[0]*Be[6] + Be[7]*CM[2]*Be[7] )*Be[8]; //b6*c0*conj(b6) + b7*c2*conj(b7)
	Matrix[22] = Matrix[22] + ( Be[6]*CM[1]*Be[5] + Be[7]*CM[2]*Be[4] )*Be[8]; //b6*c1*conj(b5) + b7*c2*conj(b4)
	Matrix[23] = Matrix[23] + ( Be[6]*CM[0]*Be[4] + Be[7]*CM[2]*Be[5] )*Be[8]; //b6*c0*conj(b4) + b7*c2*conj(b5)
	Matrix[24] = Matrix[24] + ( Be[6]*CM[1]*Be[3] + Be[7]*CM[2]*Be[2] )*Be[8]; //b6*c1*conj(b3) + b7*c2*conj(b2)
	Matrix[25] = Matrix[25] + ( Be[6]*CM[0]*Be[2] + Be[7]*CM[2]*Be[3] )*Be[8]; //b6*c0*conj(b2) + b7*c2*conj(b3)
	Matrix[26] = Matrix[26] + ( Be[6]*CM[1]*Be[1] + Be[7]*CM[2]*Be[0] )*Be[8]; //b6*c1*conj(b1) + b7*c2*conj(b0)
	Matrix[27] = Matrix[27] + ( Be[6]*CM[0]*Be[0] + Be[7]*CM[2]*Be[1] )*Be[8]; //b6*c0*conj(b0) + b7*c2*conj(b1)
	Matrix[28] = Matrix[28] + ( Be[6]*CM[2]*Be[6] + Be[7]*CM[0]*Be[7] )*Be[8]; //b6*c2*conj(b6) + b7*c0*conj(b7)
	Matrix[29] = Matrix[29] + ( Be[7]*CM[1]*Be[6] + Be[6]*CM[2]*Be[7] )*Be[8]; //b7*c1*conj(b6) + b6*c2*conj(b7)
	Matrix[30] = Matrix[30] + ( Be[6]*CM[2]*Be[4] + Be[7]*CM[0]*Be[5] )*Be[8]; //b6*c2*conj(b4) + b7*c0*conj(b5)
	Matrix[31] = Matrix[31] + ( Be[7]*CM[1]*Be[4] + Be[6]*CM[2]*Be[5] )*Be[8]; //b7*c1*conj(b4) + b6*c2*conj(b5)
	Matrix[32] = Matrix[32] + ( Be[6]*CM[2]*Be[2] + Be[7]*CM[0]*Be[3] )*Be[8]; //b6*c2*conj(b2) + b7*c0*conj(b3)
	Matrix[33] = Matrix[33] + ( Be[7]*CM[1]*Be[2] + Be[6]*CM[2]*Be[3] )*Be[8]; //b7*c1*conj(b2) + b6*c2*conj(b3)
	Matrix[34] = Matrix[34] + ( Be[6]*CM[2]*Be[0] + Be[7]*CM[0]*Be[1] )*Be[8]; //b6*c2*conj(b0) + b7*c0*conj(b1)
	Matrix[35] = Matrix[35] + ( Be[7]*CM[1]*Be[0] + Be[6]*CM[2]*Be[1] )*Be[8]; //b7*c1*conj(b0) + b6*c2*conj(b1)

/*
	printf("ESM\n");
	for (int ii=0;ii<36;ii++)
	{
		printf("%lf\n",Matrix[ii]/1E7);
	}
*/	

}

//	Calculate element stress 
void CQ4::ElementStress(double* stress, double* Displacement)
{
	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);	// Pointer to material of the element

	double Be[9];
	double Ga=0.57735;
	double strain[3];
	double CM[3];

	CM[0] = material_->E / (1- material_->Nu * material_->Nu);
	CM[1] = CM[0] * material_->Nu;
	CM[2] = (CM[0] - CM[1])/2;
	

	CalculateBe(Be,-Ga,-Ga);
	strain[0] = Be[0]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) + Be[2]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[4]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) + Be[6]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0);
	strain[1] = Be[1]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) + Be[3]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[5]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) + Be[7]*Displacement[LocationMatrix_[7]-1]*(!LocationMatrix_[7]==0);
	strain[2] = Be[1]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) + Be[3]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[5]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) + Be[7]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0) + Be[0]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) + Be[2]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[4]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) + Be[6]*Displacement[LocationMatrix_[7]-1]*(double)(!LocationMatrix_[7]==0);
	
	stress[0]=CM[0]*strain[0] + CM[1]*strain[1];
	stress[1]=CM[1]*strain[0] + CM[0]*strain[1];
	stress[2]=CM[2]*strain[2];
	
	
	CalculateBe(Be,Ga,-Ga);
	strain[0] = Be[0]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) + Be[2]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[4]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) + Be[6]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0);
	strain[1] = Be[1]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) + Be[3]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[5]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) + Be[7]*Displacement[LocationMatrix_[7]-1]*(!LocationMatrix_[7]==0);
	strain[2] = Be[1]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) + Be[3]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[5]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) + Be[7]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0) + Be[0]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) + Be[2]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[4]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) + Be[6]*Displacement[LocationMatrix_[7]-1]*(double)(!LocationMatrix_[7]==0);
	
	stress[3]=CM[0]*strain[0] + CM[1]*strain[1];
	stress[4]=CM[1]*strain[0] + CM[0]*strain[1];
	stress[5]=CM[2]*strain[2];
	
	
	CalculateBe(Be,-Ga,Ga);
	strain[0] = Be[0]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) + Be[2]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[4]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) + Be[6]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0);
	strain[1] = Be[1]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) + Be[3]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[5]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) + Be[7]*Displacement[LocationMatrix_[7]-1]*(!LocationMatrix_[7]==0);
	strain[2] = Be[1]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) + Be[3]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[5]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) + Be[7]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0) + Be[0]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) + Be[2]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[4]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) + Be[6]*Displacement[LocationMatrix_[7]-1]*(double)(!LocationMatrix_[7]==0);
	
	stress[6]=CM[0]*strain[0] + CM[1]*strain[1];
	stress[7]=CM[1]*strain[0] + CM[0]*strain[1];
	stress[8]=CM[2]*strain[2];
	
	
	CalculateBe(Be,Ga,Ga);
	strain[0] = Be[0]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) + Be[2]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[4]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) + Be[6]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0);
	strain[1] = Be[1]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) + Be[3]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[5]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) + Be[7]*Displacement[LocationMatrix_[7]-1]*(!LocationMatrix_[7]==0);
	strain[2] = Be[1]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) + Be[3]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[5]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) + Be[7]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0) + Be[0]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) + Be[2]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[4]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) + Be[6]*Displacement[LocationMatrix_[7]-1]*(double)(!LocationMatrix_[7]==0);
	
	stress[9]=CM[0]*strain[0] + CM[1]*strain[1];
	stress[10]=CM[1]*strain[0] + CM[0]*strain[1];
	stress[11]=CM[2]*strain[2];
	
}

void CQ4::ElementCoord (double* coord)
{
	double Ga=0.57735;
	double eta,psi;
	
	eta = -Ga; psi = -Ga;
	coord[0]=( nodes_[0]->XYZ[0]*(1-eta)*(1-psi) + nodes_[1]->XYZ[0]*(1-eta)*(1+psi) + nodes_[2]->XYZ[0]*(1+eta)*(1+psi) + nodes_[3]->XYZ[0]*(1+eta)*(1-psi) )/4;
	coord[1]=( nodes_[0]->XYZ[1]*(1-eta)*(1-psi) + nodes_[1]->XYZ[1]*(1-eta)*(1+psi) + nodes_[2]->XYZ[1]*(1+eta)*(1+psi) + nodes_[3]->XYZ[1]*(1+eta)*(1-psi) )/4;

	eta = Ga; psi = -Ga;
	coord[2]=( nodes_[0]->XYZ[0]*(1-eta)*(1-psi) + nodes_[1]->XYZ[0]*(1-eta)*(1+psi) + nodes_[2]->XYZ[0]*(1+eta)*(1+psi) + nodes_[3]->XYZ[0]*(1+eta)*(1-psi) )/4;
	coord[3]=( nodes_[0]->XYZ[1]*(1-eta)*(1-psi) + nodes_[1]->XYZ[1]*(1-eta)*(1+psi) + nodes_[2]->XYZ[1]*(1+eta)*(1+psi) + nodes_[3]->XYZ[1]*(1+eta)*(1-psi) )/4;

	eta = -Ga; psi = Ga;
	coord[4]=( nodes_[0]->XYZ[0]*(1-eta)*(1-psi) + nodes_[1]->XYZ[0]*(1-eta)*(1+psi) + nodes_[2]->XYZ[0]*(1+eta)*(1+psi) + nodes_[3]->XYZ[0]*(1+eta)*(1-psi) )/4;
	coord[5]=( nodes_[0]->XYZ[1]*(1-eta)*(1-psi) + nodes_[1]->XYZ[1]*(1-eta)*(1+psi) + nodes_[2]->XYZ[1]*(1+eta)*(1+psi) + nodes_[3]->XYZ[1]*(1+eta)*(1-psi) )/4;

	eta = Ga; psi = Ga;
	coord[6]=( nodes_[0]->XYZ[0]*(1-eta)*(1-psi) + nodes_[1]->XYZ[0]*(1-eta)*(1+psi) + nodes_[2]->XYZ[0]*(1+eta)*(1+psi) + nodes_[3]->XYZ[0]*(1+eta)*(1-psi) )/4;
	coord[7]=( nodes_[0]->XYZ[1]*(1-eta)*(1-psi) + nodes_[1]->XYZ[1]*(1-eta)*(1+psi) + nodes_[2]->XYZ[1]*(1+eta)*(1+psi) + nodes_[3]->XYZ[1]*(1+eta)*(1-psi) )/4;
	
}
