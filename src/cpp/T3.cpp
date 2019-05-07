/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "T3.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CT3::CT3()
{
	NEN_ = 3;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 6;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CT3::~CT3()
{
}

//	Read element data from stream Input
bool CT3::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N;	// element Number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl 
			 << "    Expected element : " << Ele + 1 << endl
			 << "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set Number
	unsigned int N1, N2, N3;	// 3 node Numbers

	Input >> N1 >> N2 >> N3 >> MSet;
    ElementMaterial_ = dynamic_cast<C4QMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	return true;
}

//	Write element data to stream
void CT3::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber << setw(9) << nodes_[1]->NodeNumber 
		<< setw(9) << nodes_[2]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//  Generate location matrix: the global equation Number that corresponding to each DOF of the element
//	Caution:  Equation Number is Numbered from 1 !
void CT3::GenerateLocationMatrix()
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
unsigned int CT3::SizeOfStiffnessMatrix() { return 21; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element


void CT3::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	
//constitutive matrix
	double CM[3];
	C4QMaterial* material_ = dynamic_cast<C4QMaterial*>(ElementMaterial_);
	CM[0] = material_->E / (1- material_->poisson * material_->poisson);
	CM[1] = CM[0] * material_->poisson;
	CM[2] = (CM[0] - CM[1])/2;
	
	double Be[6];
	Be[0] = nodes_[1]->XYZ[1] - nodes_[2]->XYZ[1];
	Be[1] = nodes_[2]->XYZ[0] - nodes_[1]->XYZ[0];
	Be[2] = nodes_[2]->XYZ[1] - nodes_[0]->XYZ[1];
	Be[3] = nodes_[0]->XYZ[0] - nodes_[2]->XYZ[0];
	Be[4] = nodes_[0]->XYZ[1] - nodes_[1]->XYZ[1];
	Be[5] = nodes_[1]->XYZ[0] - nodes_[0]->XYZ[0];
	
	double Area;
	Area = nodes_[0]->XYZ[0]*nodes_[1]->XYZ[1] + nodes_[1]->XYZ[0]*nodes_[2]->XYZ[1] + nodes_[2]->XYZ[0]*nodes_[0]->XYZ[1] - nodes_[0]->XYZ[0]*nodes_[2]->XYZ[1]- nodes_[1]->XYZ[0]*nodes_[0]->XYZ[1] - nodes_[2]->XYZ[0]*nodes_[1]->XYZ[1];
	Area = Area/2;
	
	Area = Area*material_->thick;
	
	Matrix[ 0] =Area*( Be[0]*CM[0]* Be[0]  + Be[1]*CM[2]* Be[1] );
	Matrix[ 1] =Area*( Be[0]*CM[2]* Be[0]  + Be[1]*CM[0]* Be[1] );
	Matrix[ 2] =Area*( Be[1]*CM[1]* Be[0]  + Be[0]*CM[2]* Be[1] );
	Matrix[ 3] =Area*( Be[2]*CM[0]* Be[2]  + Be[3]*CM[2]* Be[3] );
	Matrix[ 4] =Area*( Be[2]*CM[1]* Be[1]  + Be[3]*CM[2]* Be[0] );
	Matrix[ 5] =Area*( Be[2]*CM[0]* Be[0]  + Be[3]*CM[2]* Be[1] );
	Matrix[ 6] =Area*( Be[2]*CM[2]* Be[2]  + Be[3]*CM[0]* Be[3] );
	Matrix[ 7] =Area*( Be[3]*CM[1]* Be[2]  + Be[2]*CM[2]* Be[3] );
	Matrix[ 8] =Area*( Be[2]*CM[2]* Be[0]  + Be[3]*CM[0]* Be[1] );
	Matrix[ 9] =Area*( Be[3]*CM[1]* Be[0]  + Be[2]*CM[2]* Be[1] );
	Matrix[10] =Area*( Be[4]*CM[0]* Be[4]  + Be[5]*CM[2]* Be[5] );
	Matrix[11] =Area*( Be[4]*CM[1]* Be[3]  + Be[5]*CM[2]* Be[2] );
	Matrix[12] =Area*( Be[4]*CM[0]* Be[2]  + Be[5]*CM[2]* Be[3] );
	Matrix[13] =Area*( Be[4]*CM[1]* Be[1]  + Be[5]*CM[2]* Be[0] );
	Matrix[14] =Area*( Be[4]*CM[0]* Be[0]  + Be[5]*CM[2]* Be[1] );
	Matrix[15] =Area*( Be[4]*CM[2]* Be[4]  + Be[5]*CM[0]* Be[5] );
	Matrix[16] =Area*( Be[5]*CM[1]* Be[4]  + Be[4]*CM[2]* Be[5] );
	Matrix[17] =Area*( Be[4]*CM[2]* Be[2]  + Be[5]*CM[0]* Be[3] );
	Matrix[18] =Area*( Be[5]*CM[1]* Be[2]  + Be[4]*CM[2]* Be[3] );
	Matrix[19] =Area*( Be[4]*CM[2]* Be[0]  + Be[5]*CM[0]* Be[1] );
	Matrix[20] =Area*( Be[5]*CM[1]* Be[0]  + Be[4]*CM[2]* Be[1] );
	
}

//	Calculate element stress 
void CT3::ElementStress(double* stress, double* Displacement)
{
//constitutive matrix
	double CM[3];
	C4QMaterial* material_ = dynamic_cast<C4QMaterial*>(ElementMaterial_);
	CM[0] = material_->E / (1- material_->poisson * material_->poisson);
	CM[1] = CM[0] * material_->poisson;
	CM[2] = (CM[0] - CM[1])/2;
	
	double Be[6];
	Be[0] = nodes_[1]->XYZ[1] - nodes_[2]->XYZ[1];
	Be[1] = nodes_[2]->XYZ[0] - nodes_[1]->XYZ[0];
	Be[2] = nodes_[2]->XYZ[1] - nodes_[0]->XYZ[1];
	Be[3] = nodes_[0]->XYZ[0] - nodes_[2]->XYZ[0];
	Be[4] = nodes_[0]->XYZ[1] - nodes_[1]->XYZ[1];
	Be[5] = nodes_[1]->XYZ[0] - nodes_[0]->XYZ[0];
	
	double strain[3];
	strain[0] = Be[0]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) + Be[2]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[4]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0);
	strain[1] = Be[1]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) + Be[3]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[5]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0);
	strain[2] = Be[1]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) + Be[3]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[5]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) + Be[0]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) + Be[2]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[4]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) ;
	
	stress[0]=CM[0]*strain[0] + CM[1]*strain[1];
	stress[1]=CM[1]*strain[0] + CM[0]*strain[1];
	stress[2]=CM[2]*strain[2];
}

void CT3::ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* PostPositions)
{
	for (unsigned int i=0;i<NEN_;i++)
	{
		for (unsigned int j=0;j<2;j++)
		{
			PrePositions[2*i+j] = nodes_[i] ->XYZ[0];
			if (LocationMatrix_[2*i+j])
				PostPositions[2*i+j] = PrePositions[2*i+j] + Displacement[LocationMatrix_[2*i+j]-1];
			else
				PostPositions[2*i+j] = PrePositions[2*i+j];
		}	
	}
	ElementStress(stress, Displacement);
}