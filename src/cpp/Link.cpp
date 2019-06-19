/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Link.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

//	Constructor
CLink::CLink()
{
	NEN_ = 2; // Each element has 2 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CLink::~CLink()
{
	delete[]nodes_;
	delete[]LocationMatrix_;
}

//	Read element data from stream Input
bool CLink::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N; // element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl
			<< "    Expected element : " << Ele + 1 << endl
			<< "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;   // Material property set number
	unsigned int N1, N2; // Left node number and right node number
	Input >> N1 >> N2 >> MSet;
	ElementMaterial_ = static_cast<CBeamMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CLink::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber
		<< setw(12) << ElementMaterial_->nset << endl;
}

//	Write element data to stream
void CLink::WritePlot(COutPlot& output, unsigned int Ele)
{
	output << 2 << setw(9) << nodes_[0]->NodeNumber - 1 << setw(9) << nodes_[1]->NodeNumber - 1 << endl;
}

void CLink::WritePlotPost(COutPlotPost& output, unsigned int Ele)
{
	output << 2 << setw(9) << nodes_[0]->NodeNumber - 1 << setw(9) << nodes_[1]->NodeNumber - 1 << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the
//  element
//	Caution:  Equation number is numbered from 1 !
void CLink::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
	{
		for (unsigned int D = 0; D < 6; D++)
		{
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
		}
	}
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node beam element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 78 elements
unsigned int CLink::SizeOfStiffnessMatrix() { return 78; }

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CLink::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	double k=1000000000000;
	Matrix[0]=k;
	Matrix[1]=k;
	Matrix[3]=k;
	Matrix[6]=k;
	Matrix[10]=k;
	Matrix[15]=k;
	Matrix[21]=k;
	Matrix[28]=k;
	Matrix[36]=k;
	Matrix[45]=k;
	Matrix[55]=k;
	Matrix[66]=k;
	Matrix[27]=-k;
	Matrix[34]=-k;
	Matrix[42]=-k;
	Matrix[51]=-k;
	Matrix[61]=-k;
	Matrix[72]=-k;
}
//	Calculate element stress
void CLink::ElementStress(double* stress, double* Displacement)
{}

void CLink::ElementMass(double* mass)
{
	clear(mass, 1);
}

void CLink::ElementPostInfo(double* beamstress, double* Displacement, double* prePositionBeam, double* postPositionBeam)
{}

void CLink::GravityCalculation(double* ptr_force)
{
	clear(ptr_force,12);
}