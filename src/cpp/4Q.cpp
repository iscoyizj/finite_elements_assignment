/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#define singular_tolerance 0.001

#include "4Q.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include "..\..\make\4Q.H"

using namespace std;

//	Constructor
C4Q::C4Q()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 8;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
C4Q::~C4Q()
{
}

//	Read element data from stream Input
bool C4Q::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2, N3, N4;	// Four node number orderly

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = dynamic_cast<CplaneMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void C4Q::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1
		<< setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber
		<< setw(9) << nodes_[2]->NodeNumber
		<< setw(9) << nodes_[3]->NodeNumber
		<< setw(12) << ElementMaterial_->nset << endl;
		
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void C4Q::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < 2; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 4Q element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has  elements
unsigned int C4Q::SizeOfStiffnessMatrix() { return 36; }

double C4Q::GravityofElement()
{
	return 0.0;
}

//  Calculate Shape function for 2D elasticity 4Q element
void C4Q::Nmat4Q(double eta, double psi, double* Nmat) 
{
	Nmat[0] = (1 - psi)*(1 - eta) / 4;
	Nmat[1] = (1 + psi)*(1 - eta) / 4;
	Nmat[2] = (1 + psi)*(1 + eta) / 4;
	Nmat[3] = (1 - psi)*(1 + eta) / 4;
}

// Calculate strain matrix for 4Q element
double C4Q::Bmat4Q(double eta, double psi, double* Bmat) 
{
	double GN[8];
	GN[0] = (eta - 1) / 4;
	GN[1] = (1 - eta) / 4;
	GN[2] = (eta + 1) / 4;
	GN[3] = (-eta - 1) / 4;
	GN[4] = (psi - 1) / 4;
	GN[5] = (-eta - 1) / 4;
	GN[6] = (psi + 1) / 4;
	GN[7] = (1 - psi) / 4;

	double J[4];
	J[0] = GN[0] * nodes_[0]->XYZ[0] + GN[1] * nodes_[1]->XYZ[0] + GN[2] * nodes_[2]->XYZ[0] + GN[3] * nodes_[3]->XYZ[0];
	J[1] = GN[4] * nodes_[0]->XYZ[0] + GN[5] * nodes_[1]->XYZ[0] + GN[6] * nodes_[2]->XYZ[0] + GN[7] * nodes_[3]->XYZ[0];
	J[2] = GN[0] * nodes_[0]->XYZ[1] + GN[1] * nodes_[1]->XYZ[1] + GN[2] * nodes_[2]->XYZ[1] + GN[3] * nodes_[3]->XYZ[1];
	J[3] = GN[4] * nodes_[0]->XYZ[1] + GN[5] * nodes_[1]->XYZ[1] + GN[6] * nodes_[2]->XYZ[1] + GN[7] * nodes_[3]->XYZ[1];

	double detJ;
	detJ = J[0] * J[3] - J[1] * J[2];

	if (abs(detJ) < singular_tolerance) 
	{
		cerr << "Jaccobian Matrix maybe singular" << endl;
	}
	
	double inv_J[4];
	inv_J[0] = J[3] / detJ;
	inv_J[1] = -J[1] / detJ;
	inv_J[2] = -J[2] / detJ;
	inv_J[3] = J[0] / detJ;

	Bmat[0] = J[0] * GN[0] + J[2] * GN[4];
	Bmat[1] = J[0] * GN[1] + J[2] * GN[5];
	Bmat[2] = J[0] * GN[2] + J[2] * GN[6];
	Bmat[3] = J[0] * GN[3] + J[2] * GN[7];
	Bmat[4] = J[1] * GN[0] + J[3] * GN[4];
	Bmat[5] = J[1] * GN[1] + J[3] * GN[5];
	Bmat[6] = J[1] * GN[2] + J[3] * GN[6];
	Bmat[7] = J[1] * GN[3] + J[3] * GN[7];
	return detJ;
}

//  Gauss Integral at gauss node (eta,psi)
//  Element stiffness matrix is sum of gauss integral
void C4Q::Integral_gauss_node(const double& eta, const double& psi, const double& weight, double* k_e)
{
	double * B;
	double * D;
	B = new double[8];
	D = new double[4];
	CplaneMaterial* material_ = dynamic_cast<CplaneMaterial*>(ElementMaterial_);	// Pointer to material of the element
	D = material_->Dmat;
	double d1(D[0]), v(material_->mu);
	double d33((1 - v) / 2.0);
	double DetJ = Bmat4Q(eta, psi, B);
	double w = weight * abs(DetJ)*d1;
	
	k_e[0] += w * (B[0] * B[0] + d33 * B[4] * B[4]);
	k_e[1] += w * (d33 * B[0] * B[4] + v * B[0] * B[4]);
	k_e[2] += w * (d33 * B[0] * B[0] + B[4] * B[4]);
	k_e[3] += w * (B[0] * B[1] + d33 * B[4] * B[5]);
	k_e[4] += w * (v * B[1] * B[4] + d33 * B[0] * B[5]);
	k_e[5] += w * (B[1] * B[1] + d33 * B[5] * B[5]);
	k_e[6] += w * (d33 * B[1] * B[4] + v * B[0] * B[5]);
	k_e[7] += w * (d33 * B[0] * B[1] + B[4] * B[5]);
	k_e[8] += w * (d33 * B[1] * B[5] + v * B[1] * B[5]);
	k_e[9] += w * (d33 * B[1] * B[1] + B[5] * B[5]);
	k_e[10] += w * (B[0] * B[2] + d33 * B[4] * B[6]);
	k_e[11] += w * (v * B[2] * B[4] + d33 * B[0] * B[6]);
	k_e[12] += w * (B[1] * B[2] + d33 * B[5] * B[6]);
	k_e[13] += w * (v * B[2] * B[5] + d33 * B[1] * B[6]);
	k_e[14] += w * (B[2] * B[2] + d33 * B[6] * B[6]);
	k_e[15] += w * (d33 * B[2] * B[4] + v * B[0] * B[6]);
	k_e[16] += w * (d33 * B[0] * B[2] + B[4] * B[6]);
	k_e[17] += w * (d33 * B[2] * B[5] + v * B[1] * B[6]);
	k_e[18] += w * (d33 * B[1] * B[2] + B[5] * B[6]);
	k_e[19] += w * (d33 * B[2] * B[6] + v * B[2] * B[6]);
	k_e[20] += w * (d33 * B[2] * B[2] + B[6] * B[6]);
	k_e[21] += w * (B[0] * B[3] + d33 * B[4] * B[7]);
	k_e[22] += w * (v * B[3] * B[4] + d33 * B[0] * B[7]);
	k_e[23] += w * (B[1] * B[3] + d33 * B[5] * B[7]);
	k_e[24] += w * (v * B[3] * B[5] + d33 * B[1] * B[7]);
	k_e[25] += w * (B[2] * B[3] + d33 * B[6] * B[7]);
	k_e[26] += w * (v * B[3] * B[6] + d33 * B[2] * B[7]);
	k_e[27] += w * (B[3] * B[3] + d33 * B[7] * B[7]);
	k_e[28] += w * (d33 * B[3] * B[4] + v * B[0] * B[7]);
	k_e[29] += w * (d33 * B[0] * B[3] + B[4] * B[7]);
	k_e[30] += w * (d33 * B[3] * B[5] + v * B[1] * B[7]);
	k_e[31] += w * (d33 * B[1] * B[3] + B[5] * B[7]);
	k_e[32] += w * (d33 * B[3] * B[6] + v * B[2] * B[7]);
	k_e[33] += w * (d33 * B[2] * B[3] + B[6] * B[7]);
	k_e[34] += w * (d33 * B[3] * B[7] + v * B[3] * B[7]);
	k_e[35] += w * (d33 * B[3] * B[3] + B[7] * B[7]);
}


//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void C4Q::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	//	Calculate element stiffness matrix
	double g_n = 1 / std::sqrt(3.0f);
	Integral_gauss_node(g_n, -g_n, 1, Matrix);
	Integral_gauss_node(g_n, g_n, 1, Matrix);
	Integral_gauss_node(-g_n, -g_n, 1, Matrix);
	Integral_gauss_node(-g_n, g_n, 1, Matrix);	
}


//	Calculate element stress 
void C4Q::ElementStress_gauss_node(double eta, double psi, double* stress, double* Displacement)
{
	CplaneMaterial* material_ = dynamic_cast<CplaneMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double * ele_d;
	ele_d = new double[8];
	clear(ele_d, 8);

	double * ele_strain;
	ele_strain = new double[3];
	clear(ele_strain, 3);

	double * Bmat;
	Bmat = new double[8];
	clear(Bmat, 8);
	Bmat4Q(eta, psi, Bmat);
	for (unsigned int i = 0; i < ND_; i++)
	{
		if (LocationMatrix_[i])
			ele_d[i] =  Displacement[LocationMatrix_[i] - 1];
	}
	ele_strain[0] = Bmat[0] * ele_d[0] + Bmat[1] * ele_d[2] + Bmat[2] * ele_d[4] + Bmat[3] * ele_d[6];
    ele_strain[1] = Bmat[4] * ele_d[1] + Bmat[5] * ele_d[3] + Bmat[6] * ele_d[5] + Bmat[7] * ele_d[7];
    ele_strain[2] = Bmat[4] * ele_d[0] + Bmat[0] * ele_d[1] + Bmat[5] * ele_d[2] + Bmat[1] * ele_d[3]+ Bmat[6] * ele_d[4] + Bmat[2] * ele_d[5] + Bmat[7] * ele_d[6] + Bmat[3] * ele_d[7];
	double *Dmat;
	Dmat = new double[4];
	clear(Dmat, 4);
	Dmat = material_->Dmat;
	stress[0] = Dmat[0] * ele_strain[0] + Dmat[2] * ele_strain[1];
	stress[1] = Dmat[2] * ele_strain[0] + Dmat[1] * ele_strain[1];
	stress[3] = Dmat[3] * ele_strain[3];
}
void C4Q::ElementStress(double* stress, double* Displacement) 
{
	clear(stress, 8);
	double g_n = 1 / std::sqrt(3.0f);
	ElementStress_gauss_node(g_n, -g_n, stress, Displacement);
	ElementStress_gauss_node(g_n, g_n, stress + 2, Displacement);
	ElementStress_gauss_node(-g_n, -g_n, stress + 4, Displacement);
	ElementStress_gauss_node(-g_n, -g_n, stress + 6, Displacement);
}
void C4Q::Gauss_node_coordinate(double eta,double psi, double* coordinte )
{
	double *Nmat;
	Nmat4Q(eta, psi, Nmat);
	coordinte[0] = Nmat[0] * nodes_[0]->XYZ[0] + Nmat[1] * nodes_[1]->XYZ[0] + Nmat[2] * nodes_[2]->XYZ[0] + Nmat[3] * nodes_[0]->XYZ[0];
	coordinte[1] = Nmat[0] * nodes_[0]->XYZ[0] + Nmat[1] * nodes_[1]->XYZ[1] + Nmat[2] * nodes_[2]->XYZ[1] + Nmat[3] * nodes_[3]->XYZ[1];

}
void C4Q::Stress_coordinate(double* coordinate)
{
	clear(coordinate, 8);
	double g_n = 1 / std::sqrt(3.0f);
	Gauss_node_coordinate(g_n, -g_n, coordinate);
	Gauss_node_coordinate(g_n, -g_n, coordinate+2);
	Gauss_node_coordinate(g_n, -g_n, coordinate+4);
	Gauss_node_coordinate(g_n, -g_n, coordinate+6);
}
