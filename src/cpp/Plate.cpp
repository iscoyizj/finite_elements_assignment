/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Plate.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CPlate::CPlate()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = NULL;
}

//	Desconstructor
CPlate::~CPlate()
{
}

//	Read element data from stream Input
bool CPlate::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2, N3, N4;	// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = &(dynamic_cast<CPlateMaterial*>(MaterialSets))[MSet - 1];
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CPlate::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber << setw(9)
		<< nodes_[2]->NodeNumber << setw(9) << nodes_[3]->NodeNumber
		<< setw(12)
		<< ElementMaterial_->nset << endl;
}

//	Write element data to stream
void CPlate::WritePlot(COutPlot& output, unsigned int Ele)
{
	output << 4 << setw(9) << nodes_[0]->NodeNumber-1
		<< setw(9) << nodes_[1]->NodeNumber-1 << setw(9)
		<< nodes_[2]->NodeNumber-1 << setw(9) << nodes_[3]->NodeNumber-1
		<< endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CPlate::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 2; D < 5; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];

}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 4 node plate element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 78 elements
unsigned int CPlate::SizeOfStiffnessMatrix() { return 78; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CPlate::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	CPlateMaterial* material = dynamic_cast<CPlateMaterial*>(ElementMaterial_);

	double E = material->E;
	double poisson = material->poisson;
	double t = material->thick;

	double a = abs(nodes_[1]->XYZ[0] - nodes_[0]->XYZ[0]) / 2;
	double b = abs(nodes_[3]->XYZ[1] - nodes_[0]->XYZ[1]) / 2;
	double D0 = E * t*t*t / 12 / (1 - poisson * poisson);

	double k1 = 21 - 6 * poisson + 30 * b*b / a / a + 30 * a*a / b / b;
	double k2 = 8 * b*b - 8 * poisson*b*b + 40 * a*a;
	double k3 = 8 * a*a - 8 * poisson*a*a + 40 * b*b;
	double k4 = 3 * b + 12 * poisson*b + 30 * a*a / b;
	double k5 = 3 * a + 12 * poisson*a + 30 * b*b / a;
	double k6 = 30 * poisson*a*b;
	double k7 = -21 + 6 * poisson - 30 * b*b / a / a + 15 * a*a / b / b;
	double k8 = -8 * b*b + 8 * poisson*b*b + 20 * a*a;
	double k9 = -2 * a*a + 2 * poisson*a*a + 20 * b*b;
	double k10 = -3 * b - 12 * poisson*b + 15 * a*a / b;
	double k11 = 3 * a - 3 * poisson*a + 30 * b*b / a;
	double k12 = 21 - 6 * poisson - 15 * b*b / a / a - 15 * a*a / b / b;
	double k13 = 2 * b*b - 2 * poisson*b*b + 10 * a*a;
	double k14 = 2 * a*a - 2 * poisson*a*a + 10 * b*b;
	double k15 = -3 * b + 3 * poisson*b + 15 * a*a / b;
	double k16 = -3 * a + 3 * poisson*a + 15 * b*b / a;
	double k17 = -21 + 6 * poisson + 15 * b*b / a / a - 30 * a*a / b / b;
	double k18 = -2 * b*b + 2 * poisson*b*b + 20 * a*a;
	double k19 = -8 * a*a + 8 * poisson*a*a + 20 * b*b;
	double k20 = 3 * b - 3 * poisson*b + 30 * a*a / b;
	double k21 = -3 * a - 12 * poisson*a + 15 * b*b / a;

	Matrix[0] = k1;
	Matrix[1] = k2;
	Matrix[2] = k4;
	Matrix[3] = k3;
	Matrix[4] = -k6;
	Matrix[5] = -k5;
	Matrix[6] = k1;
	Matrix[7] = k11;
	Matrix[8] = k10;
	Matrix[9] = k7;
	Matrix[10] = k2;
	Matrix[11] = k4;
	Matrix[12] = 0;
	Matrix[13] = k8;
	Matrix[14] = k10;
	Matrix[15] = k3;
	Matrix[16] = k6;
	Matrix[17] = k5;
	Matrix[18] = k9;
	Matrix[19] = 0;
	Matrix[20] = -k11;
	Matrix[21] = k1;
	Matrix[22] = k21;
	Matrix[23] = -k20;
	Matrix[24] = k17;
	Matrix[25] = k16;
	Matrix[26] = -k15;
	Matrix[27] = k12;
	Matrix[28] = k2;
	Matrix[29] = -k4;
	Matrix[30] = 0;
	Matrix[31] = k18;
	Matrix[32] = k20;
	Matrix[33] = 0;
	Matrix[34] = k13;
	Matrix[35] = k15;
	Matrix[36] = k3;
	Matrix[37] = -k6;
	Matrix[38] = k5;
	Matrix[39] = k19;
	Matrix[40] = 0;
	Matrix[41] = k21;
	Matrix[42] = k14;
	Matrix[43] = 0;
	Matrix[44] = -k16;
	Matrix[45] = k1;
	Matrix[46] = -k11;
	Matrix[47] = -k10;
	Matrix[48] = k7;
	Matrix[49] = -k16;
	Matrix[50] = -k15;
	Matrix[51] = k12;
	Matrix[52] = -k21;
	Matrix[53] = -k20;
	Matrix[54] = k17;
	Matrix[55] = k2;
	Matrix[56] = -k4;
	Matrix[57] = 0;
	Matrix[58] = k8;
	Matrix[59] = -k10;
	Matrix[60] = 0;
	Matrix[61] = k13;
	Matrix[62] = k15;
	Matrix[63] = 0;
	Matrix[64] = k18;
	Matrix[65] = k20;
	Matrix[66] = k3;
	Matrix[67] = k6;
	Matrix[68] = -k5;
	Matrix[69] = k9;
	Matrix[70] = 0;
	Matrix[71] = k11;
	Matrix[72] = k14;
	Matrix[73] = 0;
	Matrix[74] = k16;
	Matrix[75] = k19;
	Matrix[76] = 0;
	Matrix[77] = -k21;

	for (unsigned i = 0; i < 78; i++)
	{
		Matrix[i] = Matrix[i] * D0 / 30 / a / b;
	}

}

//!	Calculate element stiffness matrix 
void CPlate::ElementMass(double* Mass)
{
	clear(Mass, SizeOfStiffnessMatrix());

	CPlateMaterial* material = dynamic_cast<CPlateMaterial*>(ElementMaterial_);

	double density = material->density;
	double t = material->thick;
	double a = abs(nodes_[1]->XYZ[0] - nodes_[0]->XYZ[0]) / 2;
	double b = abs(nodes_[3]->XYZ[1] - nodes_[0]->XYZ[1]) / 2;

	double m = density * t*a*b;
	Mass[0] = 0.548253968244966*m;
	Mass[1] = 0.050793650794596*b*b*m;
	Mass[2] = 0.146349206350509*b*m;
	Mass[3] = 0.050793650794596*a*a*m;
	Mass[4] = -0.040000000001326*a*b*m;
	Mass[5] = -0.146349206350509*a*m;
	Mass[6] = 0.548253968244966*m;
	Mass[7] = -0.086984126985246*a*m;
	Mass[8] = 0.063174603177151*b*m;
	Mass[9] = 0.194603174606175*m;
	Mass[10] = 0.050793650794596*b*b*m;
	Mass[11] = 0.146349206350509*b*m;
	Mass[12] = -0.026666666667762*a*b*m;
	Mass[13] = 0.025396825398050*b*b*m;
	Mass[14] = 0.063174603177151*b*m;
	Mass[15] = 0.050793650794596*a*a*m;
	Mass[16] = 0.040000000001326*a*b*m;
	Mass[17] = 0.146349206350509*a*m;
	Mass[18] = -0.038095238095595*a*a*m;
	Mass[19] = 0.026666666667762*a*b*m;
	Mass[20] = 0.086984126985246*a*m;
	Mass[21] = 0.548253968244966*m;
	Mass[22] = 0.063174603177151*a*m;
	Mass[23] = 0.086984126985246*b*m;
	Mass[24] = 0.194603174606175*m;
	Mass[25] = -0.036825396827011*a*m;
	Mass[26] = 0.036825396827011*b*m;
	Mass[27] = 0.062539682542684*m;
	Mass[28] = 0.050793650794596*b*b*m;
	Mass[29] = -0.146349206350509*b*m;
	Mass[30] = -0.026666666667762*a*b*m;
	Mass[31] = -0.038095238095595*b*b*m;
	Mass[32] = -0.086984126985246*b*m;
	Mass[33] = 0.017777777778649*a*b*m;
	Mass[34] = -0.019047619048362*b*b*m;
	Mass[35] = -0.036825396827011*b*m;
	Mass[36] = 0.050793650794596*a*a*m;
	Mass[37] = -0.040000000001326*a*b*m;
	Mass[38] = 0.146349206350509*a*m;
	Mass[39] = 0.025396825398050*a*a*m;
	Mass[40] = 0.026666666667762*a*b*m;
	Mass[41] = 0.063174603177151*a*m;
	Mass[42] = -0.019047619048362*a*a*m;
	Mass[43] = 0.017777777778649*a*b*m;
	Mass[44] = 0.036825396827011*a*m;
	Mass[45] = 0.548253968244966*m;
	Mass[46] = 0.086984126985246*a*m;
	Mass[47] = -0.063174603177151*b*m;
	Mass[48] = 0.194603174606175*m;
	Mass[49] = 0.036825396827011*a*m;
	Mass[50] = 0.036825396827011*b*m;
	Mass[51] = 0.062539682542684*m;
	Mass[52] = -0.063174603177151*a*m;
	Mass[53] = 0.086984126985246*b*m;
	Mass[54] = 0.194603174606175*m;
	Mass[55] = 0.050793650794596*b*b*m;
	Mass[56] = -0.146349206350509*b*m;
	Mass[57] = -0.026666666667762*a*b*m;
	Mass[58] = 0.025396825398050*b*b*m;
	Mass[59] = -0.063174603177151*b*m;
	Mass[60] = -0.017777777778649*a*b*m;
	Mass[61] = -0.019047619048362*b*b*m;
	Mass[62] = -0.036825396827011*b*m;
	Mass[63] = 0.026666666667762*a*b*m;
	Mass[64] = -0.038095238095595*b*b*m;
	Mass[65] = -0.086984126985246*b*m;
	Mass[66] = 0.050793650794596*a*a*m;
	Mass[67] = 0.040000000001326*a*b*m;
	Mass[68] = -0.146349206350509*a*m;
	Mass[69] = -0.038095238095595*a*a*m;
	Mass[70] = 0.026666666667762*a*b*m;
	Mass[71] = -0.086984126985246*a*m;
	Mass[72] = -0.019047619048362*a*a*m;
	Mass[73] = -0.017777777778649*a*b*m;
	Mass[74] = -0.036825396827011*a*m;
	Mass[75] = 0.025396825398050*a*a*m;
	Mass[76] = -0.026666666667762*a*b*m;
	Mass[77] = -0.063174603177151*a*m;
}

//	Calculate element stress 
void CPlate::ElementStress(double* stress, double* Displacement)
{
	clear(stress, 32);

	CPlateMaterial* material = dynamic_cast<CPlateMaterial*>(ElementMaterial_);

	double E = material->E;
	double poisson = material->poisson;
	double t = material->thick;


	double a = abs(nodes_[1]->XYZ[0] - nodes_[0]->XYZ[0]) / 2;
	double b = abs(nodes_[3]->XYZ[1] - nodes_[0]->XYZ[1]) / 2;
	double x0 = (nodes_[1]->XYZ[0] + nodes_[0]->XYZ[0]) / 2;
	double y0 = (nodes_[3]->XYZ[1] + nodes_[0]->XYZ[1]) / 2;



	//double G  = E*t*t*t/96/(1-poisson*poisson)/a/b;
	double G = E * t / 16 / (1 - poisson * poisson) / a / b;

	double gp[2] = { -0.5773502692,0.5773502692 };      // gauss point

	double psit[4] = { -1,1,1,-1 };
	double etat[4] = { -1,-1,1,1 };

	for (unsigned int i = 0; i < 2; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			double psi = gp[i];
			double eta = gp[j];
			stress[4 * i + 2 * j + 24] = x0 + a * psi;
			stress[4 * i + 2 * j + 25] = y0 + b * eta;
			double S[36];
			for (unsigned int k = 0; k < 4; k++)
			{
				double psii = psit[k];
				double etai = etat[k];
				double psi0 = psii * psi;
				double eta0 = etai * eta;
				S[9 * k] = G * (6 * b / a * psi0*(1 + eta0) + 6 * poisson*a / b * eta0*(1 + psi0));
				S[9 * k + 1] = G * (6 * poisson*b / a * psi0*(1 + eta0) + 6 * a / b * eta0*(1 + psi0));
				S[9 * k + 2] = G * ((1 - poisson)*psii*etai*(3 * psi*psi + 3 * eta*eta - 4));
				S[9 * k + 3] = G * (-2 * poisson*a*etai*(1 + psi0)*(1 + 3 * eta0));
				S[9 * k + 4] = G * (-2 * a*etai*(1 + psi0)*(1 + 3 * eta0));
				S[9 * k + 5] = G * (-(1 - poisson)*b*psii*(3 * eta*eta + 2 * eta0 - 1));
				S[9 * k + 6] = G * 2 * b*psii*(1 + 3 * psi0)*(1 + eta0);
				S[9 * k + 7] = G * 2 * poisson*b*psii*(1 + 3 * psi0)*(1 + eta0);
				S[9 * k + 8] = G * (1 - poisson)*a*etai*(3 * psi*psi + 2 * psi0 - 1);
			}
			for (unsigned int k = 0; k < 12; k++)
			{
				if (LocationMatrix_[k])
				{
					stress[6 * i + 3 * j] += S[3 * k] * Displacement[LocationMatrix_[k] - 1];
					stress[6 * i + 3 * j + 1] += S[3 * k + 1] * Displacement[LocationMatrix_[k] - 1];
					stress[6 * i + 3 * j + 2] += S[3 * k + 2] * Displacement[LocationMatrix_[k] - 1];
				}
			}
			stress[12 + 6 * i + 3 * j] = stress[6 * i + 3 * j] / (t*t*t/12);
			stress[12 + 6 * i + 3 * j + 1] = stress[6 * i + 3 * j + 1] / (t*t*t / 12);
			stress[12 + 6 * i + 3 * j + 2] = stress[6 * i + 3 * j + 2] / (t*t*t / 12);
		}
	}
}

//	Calculate element stress for plot
void CPlate::ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* newlocation)
{
	CPlateMaterial* material = dynamic_cast<CPlateMaterial*>(ElementMaterial_);
	for (unsigned int j = 0; j < NEN_; j++)
	{
		PrePositions[3 * j] = nodes_[j]->XYZ[0];
		PrePositions[1+3 * j] = nodes_[j]->XYZ[1];
		PrePositions[2+3 * j] = nodes_[j]->XYZ[2];
		newlocation[3 * j] = nodes_[j]->XYZ[0];
		newlocation[1 + 3 * j] = nodes_[j]->XYZ[1];
		if (LocationMatrix_[3 * j])
			newlocation[2 + 3 * j] = nodes_[j]->XYZ[2] + Displacement[LocationMatrix_[3 * j] - 1];
		else
			newlocation[2 + 3 * j] = nodes_[j]->XYZ[2];

	}
}

//! Calulate element gravity
void CPlate::GravityCalculation(double* ptr_force)
{
	double g = 9.8;
	CPlateMaterial* material_ = dynamic_cast<CPlateMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double density = material_->density;
	double thick = material_->thick;
	double X[4] = { nodes_[0]->XYZ[0],nodes_[1]->XYZ[0],nodes_[2]->XYZ[0],nodes_[3]->XYZ[0] };
	double Y[4] = { nodes_[0]->XYZ[1],nodes_[1]->XYZ[1],nodes_[2]->XYZ[1],nodes_[3]->XYZ[1] };
	double Z[4] = { nodes_[0]->XYZ[2],nodes_[1]->XYZ[2],nodes_[2]->XYZ[2],nodes_[3]->XYZ[2] };

	//Calculate the unit normal vector for each node
	double l3[4], m3[4], n3[4];
	l3[0] = (Z[0] * (Y[0] - Y[1]) - Z[0] * (Y[0] - Y[3]) + Z[1] * (Y[0] - Y[3]) - Z[3] * (Y[0] - Y[1])) / 4;
	m3[0] = (X[0] * (Z[0] - Z[1]) - X[0] * (Z[0] - Z[3]) + X[1] * (Z[0] - Z[3]) - X[3] * (Z[0] - Z[1])) / 4;
	n3[0] = (Y[0] * (X[0] - X[1]) - Y[0] * (X[0] - X[3]) + Y[1] * (X[0] - X[3]) - Y[3] * (X[0] - X[1])) / 4;
	l3[1] = (Z[1] * (Y[0] - Y[1]) - Z[0] * (Y[1] - Y[2]) - Z[2] * (Y[0] - Y[1]) + Z[1] * (Y[1] - Y[2])) / 4;
	m3[1] = (X[1] * (Z[0] - Z[1]) - X[0] * (Z[1] - Z[2]) - X[2] * (Z[0] - Z[1]) + X[1] * (Z[1] - Z[2])) / 4;
	n3[1] = (Y[1] * (X[0] - X[1]) - Y[0] * (X[1] - X[2]) - Y[2] * (X[0] - X[1]) + Y[1] * (X[1] - X[2])) / 4;
	l3[2] = (Z[2] * (Y[1] - Y[2]) - Z[1] * (Y[2] - Y[3]) - Z[3] * (Y[1] - Y[2]) + Z[2] * (Y[2] - Y[3])) / 4;
	m3[2] = (X[2] * (Z[1] - Z[2]) - X[1] * (Z[2] - Z[3]) - X[3] * (Z[1] - Z[2]) + X[2] * (Z[2] - Z[3])) / 4;
	n3[2] = (Y[2] * (X[1] - X[2]) - Y[1] * (X[2] - X[3]) - Y[3] * (X[1] - X[2]) + Y[2] * (X[2] - X[3])) / 4;
	l3[3] = (Z[2] * (Y[0] - Y[3]) - Z[0] * (Y[2] - Y[3]) - Z[3] * (Y[0] - Y[3]) + Z[3] * (Y[2] - Y[3])) / 4;
	m3[3] = (X[2] * (Z[0] - Z[3]) - X[0] * (Z[2] - Z[3]) - X[3] * (Z[0] - Z[3]) + X[3] * (Z[2] - Z[3])) / 4;
	n3[3] = (Y[2] * (X[0] - X[3]) - Y[0] * (X[2] - X[3]) - Y[3] * (X[0] - X[3]) + Y[3] * (X[2] - X[3])) / 4;
	for (unsigned int loop = 0; loop < 4; loop++) {
		double lmn3 = sqrt(l3[loop] * l3[loop] + m3[loop] * m3[loop] + n3[loop] * n3[loop]);
		l3[loop] /= lmn3;
		m3[loop] /= lmn3;
		n3[loop] /= lmn3;
	}
	double l1[4], m1[4], n1[4];
	double l2[4], m2[4], n2[4];
	for (unsigned int loop = 0; loop < 4; loop++) {
		double lmn1 = sqrt(n3[loop] * n3[loop] + l3[loop] * l3[loop]);
		l1[loop] = n3[loop] / lmn1;
		m1[loop] = 0;
		n1[loop] = -l3[loop] / lmn1;
		l2[loop] = -m3[loop] * l3[loop] / lmn1;
		m2[loop] = (n3[loop] * n3[loop] + l3[loop] * l3[loop]) / lmn1;
		n2[loop] = -m3[loop] * n2[loop] / lmn1;
	}

	double gausspoint[2] = { -0.57735027, 0.57735027 };	//Gauss point when ngp=2
	double weight_gauss[2] = { 1,1 };	//Gauss weight

	double volumn = 0;
	for (unsigned int i = 0; i < 2; i++) {
		double xi = gausspoint[i];
		for (unsigned int j = 0; j < 2; j++) {
			double eta = gausspoint[j];
			for (unsigned int k = 0; k < 2; k++) {
				double zeta = gausspoint[k];
				double N[4] = { (1 - xi)*(1 - eta) / 4,(1 + xi)*(1 - eta) / 4,(1 + xi)*(1 + eta) / 4,(1 - xi)*(1 + eta) / 4 };
				double GNxi[4] = { -(1 - eta) / 4,(1 - eta) / 4,(1 + eta) / 4,-(1 + eta) / 4 };		//GNxi=dN/dxi
				double GNeta[4] = { -(1 - xi) / 4,-(1 + xi) / 4,(1 + xi) / 4,(1 - xi) / 4 };		//GNeta=dN/deta
				double J[3][3] = { 0,0,0,0,0,0,0,0,0 };			//Calculate the Jocabian Matrix
				for (unsigned int loop = 0; loop < 4; loop++) {
					J[0][0] += GNxi[loop] * (X[loop] + 0.5*zeta*thick*l3[loop]);
					J[0][1] += GNxi[loop] * (Y[loop] + 0.5*zeta*thick*m3[loop]);
					J[0][2] += GNxi[loop] * (Z[loop] + 0.5*zeta*thick*n3[loop]);
					J[1][0] += GNeta[loop] * (X[loop] + 0.5*zeta*thick*l3[loop]);
					J[1][1] += GNeta[loop] * (Y[loop] + 0.5*zeta*thick*m3[loop]);
					J[1][2] += GNeta[loop] * (Z[loop] + 0.5*zeta*thick*n3[loop]);
					J[2][0] += N[loop] * 0.5*thick*l3[loop];
					J[2][1] += N[loop] * 0.5*thick*m3[loop];
					J[2][2] += N[loop] * 0.5*thick*n3[loop];
				}
				double detJ = J[0][0] * J[1][1] * J[2][2] - J[0][0] * J[1][2] * J[2][1] - J[0][1] * J[1][0] * J[2][2] + J[0][1] * J[1][2] * J[2][0] + J[0][2] * J[1][0] * J[2][1] - J[0][2] * J[1][1] * J[2][0];
				volumn += fabs(detJ)*weight_gauss[i] * weight_gauss[j] * weight_gauss[k];


			}
		}
	}
	weight = volumn * g*density;
	double weight_avg = weight / 4;
	double Xc = (X[0] + X[1] + X[2] + X[3]) / 4;
	double Yc = (Y[0] + Y[1] + Y[2] + Y[3]) / 4;
	double Zc = (Z[0] + Z[1] + Z[2] + Z[3]) / 4;
	for (unsigned int loop = 0; loop < 4; loop++) {
		ptr_force[3 * loop] = -weight_avg;
		ptr_force[3 * loop + 1] = 2 * weight_avg / 3 * (Yc - Y[loop]);
		ptr_force[3 * loop + 2] = -2 * weight_avg / 3 * (Xc - X[loop]);
	}

}


//	Recover element stress 
void CPlate::RecoverElementStress(double* Displacement, double* A)
{

	CPlateMaterial* material = dynamic_cast<CPlateMaterial*>(ElementMaterial_);

	double E = material->E;
	double poisson = material->poisson;
	double t = material->thick;

	double a = abs(nodes_[1]->XYZ[0] - nodes_[0]->XYZ[0]) / 2;
	double b = abs(nodes_[3]->XYZ[1] - nodes_[0]->XYZ[1]) / 2;
	double x0 = (nodes_[1]->XYZ[0] + nodes_[0]->XYZ[0]) / 2;
	double y0 = (nodes_[3]->XYZ[1] + nodes_[0]->XYZ[1]) / 2;
	//double G  = E*t*t*t/96/(1-poisson*poisson)/a/b;
	double G = E * t / 16 / (1 - poisson * poisson) / a / b;

	double gp[2] = { -0.5773502692,0.5773502692 };      // gauss point

	double psit[4] = { -1,1,1,-1 };
	double etat[4] = { -1,-1,1,1 };

	int Node[4];

	Node[0] = nodes_[0]->NodeNumber - 1;
	Node[1] = nodes_[1]->NodeNumber - 1;
	Node[2] = nodes_[2]->NodeNumber - 1;
	Node[3] = nodes_[3]->NodeNumber - 1;


	for (unsigned int i = 0; i < 2; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			double psi = gp[i];
			double eta = gp[j];

			A[9 * (i*i + 3 * j*j - 2 * i*j) + 73 * Node[i*i + 3 * j*j - 2 * i*j]] = x0 + a * psi - nodes_[i*i + 3 * j*j - 2 * i*j]->XYZ[0];
			A[1 + 9 * (i*i + 3 * j*j - 2 * i*j) + 73 * Node[i*i + 3 * j*j - 2 * i*j]] = y0 + b * eta - nodes_[i*i + 3 * j*j - 2 * i*j]->XYZ[1];
			A[72 + 73 * Node[i*i + 3 * j*j - 2 * i*j]] += 1;


			double S[36];
			for (unsigned int k = 0; k < 4; k++)
			{
				double psii = psit[k];
				double etai = etat[k];
				double psi0 = psii * psi;
				double eta0 = etai * eta;
				S[9 * k] = G * (6 * b / a * psi0*(1 + eta0) + 6 * poisson*a / b * eta0*(1 + psi0));
				S[9 * k + 1] = G * (6 * poisson*b / a * psi0*(1 + eta0) + 6 * a / b * eta0*(1 + psi0));
				S[9 * k + 2] = G * ((1 - poisson)*psii*etai*(3 * psi*psi + 3 * eta*eta - 4));
				S[9 * k + 3] = G * (-2 * poisson*a*etai*(1 + psi0)*(1 + 3 * eta0));
				S[9 * k + 4] = G * (-2 * a*etai*(1 + psi0)*(1 + 3 * eta0));
				S[9 * k + 5] = G * (-(1 - poisson)*b*psii*(3 * eta*eta + 2 * eta0 - 1));
				S[9 * k + 6] = G * 2 * b*psii*(1 + 3 * psi0)*(1 + eta0);
				S[9 * k + 7] = G * 2 * poisson*b*psii*(1 + 3 * psi0)*(1 + eta0);
				S[9 * k + 8] = G * (1 - poisson)*a*etai*(3 * psi*psi + 2 * psi0 - 1);
			}
			for (unsigned int k = 0; k < 12; k++)
			{
				if (LocationMatrix_[k])
				{
					A[9 * (i*i + 3 * j*j - 2 * i*j) + 73 * Node[i*i + 3 * j*j - 2 * i*j] + 3] += S[3 * k] * Displacement[LocationMatrix_[k] - 1];
					A[9 * (i*i + 3 * j*j - 2 * i*j) + 73 * Node[i*i + 3 * j*j - 2 * i*j] + 4] += S[3 * k + 1] * Displacement[LocationMatrix_[k] - 1];
					A[9 * (i*i + 3 * j*j - 2 * i*j) + 73 * Node[i*i + 3 * j*j - 2 * i*j] + 6] += S[3 * k + 2] * Displacement[LocationMatrix_[k] - 1];
				}
			}
		}
	}
}