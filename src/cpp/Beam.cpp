/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Beam.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

//	Constructor
CBeam::CBeam()
{
	NEN_ = 2; // Each element has 2 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CBeam::~CBeam() 
{
	delete[]nodes_;
	delete[]LocationMatrix_;
}

//	Read element data from stream Input
bool CBeam::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
void CBeam::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber 
		<< setw(9)<< nodes_[1]->NodeNumber 
		<< setw(12) << ElementMaterial_->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the
//  element
//	Caution:  Equation number is numbered from 1 !
void CBeam::GenerateLocationMatrix()
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
unsigned int CBeam::SizeOfStiffnessMatrix() { return 78; }

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBeam::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	//	Calculate beam length
	double DX[3]; //	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

	// the foundamental constant of the element stiffness matrix
	const CBeamMaterial& material =
		dynamic_cast<CBeamMaterial&>(*ElementMaterial_); // Pointer to material of the element
	double Iz = (material.width * material.width * material.width) * material.height / 12 -
		pow(material.width - 2 * material.t_side, 3.0) * (material.height - 2 * material.t_uplow) / 12;
	double Iy = (material.height * material.height * material.height) * material.width / 12 -
		pow(material.height - 2 * material.t_uplow, 3.0) * (material.width - 2 * material.t_side) / 12;
	double Ip = Iz + Iy;
	double material_const_1 = material.E * material.width * material.height / L;
	double material_const_2 = 12 * material.E * Iz / (L * L * L);
	double material_const_3 = 12 * material.E * Iy / (L * L * L);
	double material_const_4 = material.E * Ip / ((2 + 2 * material.mu) * L);
	double material_const_5 = 4 * material.E * Iy / L;
	double material_const_6 = 4 * material.E * Iz / L;
	double material_const_7 = 6 * material.E * Iy / (L * L);
	double material_const_8 = 6 * material.E * Iz / (L * L);

	double n[3][3];

	for (unsigned int i = 0; i < 3; i++)
		n[0][i] = DX[i] / L; //orientation of x'

	double len_orientation;
	len_orientation = sqrt(material.n_x*material.n_x + material.n_y*material.n_y + material.n_z*material.n_z);

    // n matrix is the orientation cosine matrix
	n[1][0] = material.n_x/len_orientation;
	n[1][1] = material.n_y/len_orientation;
	n[1][2] = material.n_z/len_orientation; // orientation of y'
	n[2][0] = n[0][1] * n[1][2] - n[0][2] * n[1][1];
	n[2][1] = n[0][2] * n[1][0] - n[0][0] * n[1][2];
	n[2][2] = n[0][0] * n[1][1] - n[0][1] * n[1][0]; //orientation of z'

	// the set of product results from orientation cosine matrix
	// A'K A
	double orient_pconst1 = n[0][0] * n[0][0];
	double orient_pconst2 = n[0][1] * n[0][1];
	double orient_pconst3 = n[0][2] * n[0][2];
	double orient_pconst4 = n[1][0] * n[1][0];
	double orient_pconst5 = n[1][1] * n[1][1];
	double orient_pconst6 = n[1][2] * n[1][2];
	double orient_pconst7 = n[2][0] * n[2][0];
	double orient_pconst8 = n[2][1] * n[2][1];
	double orient_pconst9 = n[2][2] * n[2][2];
	double orient_pconst10 = n[0][0] * n[0][1];
	double orient_pconst11 = n[1][0] * n[1][1];
	double orient_pconst12 = n[2][0] * n[2][1];
	double orient_pconst13 = n[0][1] * n[0][2];
	double orient_pconst14 = n[1][1] * n[1][2];
	double orient_pconst15 = n[2][1] * n[2][2];
	double orient_pconst16 = n[0][0] * n[0][2];
	double orient_pconst17 = n[1][0] * n[1][2];
	double orient_pconst18 = n[2][0] * n[2][2];
	double orient_pconst19 = n[1][2] * n[2][0];
	double orient_pconst20 = n[1][0] * n[2][2];
	double orient_pconst21 = n[1][1] * n[2][0];
	double orient_pconst22 = n[1][0] * n[2][1];
	double orient_pconst23 = n[1][2] * n[2][1];
	double orient_pconst24 = n[1][1] * n[2][2];
	double orient_pconst25 = n[1][0] * n[2][0];
	double orient_pconst26 = n[1][1] * n[2][1];
	double orient_pconst27 = n[1][2] * n[2][2];

	// element stiffness matrix
	Matrix[0] = material_const_1 * orient_pconst1 + material_const_2 * orient_pconst4 + material_const_3 * orient_pconst7;
	Matrix[1] = material_const_1 * orient_pconst2 + material_const_2 * orient_pconst5 + material_const_3 * orient_pconst8;
	Matrix[2] = material_const_1 * orient_pconst10 + material_const_2 * orient_pconst11 + material_const_3 * orient_pconst12;
	Matrix[3] = material_const_1 * orient_pconst3 + material_const_2 * orient_pconst6 + material_const_3 * orient_pconst9;
	Matrix[4] = material_const_1 * orient_pconst13 + material_const_2 * orient_pconst14 + material_const_3 * orient_pconst15;
	Matrix[5] = material_const_1 * orient_pconst16 + material_const_2 * orient_pconst17 + material_const_3 * orient_pconst18;
	Matrix[6] = material_const_4 * orient_pconst1 + material_const_5 * orient_pconst4 + material_const_6 * orient_pconst7;
	Matrix[7] = material_const_8 * orient_pconst19 - material_const_7 * orient_pconst20;
	Matrix[8] = material_const_8 * orient_pconst21 - material_const_7 * orient_pconst22;
	Matrix[9] = material_const_8 * orient_pconst25 - material_const_7 * orient_pconst25;
	Matrix[10] = material_const_4 * orient_pconst2 + material_const_5 * orient_pconst5 + material_const_6 * orient_pconst8;
	Matrix[11] = material_const_4 * orient_pconst10 + material_const_5 * orient_pconst11 + material_const_6 * orient_pconst12;
	Matrix[12] = material_const_8 * orient_pconst23 - material_const_7 * orient_pconst24;
	Matrix[13] = material_const_8 * orient_pconst26 - material_const_7 * orient_pconst26;
	Matrix[14] = material_const_8 * orient_pconst22 - material_const_7 * orient_pconst21;
	Matrix[15] = material_const_4 * orient_pconst3 + material_const_5 * orient_pconst6 + material_const_6 * orient_pconst9;
	Matrix[16] = material_const_4 * orient_pconst13 + material_const_5 * orient_pconst14 + material_const_6 * orient_pconst15;
	Matrix[17] = material_const_4 * orient_pconst16 + material_const_5 * orient_pconst17 + material_const_6 * orient_pconst18;
	Matrix[18] = material_const_8 * orient_pconst27 - material_const_7 * orient_pconst27;
	Matrix[19] = material_const_8 * orient_pconst24 - material_const_7 * orient_pconst23;
	Matrix[20] = material_const_8 * orient_pconst20 - material_const_7 * orient_pconst19;
	Matrix[21] = material_const_1 * orient_pconst1 + material_const_2 * orient_pconst4 + material_const_3 * orient_pconst7;
	Matrix[22] = material_const_7 * orient_pconst19 - material_const_8 * orient_pconst20;
	Matrix[23] = material_const_7 * orient_pconst21 - material_const_8 * orient_pconst22;
	Matrix[24] = material_const_7 * orient_pconst25 - material_const_8 * orient_pconst25;
	Matrix[25] = -material_const_1 * orient_pconst16 - material_const_2 * orient_pconst17 - material_const_3 * orient_pconst18;
	Matrix[26] = -material_const_1 * orient_pconst10 - material_const_2 * orient_pconst11 - material_const_3 * orient_pconst12;
	Matrix[27] = -material_const_1 * orient_pconst1 - material_const_2 * orient_pconst4 - material_const_3 * orient_pconst7;
	Matrix[28] = material_const_1 * orient_pconst2 + material_const_2 * orient_pconst5 + material_const_3 * orient_pconst8;
	Matrix[29] = material_const_1 * orient_pconst10 + material_const_2 * orient_pconst11 + material_const_3 * orient_pconst12;
	Matrix[30] = material_const_7 * orient_pconst23 - material_const_8 * orient_pconst24;
	Matrix[31] = material_const_7 * orient_pconst26 - material_const_8 * orient_pconst26;
	Matrix[32] = material_const_7 * orient_pconst22 - material_const_8 * orient_pconst21;
	Matrix[33] = -material_const_1 * orient_pconst13 - material_const_2 * orient_pconst14 - material_const_3 * orient_pconst15;
	Matrix[34] = -material_const_1 * orient_pconst2 - material_const_2 * orient_pconst5 - material_const_3 * orient_pconst8;
	Matrix[35] = -material_const_1 * orient_pconst10 - material_const_2 * orient_pconst11 - material_const_3 * orient_pconst12;
	Matrix[36] = material_const_1 * orient_pconst3 + material_const_2 * orient_pconst6 + material_const_3 * orient_pconst9;
	Matrix[37] = material_const_1 * orient_pconst13 + material_const_2 * orient_pconst14 + material_const_3 * orient_pconst15;
	Matrix[38] = material_const_1 * orient_pconst16 + material_const_2 * orient_pconst17 + material_const_3 * orient_pconst18;
	Matrix[39] = material_const_7 * orient_pconst27 - material_const_8 * orient_pconst27;
	Matrix[40] = material_const_7 * orient_pconst24 - material_const_8 * orient_pconst23;
	Matrix[41] = material_const_7 * orient_pconst20 - material_const_8 * orient_pconst19;
	Matrix[42] = -material_const_1 * orient_pconst3 - material_const_2 * orient_pconst6 - material_const_3 * orient_pconst9;
	Matrix[43] = -material_const_1 * orient_pconst13 - material_const_2 * orient_pconst14 - material_const_3 * orient_pconst15;
	Matrix[44] = -material_const_1 * orient_pconst16 - material_const_2 * orient_pconst17 - material_const_3 * orient_pconst18;
	Matrix[45] = material_const_4 * orient_pconst1 + material_const_5 * orient_pconst4 + material_const_6 * orient_pconst7;
	Matrix[46] = material_const_7 * orient_pconst20 - material_const_8 * orient_pconst19;
	Matrix[47] = material_const_7 * orient_pconst22 - material_const_8 * orient_pconst21;
	Matrix[48] = material_const_7 * orient_pconst25 - material_const_8 * orient_pconst25;
	Matrix[49] = (material_const_5 * orient_pconst17) / 2 - material_const_4 * orient_pconst16 + (material_const_6 * orient_pconst18) / 2;
	Matrix[50] = (material_const_5 * orient_pconst11) / 2 - material_const_4 * orient_pconst10 + (material_const_6 * orient_pconst12) / 2;
	Matrix[51] = -material_const_4 * orient_pconst1 + (material_const_5 * orient_pconst4) / 2 + (material_const_6 * orient_pconst7) / 2;
	Matrix[52] = material_const_8 * orient_pconst19 - material_const_7 * orient_pconst20;
	Matrix[53] = material_const_8 * orient_pconst21 - material_const_7 * orient_pconst22;
	Matrix[54] = material_const_8 * orient_pconst25 - material_const_7 * orient_pconst25;
	Matrix[55] = material_const_4 * orient_pconst2 + material_const_5 * orient_pconst5 + material_const_6 * orient_pconst8;
	Matrix[56] = material_const_4 * orient_pconst10 + material_const_5 * orient_pconst11 + material_const_6 * orient_pconst12;
	Matrix[57] = material_const_7 * orient_pconst24 - material_const_8 * orient_pconst23;
	Matrix[58] = material_const_7 * orient_pconst26 - material_const_8 * orient_pconst26;
	Matrix[59] = material_const_7 * orient_pconst21 - material_const_8 * orient_pconst22;
	Matrix[60] = (material_const_5 * orient_pconst14) / 2 - material_const_4 * orient_pconst13 + (material_const_6 * orient_pconst15) / 2;
	Matrix[61] = -material_const_4 * orient_pconst2 + (material_const_5 * orient_pconst5) / 2 + (material_const_6 * orient_pconst8) / 2;
	Matrix[62] = (material_const_5 * orient_pconst11) / 2 - material_const_4 * orient_pconst10 + (material_const_6 * orient_pconst12) / 2;
	Matrix[63] = material_const_8 * orient_pconst23 - material_const_7 * orient_pconst24;
	Matrix[64] = material_const_8 * orient_pconst26 - material_const_7 * orient_pconst26;
	Matrix[65] = material_const_8 * orient_pconst22 - material_const_7 * orient_pconst21;
	Matrix[66] = material_const_4 * orient_pconst3 + material_const_5 * orient_pconst6 + material_const_6 * orient_pconst9;
	Matrix[67] = material_const_4 * orient_pconst13 + material_const_5 * orient_pconst14 + material_const_6 * orient_pconst15;
	Matrix[68] = material_const_4 * orient_pconst16 + material_const_5 * orient_pconst17 + material_const_6 * orient_pconst18;
	Matrix[69] = material_const_7 * orient_pconst27 - material_const_8 * orient_pconst27;
	Matrix[70] = material_const_7 * orient_pconst23 - material_const_8 * orient_pconst24;
	Matrix[71] = material_const_7 * orient_pconst19 - material_const_8 * orient_pconst20;
	Matrix[72] = -material_const_4 * orient_pconst3 + (material_const_5 * orient_pconst6) / 2 + (material_const_6 * orient_pconst9) / 2;
	Matrix[73] = (material_const_5 * orient_pconst14) / 2 - material_const_4 * orient_pconst13 + (material_const_6 * orient_pconst15) / 2;
	Matrix[74] = (material_const_5 * orient_pconst17) / 2 - material_const_4 * orient_pconst16 + (material_const_6 * orient_pconst18) / 2;
	Matrix[75] = material_const_8 * orient_pconst27 - material_const_7 * orient_pconst27;
	Matrix[76] = material_const_8 * orient_pconst24 - material_const_7 * orient_pconst23;
	Matrix[77] = material_const_8 * orient_pconst20 - material_const_7 * orient_pconst19;
}
//	Calculate element stress
void CBeam::ElementStress(double* stress, double* Displacement)
{
	const CBeamMaterial& material =
		static_cast<CBeamMaterial&>(*ElementMaterial_); // Pointer to material of the element
	clear(stress, 3);
	double DX[3]; //	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L = 0; //	 beam length

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
	}

	L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

	// stress matrix
	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * DX[i] * material.E / (L * L * L);
		S[i + 3] = -S[i];
	}

	for (unsigned int i = 0; i < 2; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			if (LocationMatrix_[i * 6 + j])
				stress[j] += S[i * 3 + j] * Displacement[LocationMatrix_[i * 6 + j] - 1];
		}
	}
}

void CBeam::ElementMass(double* mass)
{
	clear(mass, 1);
}

void CBeam::ElementPostInfo(double* beamstress, double* Displacement, double* prePositionBeam, double* postPositionBeam)
{
	const CBeamMaterial& material =
		dynamic_cast<CBeamMaterial&>(*ElementMaterial_); // Pointer to material of the element
	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++) {
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i] * DX[i];
	}

	double n[3][3];
	double L = sqrt(L2);
	for (unsigned int i = 0; i < 3; i++)
		n[0][i] = DX[i] / L; //orientation of x'

	double len_orientation;
	len_orientation = sqrt(material.n_x*material.n_x + material.n_y*material.n_y + material.n_z*material.n_z);

	n[1][0] = material.n_x/len_orientation;
	n[1][1] = material.n_y/len_orientation;
	n[1][2] = material.n_z/len_orientation; // orientation of y'
	n[2][0] = n[0][1] * n[1][2] - n[0][2] * n[1][1];
	n[2][1] = n[0][2] * n[1][0] - n[0][0] * n[1][2];
	n[2][2] = n[0][0] * n[1][1] - n[0][1] * n[1][0]; //orientation of z'

	double pos_origin[2][3]; // preposition
	double local_diag[4][3];//vector from center of rectangle to four angles in the local coordinate
	double global_diag[4][3];//vector from center of rectangle to four angles in the main coordinate
	double Iz = (material.width * material.width * material.width) * material.height / 12 -
		pow(material.width - 2 * material.t_side, 3.0) *
		(material.height - 2 * material.t_uplow) / 12;
	double Iy = (material.height * material.height * material.height) * material.width / 12 -
		pow(material.height - 2 * material.t_uplow, 3.0) *
		(material.width - 2 * material.t_side) / 12;
	double Ip = Iz + Iy;

	// Define the scale of plot
	double magCodim = 1E-1;
	local_diag[0][0] = 0;
	local_diag[1][0] = 0;
	local_diag[2][0] = 0;
	local_diag[3][0] = 0;
	local_diag[0][1] = -magCodim * material.width;
	local_diag[1][1] = magCodim * material.width;
	local_diag[2][1] = magCodim * material.width;
	local_diag[3][1] = -magCodim * material.width;
	local_diag[0][2] = magCodim * material.height;
	local_diag[1][2] = magCodim * material.height;
	local_diag[2][2] = -magCodim * material.height;
	local_diag[3][2] = -magCodim * material.height;

	for (unsigned int i = 0; i < 4; i++) 
	{
		global_diag[i][0] = n[0][0] * local_diag[i][0] + n[1][0] * local_diag[i][1] + n[2][0] * local_diag[i][2];
		global_diag[i][1] = n[0][1] * local_diag[i][0] + n[1][1] * local_diag[i][1] + n[2][1] * local_diag[i][2];
		global_diag[i][2] = n[0][2] * local_diag[i][0] + n[1][2] * local_diag[i][1] + n[2][2] * local_diag[i][2];
	}

	double global_d[2][3]; //displacement in the global coordinate
	double element_d[2][3]; //displacement in the local coordinate
	double global_theta[2][3]; //rotation in the global coordinate

	// Loop over displacements
	for (unsigned int i = 0; i < 3; i++) 
	{

		if (LocationMatrix_[i]) {
			global_d[0][i] = Displacement[LocationMatrix_[i] - 1];
			pos_origin[0][i] = nodes_[0]->XYZ[i];
		}
		else {
			global_d[0][i] = 0.0;
			pos_origin[0][i] = nodes_[0]->XYZ[i];
		}

		if (LocationMatrix_[i + 6]) {
			global_d[1][i] = Displacement[LocationMatrix_[i + 6] - 1];
			pos_origin[1][i] = nodes_[1]->XYZ[i];
		}
		else {
			global_d[1][i] = 0.0;
			pos_origin[1][i] = nodes_[1]->XYZ[i];
		}
	}

	// Loop over rotations
	for (unsigned int i = 3; i < 6; i++) {

		if (LocationMatrix_[i]) {
			global_theta[0][i - 3] = Displacement[LocationMatrix_[i] - 1];
		}
		else {
			global_theta[0][i - 3] = nodes_[0]->XYZ[i];
		}

		if (LocationMatrix_[i + 6]) {
			global_theta[1][i - 3] = Displacement[LocationMatrix_[i + 6] - 1];
		}
		else {
			global_theta[1][i - 3] = nodes_[1]->XYZ[i];
		}
	}

	for (unsigned int i = 0; i < 2; i++) {
		for (unsigned int j = 0; j < 4; j++) 
		{
			postPositionBeam[(i * 4 + j) * 3] = pos_origin[i][0] + global_d[i][0] + global_diag[j][0] + global_diag[j][2] * global_theta[i][1] - global_diag[j][1] * global_theta[i][2];
			postPositionBeam[(i * 4 + j) * 3 + 1] = pos_origin[i][1] + global_d[i][1] + global_diag[j][1] + global_diag[j][0] * global_theta[i][2] - global_diag[j][2] * global_theta[i][0];
			postPositionBeam[(i * 4 + j) * 3 + 2] = pos_origin[i][2] + global_d[i][2] + global_diag[j][2] + global_diag[j][1] * global_theta[i][0] - global_diag[j][0] * global_theta[i][1];
			prePositionBeam[(i * 4 + j) * 3] = pos_origin[i][0] + global_diag[j][0];
			prePositionBeam[(i * 4 + j) * 3 + 1] = pos_origin[i][1] + global_diag[j][1];
			prePositionBeam[(i * 4 + j) * 3 + 2] = pos_origin[i][2] + global_diag[j][2];
		}
	}

	double element_theta[2][3]; //rotation in the local coordinate
	//coordinate conversion
	for (unsigned int i = 0; i < 2; i++) 
	{
		element_d[i][0] = n[0][0] * global_d[i][0] + n[1][0] * global_d[i][1] + n[2][0] * global_d[i][2];
		element_d[i][1] = n[0][1] * global_d[i][0] + n[1][1] * global_d[i][1] + n[2][1] * global_d[i][2];
		element_d[i][2] = n[0][2] * global_d[i][0] + n[1][2] * global_d[i][1] + n[2][2] * global_d[i][2];
		element_theta[i][0] = n[0][0] * global_theta[i][0] + n[1][0] * global_theta[i][1] + n[2][0] * global_theta[i][2];
		element_theta[i][1] = n[0][1] * global_theta[i][0] + n[1][1] * global_theta[i][1] + n[2][1] * global_theta[i][2];
		element_theta[i][2] = n[0][2] * global_theta[i][0] + n[1][2] * global_theta[i][1] + n[2][2] * global_theta[i][2];
	}

	double dtheta[2][3];//rate of the change of the corner
	dtheta[0][0] = (element_theta[1][2] - element_theta[0][2]) / L;
	dtheta[1][0] = dtheta[0][0];
	dtheta[0][1] = (element_d[0][2] - element_d[1][2]) * 6 / L2 - element_theta[0][1] * 4 / L + element_theta[1][1] * 2 / L;
	dtheta[1][1] = (element_d[1][2] - element_d[0][2]) * 6 / L2 + element_theta[0][1] * 2 / L - element_theta[1][1] * 4 / L;
	dtheta[0][2] = (element_d[1][1] - element_d[0][1]) * 6 / L2 - element_theta[0][2] * 4 / L + element_theta[1][2] * 2 / L;
	dtheta[1][2] = (element_d[0][1] - element_d[1][1]) * 6 / L2 + element_theta[0][2] * 2 / L - element_theta[1][2] * 4 / L;

	double sigma1;       //Normal stress caused by strech
	double sigma2[2][2]; //Normal stress caused by bending(z-bending and y-bending)
	double tau_xy;       //Shearing Strength
	double tau_xz;

	sigma1 = material.E * (element_d[1][0] - element_d[0][0]) / L;
	tau_xy = material.E * material.height * dtheta[0][0] / (4 + 4 * material.mu);
	tau_xz = material.E * material.width * dtheta[0][0] / (4 + 4 * material.mu);

	// beam stress is aligned in the following way
	// 1 left up  2 left down  3 right down  4 right up
	// sigma xx yy zz xy yz xz
	for (unsigned int i = 0; i < 2; i++) 
	{
		sigma2[i][0] = material.E * material.width * dtheta[i][2] / 2;
		sigma2[i][1] = material.E * material.height * dtheta[i][1] / 2;

		beamstress[i * 24] = sigma1 - sigma2[i][0] + sigma2[i][1];
		beamstress[i * 24 + 1] = 0;
		beamstress[i * 24 + 2] = 0;
		beamstress[i * 24 + 3] = -tau_xy;
		beamstress[i * 24 + 4] = 0;
		beamstress[i * 24 + 5] = tau_xz;
		beamstress[i * 24 + 6] = sigma1 - sigma2[i][0] - sigma2[i][1];
		beamstress[i * 24 + 7] = 0;
		beamstress[i * 24 + 8] = 0;
		beamstress[i * 24 + 9] = tau_xy;
		beamstress[i * 24 + 10] = 0;
		beamstress[i * 24 + 11] = tau_xz;
		beamstress[i * 24 + 12] = sigma1 + sigma2[i][0] - sigma2[i][1];
		beamstress[i * 24 + 13] = 0;
		beamstress[i * 24 + 14] = 0;
		beamstress[i * 24 + 15] = tau_xy;
		beamstress[i * 24 + 16] = 0;
		beamstress[i * 24 + 17] = -tau_xz;
		beamstress[i * 24 + 18] = sigma1 + sigma2[i][0] + sigma2[i][1];
		beamstress[i * 24 + 19] = 0;
		beamstress[i * 24 + 20] = 0;
		beamstress[i * 24 + 21] = -tau_xy;
		beamstress[i * 24 + 22] = 0;
		beamstress[i * 24 + 23] = -tau_xz;
	}
}

void CBeam::GravityCalculation(double* ptr_force)
{
	double g = 9.8;
	CBeamMaterial* material_ = dynamic_cast<CBeamMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double density = material_->density;
	double DX[3]; //	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
	double L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);
	double A = material_->width*material_->height - (material_->width - 2 * material_->t_side)*(material_->height - 2 * material_->t_uplow);
	weight = density * L*A*g;
	ptr_force[0] = -weight / 2;
	ptr_force[3] = ptr_force[0];
	ptr_force[1] = L * g*density*A*DX[0];
	ptr_force[2] = -L * g*density*A*DX[1];
	ptr_force[4] = -L * g*density*A*DX[0];
	ptr_force[5] = L * g*density*A*DX[1];
}