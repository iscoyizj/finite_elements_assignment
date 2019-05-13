/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Subpara.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CSubpara::CSubpara()
{
	NEN_ = 9;	// Each element has 9 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 18;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CSubpara::~CSubpara()
{
}

//	Read element data from stream Input
bool CSubpara::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8, N9;	// nine node numbers

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> N9 >> MSet;
    ElementMaterial_ = dynamic_cast<CSubparaMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[4] = &NodeList[N5 - 1];
	nodes_[5] = &NodeList[N6 - 1];
	nodes_[6] = &NodeList[N7 - 1];
	nodes_[7] = &NodeList[N8 - 1];
	nodes_[8] = &NodeList[N9 - 1];
	for (unsigned int i = 0; i < 2; i++)                                                  //计算各边中间点的坐标
	{
		nodes_[4]->XYZ[i] = 0.5*(nodes_[0]->XYZ[i] + nodes_[1]->XYZ[i]);
		nodes_[5]->XYZ[i] = 0.5*(nodes_[1]->XYZ[i] + nodes_[2]->XYZ[i]);
		nodes_[6]->XYZ[i] = 0.5*(nodes_[2]->XYZ[i] + nodes_[3]->XYZ[i]);
		nodes_[7]->XYZ[i] = 0.5*(nodes_[3]->XYZ[i] + nodes_[0]->XYZ[i]);
		nodes_[8]->XYZ[i] = 0.5*(nodes_[5]->XYZ[i] + nodes_[7]->XYZ[i]);
	}
	return true;
}

//	Write element data to stream
void CSubpara::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1;
	for (unsigned int i = 0; i < 9; i++)
	{
		output << setw(9) << nodes_[i]->NodeNumber;
	}
	output << setw(12) << ElementMaterial_->nset << endl;
}

void CSubpara::WritePlot(COutPlot& output, unsigned int Ele)
{
	output << 4 << setw(11) << nodes_[0]->NodeNumber - 1 << nodes_[1]->NodeNumber - 1 << setw(9) << nodes_[2]->NodeNumber - 1 << setw(9)
		<< nodes_[3]->NodeNumber - 1 << endl;
}

void CSubpara::WritePlotPost(COutPlotPost& output, unsigned int Ele)
{
	output << 4 << setw(11) << nodes_[0]->NodeNumber - 1 << nodes_[1]->NodeNumber - 1 << setw(9) << nodes_[2]->NodeNumber - 1 << setw(9)
		<< nodes_[3]->NodeNumber - 1 << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CSubpara::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 2; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For NEN_ node transition element, element stiffness is a 2NEN_x2NEN_ matrix, whose upper triangular part
//	has NEN_*(2*NEN_+1) elements
unsigned int CSubpara::SizeOfStiffnessMatrix() { return NEN_*(2*NEN_+1); }


//  Calculate shape function N at (xi, eta)
void CSubpara:: calN(double N[], double xi, double eta)
{
	double Ntemp[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
	Ntemp[0] = 0.25*(1. - xi)*(1. - eta);
	Ntemp[1] = 0.25*(1. + xi)*(1. - eta);
	Ntemp[2] = 0.25*(1. + xi)*(1. + eta);
	Ntemp[3] = 0.25*(1. - xi)*(1. + eta);

	Ntemp[8] = (1. - xi*xi)*(1. - eta*eta);
	for (unsigned int i = 0; i < 4; i++)
	{
		Ntemp[i] -= 0.25*Ntemp[8];
	}


	Ntemp[4] = 0.5*((1. - xi*xi)*(1. - eta)-Ntemp[8]);

	Ntemp[5] = 0.5*((1. - eta*eta)*(1. + xi)-Ntemp[8]);

	Ntemp[6] = 0.5*((1. - xi*xi)*(1. + eta)-Ntemp[8]);

	Ntemp[7] = 0.5*((1. - eta*eta)*(1. - xi)-Ntemp[8]);
	Ntemp[0] -= 0.5*(Ntemp[4] + Ntemp[7]);
	Ntemp[1] -= 0.5*(Ntemp[4] + Ntemp[5]);
	Ntemp[2] -= 0.5*(Ntemp[5] + Ntemp[6]);
	Ntemp[3] -= 0.5*(Ntemp[6] + Ntemp[7]);
	for(unsigned int i = 0; i < 9; i++ )
	{
		N[i] = Ntemp[i];
	}
}

//  Calculate invJ at (xi, eta) 
double CSubpara::detinvJ(double invJ[][2], double GB[], double xi, double eta)
{
	double J[2][2];  double DetJ;
	J[0][0] = 0.; J[0][1] = 0.; J[1][0] = 0.; J[1][1] = 0.; 
	calGB(GB, xi ,eta);
	for (unsigned i = 0; i < NEN_; i++ )
	{
		J[0][0] += GB[2*i]*nodes_[i]->XYZ[0]; J[0][1] += GB[2*i]*nodes_[i]->XYZ[1];
		J[1][0] += GB[2*i+1]*nodes_[i]->XYZ[0]; J[1][1] += GB[2*i+1]*nodes_[i]->XYZ[1];
	}
	DetJ = J[0][0]*J[1][1] - J[1][0]*J[0][1];
	invJ[0][0] = J[1][1]/DetJ;   invJ[0][1] = -J[0][1]/DetJ;
	invJ[1][0] = -J[1][0]/DetJ;   invJ[1][1] = J[0][0]/DetJ;
	return DetJ;
}

//  Calculate shape function N at (xi, eta)
void CSubpara:: calGB(double GB[], double xi, double eta)
{
	double Btemp[18] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	Btemp[0] = -0.25*(1. - eta); Btemp[1] = -0.25*(1. - xi); 
	Btemp[2] = 0.25*(1. - eta);  Btemp[3] = -0.25*(1. + xi);
	Btemp[4] = 0.25*(1. + eta);  Btemp[5] = 0.25*(1. + xi);
	Btemp[6] = -0.25*(1. + eta); Btemp[7] = 0.25*(1. - xi);


	Btemp[16] = -2.*xi*(1. - eta*eta); Btemp[17] = -2.*eta*(1. - xi*xi);
	for (unsigned int i = 0; i < 4; i++)
	{
		Btemp[2*i] += 0.5*xi*(1. - eta*eta); Btemp[2*i+1] += 0.5*eta*(1. - xi*xi);
	}



	Btemp[8] = 0.5*(-2.*xi*(1. - eta)-Btemp[16]);
	Btemp[9] = 0.5*(-(1. - xi*xi)-Btemp[17]);

	
	Btemp[10] = 0.5*((1. - eta*eta)-Btemp[16]);
	Btemp[11] = 0.5*(-2.*eta*(1. + xi)-Btemp[17]);

	
	Btemp[12] = 0.5*(-2.*xi*(1. + eta)-Btemp[16]);
	Btemp[13] = 0.5*((1. - xi*xi)-Btemp[17]);


	Btemp[14] = 0.5*(-(1 - eta*eta)-Btemp[16]);
	Btemp[15] = 0.5*(-2.*eta*(1. - xi)-Btemp[17]);


	Btemp[0] -= 0.5*(Btemp[8] + Btemp[14]); 
	Btemp[1] -= 0.5*(Btemp[9] + Btemp[15]);
	Btemp[2] -= 0.5*(Btemp[8] + Btemp[10]);
	Btemp[3] -= 0.5*(Btemp[9] + Btemp[11]);
	Btemp[4] -= 0.5*(Btemp[10] + Btemp[12]); 
	Btemp[5] -= 0.5*(Btemp[11] + Btemp[13]);
	Btemp[6] -= 0.5*(Btemp[12] + Btemp[14]);
	Btemp[7] -= 0.5*(Btemp[13] + Btemp[15]);

	for(unsigned int i = 0; i < 9; i++ )
	{
			GB[2*i] = Btemp[2*i];
			GB[2*i+1]   = Btemp[2*i+1];
	}
}

// Calculate B at (xi,eta)
void CSubpara:: calB(double B[], double xi, double eta)
{
	double invJ[2][2];
	double* GB = new double [2*NEN_];
	detJ = detinvJ(invJ, GB, xi, eta);
	for (unsigned int i = 0; i < NEN_; i ++)
	{
		B[2*i] = invJ[0][0]*GB[2*i] + invJ[0][1]*GB[2*i+1];
		B[2*i+1] = invJ[1][0]*GB[2*i] + invJ[1][1]*GB[2*i+1];
	}
//  delete [] GB;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CSubpara::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
// Gauss points and their weights respectively
	gspos[0] = -sqrt(3./5.); gspos[1] = 0.; gspos[2] = -gspos[0];
	w[0] = 5./9.; w[1] = 8./9.; w[2] = w[0];
	CSubparaMaterial* material_ = dynamic_cast<CSubparaMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double temp = material_->E/(1.-material_->poisson*material_->poisson);
	D[0] = temp; D[1] = material_->poisson*temp; D[2] = temp*(1.-material_->poisson)/2.;
	for (unsigned int ii = 0; ii < 3; ii ++)
	{
		for (unsigned int jj = 0; jj < 3; jj++ )
		{
			double* B = new double [2*NEN_];
			calB(B, gspos[ii], gspos[jj]);

			double* DB = new double [6*NEN_];
			double temp1 = w[ii]*w[jj]*detJ;
			for (unsigned int i = 0; i < NEN_; i++)
			{
				DB[6*i] = D[0]*B[2*i];
				DB[6*i+1] =  D[1]*B[2*i];
				DB[6*i+2] =  D[2]*B[2*i+1];
				DB[6*i+3] =  D[1]*B[2*i+1];
				DB[6*i+4] =  D[0]*B[2*i+1];
				DB[6*i+5] =  D[2]*B[2*i];
			}
			for (unsigned int i = 0; i < NEN_; i++)
			{
				Matrix[i*(2*i+1)] += (B[2*i]*DB[6*i]+B[2*i+1]*DB[6*i+2])*temp1;
				for (unsigned int j = 2*i+1; j < 2*NEN_; j++)
				{
					Matrix[j*(j+3)/2-2*i] += (B[2*i]*DB[3*j]+B[2*i+1]*DB[3*j+2])*temp1;
					Matrix[j*(j+3)/2-2*i-1] += (B[2*i]*DB[3*j+2]+B[2*i+1]*DB[3*j+1])*temp1;
				}
			}
			delete [] B; delete [] DB;
		}
	}

}

//	Calculate element stress 
void CSubpara::ElementStress(double* stress, double* Displacement)
{
	clear(stress, 45);
	CSubparaMaterial* material_ = dynamic_cast<CSubparaMaterial*>(ElementMaterial_);	// Pointer to material of the element
	for (unsigned int ii = 0; ii < 3; ii ++)
	{
		for (unsigned int jj = 0; jj < 3; jj ++)
		{
			double* B = new double [2*NEN_];
			calB(B, gspos[ii], gspos[jj]);
			double* DB = new double [6*NEN_];
			double temp1 = w[ii]*w[jj]*detJ;
			for (unsigned int i = 0; i < NEN_; i++)
			{
				DB[6*i] = D[0]*B[2*i];
				DB[6*i+1] =  D[1]*B[2*i];
				DB[6*i+2] =  D[2]*B[2*i+1];
				DB[6*i+3] =  D[1]*B[2*i+1];
				DB[6*i+4] =  D[0]*B[2*i+1];
				DB[6*i+5] =  D[2]*B[2*i];
			}
			unsigned int temp2 = 9*ii+3*jj;
			for (unsigned int i = 0; i < 2*NEN_; i++)
			{
				if (LocationMatrix_[i])
				{
					stress[temp2] += DB[3*i]*Displacement[LocationMatrix_[i]-1];
					stress[temp2+1] += DB[3*i+1]*Displacement[LocationMatrix_[i]-1];
					stress[temp2+2] += DB[3*i+2]*Displacement[LocationMatrix_[i]-1];
				}
			}
			double* N = new double [NEN_];
			calN(N, gspos[ii], gspos[jj]);
			unsigned int temp3 = 27 + 6*ii + 2*jj;
			for(unsigned int i = 0; i < NEN_; i++)
			{
				stress[temp3] += N[i]*nodes_[i]->XYZ[0];
				stress[temp3+1] += N[i]*nodes_[i]->XYZ[1];
			}
			delete [] B; delete [] DB; delete [] N;
		}
	}
}

void CSubpara::RecoverElementStress(double* Displacement, double* A)
{

}


// Caculate Gravity of Elements
void CSubpara::GravityCalculation(double* ptr_force)
{

}

//	Calculate element stress for plot
void CSubpara::ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* PostPositions)
{
	CSubparaMaterial* material_ = dynamic_cast<CSubparaMaterial*>(ElementMaterial_);	// Pointer to material of the element
	for(unsigned int j=0;j<NEN_;j++)
	{for (unsigned int i = 0; i < 2; i++)
	{
		PrePositions[i+3*j] = nodes_[j]->XYZ[i];
		if (LocationMatrix_[i+2*j])
			PostPositions[i+3*j]=nodes_[j]->XYZ[i]+Displacement[LocationMatrix_[i+2*j]-1];
		else
			PostPositions[i+3*j]=nodes_[j]->XYZ[i];

	}
	PrePositions[2 + 3 * j] = nodes_[j]->XYZ[2];
	PostPositions[2 + 3 * j] = nodes_[j]->XYZ[2];
	for (unsigned int i = 0; i < 6; i++)
	{
		stress[i+6*j]=nodes_[j]->stress_node[i];
	}
	}
}