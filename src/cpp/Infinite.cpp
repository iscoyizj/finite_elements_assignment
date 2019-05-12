/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Infinite.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CInfi::CInfi()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 8;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ =  nullptr;
}

//	Desconstructor
CInfi::~CInfi()
{
}

//	Read element data from stream Input
bool CInfi::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2 ,N3 ,N4;	// 4 node numbers which is counterclock

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = dynamic_cast<CInfiMaterial*>(MaterialSets)+MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CInfi::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber 
		   << setw(9) << nodes_[1]->NodeNumber <<setw(9) << nodes_[2]->NodeNumber <<setw(9) << nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

void CInfi::WritePlot(COutPlot& output, unsigned int Ele)
{
	output << 4 << setw(11) << nodes_[0]->NodeNumber - 1 << nodes_[1]->NodeNumber - 1 << setw(9) << nodes_[2]->NodeNumber - 1 << setw(9)
		<< nodes_[3]->NodeNumber - 1 << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CInfi::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 2; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 4 node Infinite element, element stiffness is a 8x8 matrix, whose upper triangular part
//	has 36 elements
unsigned int CInfi::SizeOfStiffnessMatrix() { return 36; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CInfi::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	CInfiMaterial* material_ = dynamic_cast<CInfiMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double E1;
	double poi;
	double D1;
	if(material_->etype==1)
//	Calculate D1
	{   E1=material_->E;
	    poi=material_->poisson;
		D1 = E1/(1-poi*poi);}
	else if(material_->etype==2)
		{   E1=material_->E/(1-material_->poisson*material_->poisson);
	    poi=material_->poisson/(1-material_->poisson);
		D1 = E1/(1-poi*poi);}
	else if(material_->etype==3)
		{   E1=material_->E/(1-material_->poisson*material_->poisson);
	   poi=material_->poisson/(1-material_->poisson);
	  D1 = E1/(1-poi*poi);}

	if (material_->etype==1||material_->etype==2)
{//	Calculate element Jacobian matrix
	double J[4];
	double GN[8];
	double detJ;
	double JJ[4];
	double B[8];
	double M[36];
	for(unsigned m=0;m<36;m++)
		M[m]=0;

	for (unsigned i=0;i<4;)
		{double psi = 1/sqrt(3);
	     double eta = 1/sqrt(3);
	    if(i==0)
		{psi=psi;eta=eta;}
		else if(i==1)
		{psi=psi;eta=-eta;}
		else if(i==2)
		{psi=-psi;eta=eta;}
		else if(i==3)
		{psi=-psi;eta=-eta;}
	
	GN[0]=0.5*(1-eta)*(-2)/((1-psi)*(1-psi));
	GN[1]=0.5*(1-eta)*(2)/((1-psi)*(1-psi));
	GN[2]=0.5*(1+eta)*(2)/((1-psi)*(1-psi));
	GN[3]=0.5*(1+eta)*(-2)/((1-psi)*(1-psi));
	GN[4]=-0.5*(-2)*psi/(1-psi);
	GN[5]=-0.5*(1+(2)*psi/(1-psi));
	GN[6]=0.5*(1+(2)*psi/(1-psi));
	GN[7]=0.5*(-2)*psi/(1-psi);
	
	J[0]=0;J[1]=0;J[2]=0;J[3]=0;
	for (unsigned int j = 0; j < 4; j++)
	{J[0]+=GN[j]*nodes_[j]->XYZ[0];
	J[1]+=GN[j]*nodes_[j]->XYZ[1];
	J[2]+=GN[j+4]*nodes_[j]->XYZ[0];
	J[3]+=GN[j+4]*nodes_[j]->XYZ[1];
	}
	
	detJ=J[0]*J[3]-J[1]*J[2];
	
	JJ[0]=J[3]/detJ;
	JJ[1]=-J[1]/detJ;
	JJ[2]=-J[2]/detJ;
	JJ[3]=J[0]/detJ;
	
	B[0]=JJ[0]*GN[0]+JJ[1]*GN[4];
	B[1]=JJ[0]*GN[1]+JJ[1]*GN[5];
	B[2]=JJ[0]*GN[2]+JJ[1]*GN[6];
	B[3]=JJ[0]*GN[3]+JJ[1]*GN[7];
	B[4]=JJ[2]*GN[0]+JJ[3]*GN[4];
	B[5]=JJ[2]*GN[1]+JJ[3]*GN[5];
	B[6]=JJ[2]*GN[2]+JJ[3]*GN[6];
	B[7]=JJ[2]*GN[3]+JJ[3]*GN[7];
	

	M[0] += D1*(B[0]*B[0]+B[4]*B[4]*(1-poi)/2)*fabs(detJ);
	M[1] += D1*(B[4]*B[4]+B[0]*B[0]*(1-poi)/2)*fabs(detJ);
	M[2] += D1*(B[0]*B[4]*poi+B[4]*B[0]*(1-poi)/2)*fabs(detJ);
	M[3] += D1*(B[1]*B[1]+B[5]*B[5]*(1-poi)/2)*fabs(detJ);
	M[4] += D1*(B[4]*B[1]*poi+B[0]*B[5]*(1-poi)/2)*fabs(detJ);
	M[5] += D1*(B[0]*B[1]+B[4]*B[5]*(1-poi)/2)*fabs(detJ);
	M[6] += D1*(B[5]*B[5]+B[1]*B[1]*(1-poi)/2)*fabs(detJ);
	M[7] += D1*(B[1]*B[5]*poi+B[5]*B[1]*(1-poi)/2)*fabs(detJ);
	M[8] += D1*(B[4]*B[5]+B[0]*B[1]*(1-poi)/2)*fabs(detJ);
	M[9] += D1*(B[0]*B[5]*poi+B[4]*B[1]*(1-poi)/2)*fabs(detJ);
	M[10] += D1*(B[2]*B[2]+B[6]*B[6]*(1-poi)/2)*fabs(detJ);
	M[11] += D1*(B[5]*B[2]*poi+B[1]*B[6]*(1-poi)/2)*fabs(detJ);
	M[12] += D1*(B[1]*B[2]+B[5]*B[6]*(1-poi)/2)*fabs(detJ);
	M[13] += D1*(B[4]*B[2]*poi+B[0]*B[6]*(1-poi)/2)*fabs(detJ);
	M[14] += D1*(B[0]*B[2]+B[4]*B[6]*(1-poi)/2)*fabs(detJ);
	M[15] += D1*(B[6]*B[6]+B[2]*B[2]*(1-poi)/2)*fabs(detJ);
	M[16] += D1*(B[2]*B[6]*poi+B[6]*B[2]*(1-poi)/2)*fabs(detJ);
	M[17] += D1*(B[5]*B[6]+B[1]*B[2]*(1-poi)/2)*fabs(detJ);
	M[18] += D1*(B[1]*B[6]*poi+B[5]*B[2]*(1-poi)/2)*fabs(detJ);
	M[19] += D1*(B[4]*B[6]+B[0]*B[2]*(1-poi)/2)*fabs(detJ);
	M[20] += D1*(B[0]*B[6]*poi+B[4]*B[2]*(1-poi)/2)*fabs(detJ);
	M[21] += D1*(B[3]*B[3]+B[7]*B[7]*(1-poi)/2)*fabs(detJ);
	M[22] += D1*(B[6]*B[3]*poi+B[2]*B[7]*(1-poi)/2)*fabs(detJ);
	M[23] += D1*(B[2]*B[3]+B[6]*B[7]*(1-poi)/2)*fabs(detJ);
	M[24] += D1*(B[5]*B[3]*poi+B[1]*B[7]*(1-poi)/2)*fabs(detJ);
	M[25] += D1*(B[1]*B[3]+B[5]*B[7]*(1-poi)/2)*fabs(detJ);
	M[26] += D1*(B[4]*B[3]*poi+B[0]*B[7]*(1-poi)/2)*fabs(detJ);
	M[27] += D1*(B[0]*B[3]+B[4]*B[7]*(1-poi)/2)*fabs(detJ);
	M[28] += D1*(B[7]*B[7]+B[3]*B[3]*(1-poi)/2)*fabs(detJ);
	M[29] += D1*(B[3]*B[7]*poi+B[7]*B[3]*(1-poi)/2)*fabs(detJ);
	M[30] += D1*(B[6]*B[7]+B[2]*B[3]*(1-poi)/2)*fabs(detJ);
	M[31] += D1*(B[2]*B[7]*poi+B[6]*B[3]*(1-poi)/2)*fabs(detJ);
	M[32] += D1*(B[5]*B[7]+B[1]*B[3]*(1-poi)/2)*fabs(detJ);
	M[33] += D1*(B[1]*B[7]*poi+B[5]*B[3]*(1-poi)/2)*fabs(detJ);
	M[34] += D1*(B[4]*B[7]+B[0]*B[3]*(1-poi)/2)*fabs(detJ);
	M[35] += D1*(B[0]*B[7]*poi+B[4]*B[3]*(1-poi)/2)*fabs(detJ);
	i=i+1;
	};
		Matrix[0]=M[0];
		Matrix[1]=M[1];
		Matrix[2]=M[2];
		Matrix[3]=M[3];
		Matrix[4]=M[4];
		Matrix[5]=M[5];
		Matrix[6]=M[6];
		Matrix[7]=M[7];
		Matrix[8]=M[8];
		Matrix[9]=M[9];
		Matrix[10]=M[10];
		Matrix[11]=M[11];
		Matrix[12]=M[12];
		Matrix[13]=M[13];
		Matrix[14]=M[14];
		Matrix[15]=M[15];
		Matrix[16]=M[16];
		Matrix[17]=M[17];
		Matrix[18]=M[18];
		Matrix[19]=M[19];
		Matrix[20]=M[20];
		Matrix[21]=M[21];
		Matrix[22]=M[22];
		Matrix[23]=M[23];
		Matrix[24]=M[24];
		Matrix[25]=M[25];
		Matrix[26]=M[26];
		Matrix[27]=M[27];
		Matrix[28]=M[28];
		Matrix[29]=M[29];
		Matrix[30]=M[30];
		Matrix[31]=M[31];
		Matrix[32]=M[32];
		Matrix[33]=M[33];
		Matrix[34]=M[34];
		Matrix[35]=M[35];
}
else if(material_->etype==3)
   {//	Calculate element Jacobian matrix
	double J[4];
	double GN[8];
	double detJ;
	double JJ[4];
	double B[8];
	double M[36];
	double N[4];
	double x;
	double Nr[4];
	for(unsigned m=0;m<36;m++)
		M[m]=0;

	for (unsigned i=0;i<4;)
		{double psi = 1/sqrt(3);
	     double eta = 1/sqrt(3);
	    if(i==0)
		{psi=psi;eta=eta;}
		else if(i==1)
		{psi=psi;eta=-eta;}
		else if(i==2)
		{psi=-psi;eta=eta;}
		else if(i==3)
		{psi=-psi;eta=-eta;}
	    x=0;
		N[0]=0.25*(1-psi)*(1-eta);
		N[2]=0.25*(1+psi)*(1+eta);
		N[1]=0.25*(1+psi)*(1-eta);
		N[3]=0.25*(1-psi)*(1+eta);
		for(unsigned int p=0;p<4;p++)
		{x+=N[p]*nodes_[p]->XYZ[0];
		}
		for(unsigned int q=0;q<4;q++)
		{Nr[q]=N[q]/x;
		}

	GN[0]=0.25*(eta-1);
	GN[1]=0.25*(1-eta);
	GN[2]=0.25*(eta+1);
	GN[3]=0.25*(-eta-1);
	GN[4]=0.25*(psi-1);
	GN[5]=0.25*(-psi-1);
	GN[6]=0.25*(psi+1);
	GN[7]=0.25*(1-psi);
	
	J[0]=0;J[1]=0;J[2]=0;J[3]=0;
	for (unsigned int j = 0; j < 4; j++)
	{J[0]+=GN[j]*nodes_[j]->XYZ[0];
	J[1]+=GN[j]*nodes_[j]->XYZ[1];
	J[2]+=GN[j+4]*nodes_[j]->XYZ[0];
	J[3]+=GN[j+4]*nodes_[j]->XYZ[1];
	}
	
	detJ=J[0]*J[3]-J[1]*J[2];
	
	JJ[0]=J[3]/detJ;
	JJ[1]=-J[1]/detJ;
	JJ[2]=-J[2]/detJ;
	JJ[3]=J[0]/detJ;
	
	B[0]=JJ[0]*GN[0]+JJ[1]*GN[4];
	B[1]=JJ[0]*GN[1]+JJ[1]*GN[5];
	B[2]=JJ[0]*GN[2]+JJ[1]*GN[6];
	B[3]=JJ[0]*GN[3]+JJ[1]*GN[7];
	B[4]=JJ[2]*GN[0]+JJ[3]*GN[4];
	B[5]=JJ[2]*GN[1]+JJ[3]*GN[5];
	B[6]=JJ[2]*GN[2]+JJ[3]*GN[6];
	B[7]=JJ[2]*GN[3]+JJ[3]*GN[7];
	

	M[0] += D1*(B[0]*B[0]+B[4]*B[4]*(1-poi)/2+2*Nr[0]*B[0]*poi+Nr[0]*Nr[0])*fabs(detJ)*x*2*3.1415;
	M[1] += D1*(B[4]*B[4]+B[0]*B[0]*(1-poi)/2)*fabs(detJ)*x*2*3.1415;
	M[2] += D1*(B[0]*B[4]*poi+B[4]*B[0]*(1-poi)/2+B[4]*Nr[0]*poi)*fabs(detJ)*x*2*3.1415;
	M[3] += D1*(B[1]*B[1]+B[5]*B[5]*(1-poi)/2+2*B[1]*Nr[1]*poi+Nr[1]*Nr[1])*fabs(detJ)*x*2*3.1415;
	M[4] += D1*(B[4]*B[1]*poi+B[0]*B[5]*(1-poi)/2+B[4]*poi*Nr[1])*fabs(detJ)*x*2*3.1415;
	M[5] += D1*(B[0]*B[1]+B[4]*B[5]*(1-poi)/2+B[1]*Nr[0]*poi+B[0]*Nr[1]*poi+Nr[0]*Nr[1])*fabs(detJ)*x*2*3.1415;
	M[6] += D1*(B[5]*B[5]+B[1]*B[1]*(1-poi)/2)*fabs(detJ)*x*2*3.1415;
	M[7] += D1*(B[1]*B[5]*poi+B[5]*B[1]*(1-poi)/2+Nr[1]*poi*B[5])*fabs(detJ)*x*2*3.1415;
	M[8] += D1*(B[4]*B[5]+B[0]*B[1]*(1-poi)/2)*fabs(detJ)*x*2*3.1415;
	M[9] += D1*(B[0]*B[5]*poi+B[4]*B[1]*(1-poi)/2+Nr[0]*poi*B[5])*fabs(detJ)*x*2*3.1415;
	M[10] += D1*(B[2]*B[2]+B[6]*B[6]*(1-poi)/2+Nr[2]*poi*B[2]*2+Nr[2]*Nr[2])*fabs(detJ)*x*2*3.1415;
	M[11] += D1*(B[5]*B[2]*poi+B[1]*B[6]*(1-poi)/2+B[5]*poi*Nr[2])*fabs(detJ)*x*2*3.1415;
	M[12] += D1*(B[1]*B[2]+B[5]*B[6]*(1-poi)/2+Nr[1]*poi*B[2]+Nr[2]*poi*B[1]+Nr[1]*Nr[2])*fabs(detJ)*x*2*3.1415;
	M[13] += D1*(B[4]*B[2]*poi+B[0]*B[6]*(1-poi)/2+Nr[2]*B[4]*poi)*fabs(detJ)*x*2*3.1415;
	M[14] += D1*(B[0]*B[2]+B[4]*B[6]*(1-poi)/2+Nr[0]*poi*B[2]+Nr[2]*B[0]*poi+Nr[0]*Nr[2])*fabs(detJ)*x*2*3.1415;
	M[15] += D1*(B[6]*B[6]+B[2]*B[2]*(1-poi)/2)*fabs(detJ)*x*2*3.1415;
	M[16] += D1*(B[2]*B[6]*poi+B[6]*B[2]*(1-poi)/2+Nr[2]*poi*B[6])*fabs(detJ)*x*2*3.1415;
	M[17] += D1*(B[5]*B[6]+B[1]*B[2]*(1-poi)/2)*fabs(detJ)*x*2*3.1415;
	M[18] += D1*(B[1]*B[6]*poi+B[5]*B[2]*(1-poi)/2+Nr[1]*poi*B[6])*fabs(detJ)*x*2*3.1415;
	M[19] += D1*(B[4]*B[6]+B[0]*B[2]*(1-poi)/2)*fabs(detJ)*x*2*3.1415;
	M[20] += D1*(B[0]*B[6]*poi+B[4]*B[2]*(1-poi)/2+Nr[0]*poi*B[6])*fabs(detJ)*x*2*3.1415;
	M[21] += D1*(B[3]*B[3]+B[7]*B[7]*(1-poi)/2+Nr[3]*poi*B[3]*2+Nr[3]*Nr[3])*fabs(detJ)*x*2*3.1415;
	M[22] += D1*(B[6]*B[3]*poi+B[2]*B[7]*(1-poi)/2+B[6]*poi*Nr[3])*fabs(detJ)*x*2*3.1415;
	M[23] += D1*(B[2]*B[3]+B[6]*B[7]*(1-poi)/2+Nr[2]*poi*B[3]+Nr[3]*poi*B[2]+Nr[2]*Nr[3])*fabs(detJ)*x*2*3.1415;
	M[24] += D1*(B[5]*B[3]*poi+B[1]*B[7]*(1-poi)/2+B[5]*poi*Nr[3])*fabs(detJ)*x*2*3.1415;
	M[25] += D1*(B[1]*B[3]+B[5]*B[7]*(1-poi)/2+Nr[1]*poi*B[3]+Nr[3]*poi*B[1]+Nr[1]*Nr[3])*fabs(detJ)*x*2*3.1415;
	M[26] += D1*(B[4]*B[3]*poi+B[0]*B[7]*(1-poi)/2+B[4]*poi*Nr[3])*fabs(detJ)*x*2*3.1415;
	M[27] += D1*(B[0]*B[3]+B[4]*B[7]*(1-poi)/2+Nr[0]*poi*B[3]+Nr[3]*poi*B[0]+Nr[0]*Nr[3])*fabs(detJ)*x*2*3.1415;
	M[28] += D1*(B[7]*B[7]+B[3]*B[3]*(1-poi)/2)*fabs(detJ)*x*2*3.1415;
	M[29] += D1*(B[3]*B[7]*poi+B[7]*B[3]*(1-poi)/2+B[7]*Nr[3]*poi)*fabs(detJ)*x*2*3.1415;
	M[30] += D1*(B[6]*B[7]+B[2]*B[3]*(1-poi)/2)*fabs(detJ)*x*2*3.1415;
	M[31] += D1*(B[2]*B[7]*poi+B[6]*B[3]*(1-poi)/2+B[7]*Nr[2]*poi)*fabs(detJ)*x*2*3.1415;
	M[32] += D1*(B[5]*B[7]+B[1]*B[3]*(1-poi)/2)*fabs(detJ)*x*2*3.1415;
	M[33] += D1*(B[1]*B[7]*poi+B[5]*B[3]*(1-poi)/2+B[7]*Nr[1]*poi)*fabs(detJ)*x*2*3.1415;
	M[34] += D1*(B[4]*B[7]+B[0]*B[3]*(1-poi)/2)*fabs(detJ)*x*2*3.1415;
	M[35] += D1*(B[0]*B[7]*poi+B[4]*B[3]*(1-poi)/2+B[7]*Nr[0]*poi)*fabs(detJ)*x*2*3.1415;
	i=i+1;
	};
		Matrix[0]=M[0];
		Matrix[1]=M[1];
		Matrix[2]=M[2];
		Matrix[3]=M[3];
		Matrix[4]=M[4];
		Matrix[5]=M[5];
		Matrix[6]=M[6];
		Matrix[7]=M[7];
		Matrix[8]=M[8];
		Matrix[9]=M[9];
		Matrix[10]=M[10];
		Matrix[11]=M[11];
		Matrix[12]=M[12];
		Matrix[13]=M[13];
		Matrix[14]=M[14];
		Matrix[15]=M[15];
		Matrix[16]=M[16];
		Matrix[17]=M[17];
		Matrix[18]=M[18];
		Matrix[19]=M[19];
		Matrix[20]=M[20];
		Matrix[21]=M[21];
		Matrix[22]=M[22];
		Matrix[23]=M[23];
		Matrix[24]=M[24];
		Matrix[25]=M[25];
		Matrix[26]=M[26];
		Matrix[27]=M[27];
		Matrix[28]=M[28];
		Matrix[29]=M[29];
		Matrix[30]=M[30];
		Matrix[31]=M[31];
		Matrix[32]=M[32];
		Matrix[33]=M[33];
		Matrix[34]=M[34];
		Matrix[35]=M[35];};





}


//	Calculate element stress 
void CInfi::ElementStress(double* stress, double* Displacement)
{
	CInfiMaterial* material_ = dynamic_cast<CInfiMaterial*>(ElementMaterial_);	// Pointer to material of the element
	//	Calculate D1
	double D1 = material_->E/(1-material_->poisson*material_->poisson);
	double eta;
	double psi;
	for(unsigned j=0;j<4;j++)
	{if(j==0)
		{psi = 1/sqrt(3);
	eta = 1/sqrt(3);}
	else if (j==1)
		{ psi = 1/sqrt(3);
	 eta = -1/sqrt(3);}
	else if (j==2)
	{ psi = -1/sqrt(3);
	 eta = 1/sqrt(3);}
	else if(j==3)
    { psi = -1/sqrt(3);
	 eta = -1/sqrt(3);}
	double GN1[8];
		GN1[0]=0.5*(1-eta)*(-2)/((1-psi)*(1-psi));
	GN1[1]=0.5*(1-eta)*(2)/((1-psi)*(1-psi));
	GN1[2]=0.5*(1+eta)*(2)/((1-psi)*(1-psi));
	GN1[3]=0.5*(1+eta)*(-2)/((1-psi)*(1-psi));
	GN1[4]=-0.5*(-2)*psi/(1-psi);
	GN1[5]=-0.5*(1+(2)*psi/(1-psi));
	GN1[6]=0.5*(1+(2)*psi/(1-psi));
	GN1[7]=0.5*(-2)*psi/(1-psi);
	double N1[4];
	N1[0]=(-2)*psi/(1-psi)*(0.5-0.5*eta);
	N1[1]=(1+(2)*psi/(1-psi))*(0.5-0.5*eta);
	N1[2]=(1+(2)*psi/(1-psi))*(0.5+0.5*eta);
	N1[3]=(-2)*psi/(1-psi)*(0.5+0.5*eta);
	double J1[4];
	J1[0]=0;
	J1[1]=0;
	J1[2]=0;
	J1[3]=0;
	for (unsigned int i = 0; i < 4; i++)
	{J1[0]=J1[0]+GN1[i]*nodes_[i]->XYZ[0];
	J1[1]=J1[1]+GN1[i]*nodes_[i]->XYZ[1];
	J1[2]=J1[2]+GN1[i+4]*nodes_[i]->XYZ[0];
	J1[3]=J1[3]+GN1[i+4]*nodes_[i]->XYZ[1];
	}
	double detJ1=J1[0]*J1[3]-J1[1]*J1[2];
	double JJ1[4];
	JJ1[0]=J1[3]/detJ1;
	JJ1[1]=-J1[1]/detJ1;
	JJ1[2]=-J1[2]/detJ1;
	JJ1[3]=J1[0]/detJ1;
	double B1[8];
	for (unsigned int i = 0; i < 4; i++)
		{B1[i]=JJ1[0]*GN1[i]+JJ1[1]*GN1[i+4];
	     B1[i+4]=JJ1[2]*GN1[i]+JJ1[3]*GN1[i+4];
	     }
	for(unsigned int i = 0; i < 4; i++)
	    {stress[0+6*j] += N1[i]*nodes_[i]->XYZ[0];
	     stress[1+6*j] += N1[i]*nodes_[i]->XYZ[1];
		 stress[2+6*j] += N1[i]*nodes_[i]->XYZ[2];}
	double d[8];
	for(unsigned int i = 0; i < 8; i++)
	{if (LocationMatrix_[i])
	    d[i]= Displacement[LocationMatrix_[i]-1];
	else if(LocationMatrix_[i]==0)
		d[i]=0;
	}
	double s[3];
	s[0]=B1[0]*d[0]+B1[1]*d[2]+B1[2]*d[4]+B1[3]*d[6];
	s[1]=B1[4]*d[1]+B1[5]*d[3]+B1[6]*d[5]+B1[7]*d[7];
	s[2]=B1[4]*d[0]+B1[0]*d[1]+B1[5]*d[2]+B1[1]*d[3]+B1[6]*d[4]+B1[2]*d[5]+B1[7]*d[6]+B1[3]*d[7];
	stress[3+6*j]+=D1*(s[0]+material_->poisson*s[1]);
	stress[4+6*j]+=D1*(s[0]*material_->poisson+s[1]);
	stress[5+6*j]+=D1*(1-material_->poisson)/2*s[2];
	}
}
//	Calculate element stress for plot
void CInfi::ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* PostPositions)
{
	CInfiMaterial* material_ = dynamic_cast<CInfiMaterial*>(ElementMaterial_);	// Pointer to material of the element
	for(unsigned int j=0;j<NEN_;j++)
	{
		for (unsigned int i = 0; i < 2; i++)
	    {
			PrePositions[i + 3 * j] = nodes_[j]->XYZ[i];
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

void CInfi::RecoverElementStress(double* Displacement, double* A)
{

}

// Caculate Gravity of Elements
void CInfi::GravityCalculation(double* ptr_force)
{

}