#include "4Q.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
C4Q::C4Q()
{
	NEN_ = 4;    // Each element has 4 nodes
	nodes_=new CNode*[NEN_];

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
	unsigned int N1, N2,N3,N4;	// The first node from lower left, and the order is counterclock

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	return true;
}
//	Write element data to stream
void CBar::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber << setw(9) << nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
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
//	For 4 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 36 elements
unsigned int CBar::SizeOfStiffnessMatrix() { return 36; }
double * CBar::NMatrix(double eta,double psi)
{
	double * N;
	N[0] = 0.25*(1-psi)*(1-eta); 
	N[1] = 0.25*(1+psi)*(1-eta); 
	N[2] = 0.25*(1+psi)*(1+eta); 
	N[3] = 0.25*(1-psi)*(1+eta);
	return N;
}
void C4Q::ElementStiffness(double* Matrix){
	clear(Matrix, SizeOfStiffnessMatrix());
	C4QMaterial* material_ = dynamic_cast<C4QMaterial*>(ElementMaterial_);
	E1=material_->E;
	poisson1=material_->poisson;
	double  B[8];
	double  GN[8];
	double	J[4];
	double detj;
	double inv_j[4];
	double D[3];
	int ngp=6;
	double gp[ngp]=[-0.9324695142,-0.6612093865,-0.2386191861,0.2386191861,0.6612093865,0.9324695142];
    double w[ngp]=[0.1713244924,0.3607615730,0.4679139346,0.4679139346,0.3607615730,0.1713244924];
	D[0]=E1/(1-poisson1*poisson1);
	D[1]=poisson1*E1/(1-poisson1*poisson1);
	D[2]=E1/(2*(1+poisson1));
	for(unsigned int i=0,i++,i<ngp){ 
   		for(unsigned int j=0,j++,j<ngp){
    		eta = gp(i);             
    		psi = gp(j);
			GN[0]=0.25*(eta-1);
			GN[1]=0.25*(1-eta);
			GN[2]=0.25*(eta+1);
			GN[3]=0.25*(-eta-1);
			GN[4]=0.25*(psi-1);
			GN[5]=0.25*(-psi-1);
			GN[6]=0.25*(psi+1);
			GN[7]=0.25*(1-psi);
			for (unsigned int i=0; i < 4; i++){
				J[i]=0;
			}
			for (unsigned int j = 0; j < 4; j++){
						J[0]+=GN[j]*nodes_[j]->XYZ[0];
						J[1]+=GN[j]*nodes_[j]->XYZ[1];
						J[2]+=GN[j + 4]*nodes_[j]->XYZ[0];
						J[3]+=GN[j + 4]*nodes_[j]->XYZ[1];
			}
			detj=J[0]*J[3]-J[1]*J[2];
			inv_j[0]=J[3]/detJ;
			inv_j[1]=-J[1]/detJ;
			inv_j[2]=-J[2]/detJ;
			inv_j[3]=J[0]/detJ;
			for (unsigned int i=0; i < 2; i++){
				for (unsigned int j = 0; j < 4; j++){
					B[4*i+j]=inv_j[2*i]*GN[j]+inv_j[2*i+1]*GN[4+j];
				}
			}
			Matrix[0]+=w[i]*w[j]*(B[0]*D[0]*B[0] + B[4]*D[2]*B[4])*fabs(detj);
			Matrix[1]+=w[i]*w[j]*(B[0]*D[0]*B[0] + B[4]*D[2]*B[4])*fabs(detj);
			Matrix[2]+=w[i]*w[j]*(B[0]*D[1]*B[4] + B[4]*D[2]*B[0])*fabs(detj);
			Matrix[3]+=w[i]*w[j]*(B[0]*D[1]*B[5] + B[4]*D[2]*B[1])*fabs(detj);
			Matrix[4]+=w[i]*w[j]*(B[0]*D[0]*B[3] + B[4]*D[2]*B[7])*fabs(detj);
			Matrix[5]+=w[i]*w[j]*(B[4]*D[1]*B[1] + B[0]*D[2]*B[5])*fabs(detj);
			Matrix[6]+=w[i]*w[j]*(B[0]*D[2]*B[3] + B[4]*D[0]*B[7])*fabs(detj);
			Matrix[7]+=w[i]*w[j]*(B[1]*D[1]*B[6] + B[5]*D[2]*B[2])*fabs(detj);
			Matrix[8]+=w[i]*w[j]*(B[0]*D[1]*B[4] + B[4]*D[2]*B[0])*fabs(detj);
			Matrix[9]+=w[i]*w[j]*(B[0]*D[0]*B[1] + B[4]*D[2]*B[5])*fabs(detj);
			Matrix[10]+=w[i]*w[j]*(B[0]*D[0]*B[2] + B[4]*D[2]*B[6])*fabs(detj);
			Matrix[11]+=w[i]*w[j]*(B[0]*D[1]*B[7] + B[4]*D[2]*B[3])*fabs(detj);
			Matrix[12]+=w[i]*w[j]*(B[0]*D[2]*B[1] + B[4]*D[0]*B[5])*fabs(detj);
			Matrix[13]+=w[i]*w[j]*(B[1]*D[0]*B[0] + B[5]*D[2]*B[4])*fabs(detj);
			Matrix[14]+=w[i]*w[j]*(B[1]*D[0]*B[3] + B[5]*D[2]*B[7])*fabs(detj);
			Matrix[15]+=w[i]*w[j]*(B[0]*D[1]*B[5] + B[4]*D[2]*B[1])*fabs(detj);
			Matrix[16]+=w[i]*w[j]*(B[0]*D[1]*B[6] + B[4]*D[2]*B[2])*fabs(detj);
			Matrix[17]+=w[i]*w[j]*(B[4]*D[1]*B[0] + B[0]*D[2]*B[4])*fabs(detj);
			Matrix[18]+=w[i]*w[j]*(B[4]*D[1]*B[2] + B[0]*D[2]*B[6])*fabs(detj);
			Matrix[19]+=w[i]*w[j]*(B[1]*D[1]*B[4] + B[5]*D[2]*B[0])*fabs(detj);
			Matrix[20]+=w[i]*w[j]*(B[1]*D[1]*B[7] + B[5]*D[2]*B[3])*fabs(detj);
			Matrix[21]+=w[i]*w[j]*(B[0]*D[0]*B[3] + B[4]*D[2]*B[7])*fabs(detj);
			Matrix[22]+=w[i]*w[j]*(B[0]*D[2]*B[0] + B[4]*D[0]*B[4])*fabs(detj);
			Matrix[23]+=w[i]*w[j]*(B[0]*D[2]*B[2] + B[4]*D[0]*B[6])*fabs(detj);
			Matrix[24]+=w[i]*w[j]*(B[1]*D[0]*B[1] + B[5]*D[2]*B[5])*fabs(detj);
			Matrix[25]+=w[i]*w[j]*(B[5]*D[1]*B[0] + B[1]*D[2]*B[4])*fabs(detj);
			Matrix[26]+=w[i]*w[j]*(B[4]*D[1]*B[1] + B[0]*D[2]*B[5])*fabs(detj);
			Matrix[27]+=w[i]*w[j]*(B[4]*D[1]*B[3] + B[0]*D[2]*B[7])*fabs(detj);
			Matrix[28]+=w[i]*w[j]*(B[1]*D[1]*B[5] + B[5]*D[2]*B[1])*fabs(detj);
			Matrix[29]+=w[i]*w[j]*(B[1]*D[2]*B[0] + B[5]*D[0]*B[4])*fabs(detj);
			Matrix[30]+=w[i]*w[j]*(B[0]*D[2]*B[3] + B[4]*D[0]*B[7])*fabs(detj);
			Matrix[31]+=w[i]*w[j]*(B[1]*D[0]*B[2] + B[5]*D[2]*B[6])*fabs(detj);
			Matrix[32]+=w[i]*w[j]*(B[5]*D[1]*B[1] + B[1]*D[2]*B[5])*fabs(detj);
			Matrix[33]+=w[i]*w[j]*(B[1]*D[1]*B[6] + B[5]*D[2]*B[2])*fabs(detj);
			Matrix[34]+=w[i]*w[j]*(B[1]*D[2]*B[1] + B[5]*D[0]*B[5])*fabs(detj);
			Matrix[35]+=w[i]*w[j]*(B[5]*D[1]*B[2] + B[1]*D[2]*B[6])*fabs(detj);
		}
	}

}
