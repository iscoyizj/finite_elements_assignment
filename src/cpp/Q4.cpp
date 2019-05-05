#include "Q4.h"
#include<iostream>
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

CQ4::~CQ4(){}

bool CQ4::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList){
	unsigned int N;
	Input>>N;// element number
	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl 
			 << "    Expected element : " << Ele + 1 << endl
			 << "    Provided element : " << N << endl;

		return false;
	}
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// Node number in counter clockwise order
	Input >> N1 >> N2 >> N3 >> N4>> MSet;
	ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	return true;
}

void CQ4::Write(COutputter& output, unsigned int Ele){
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber <<setw(9)<<nodes_[2]->NodeNumber<<setw(9)<<nodes_[3]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

void CQ4::GenerateLocationMatrix(){
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 2; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

unsigned int CQ4::SizeOfStiffnessMatrix() { return 36; }

void CQ4::ElementStiffness(double* Matrix){
	clear(Matrix, SizeOfStiffnessMatrix());
	double X[4]={nodes_[0]->XYZ[0],nodes_[1]->XYZ[0],nodes_[2]->XYZ[0],nodes_[3]->XYZ[0]};
	double Y[4]={nodes_[0]->XYZ[1],nodes_[1]->XYZ[1],nodes_[2]->XYZ[1],nodes_[3]->XYZ[1]};
	double gausspoint[2]={-0.57735027, 0.57735027};	//Gauss point when ngp=2
	double weight[2]={1,1};	//Gauss weight

	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);
	double E_=material_->E;
	double nu_=material_->nu;
	double D1=E_/(1.0-nu_*nu_);
	double D2=E_*nu_/(1.0-nu_*nu_);
	double D3=E_/2.0/(1.0+nu_);
	//D=[D1,D2,0;D2,D1,0;0,0,D3]

	for(unsigned int i=0;i<2;i++){
		double eta=gausspoint[i];
		for(unsigned int j=0;j<2;j++){
			double psi=gausspoint[j];
	double GN[2][4]={(eta-1)/4,(1-eta)/4,(1+eta)/4,(-eta-1)/4,(psi-1)/4,(-psi-1)/4,(1+psi)/4,(1-psi)/4};
	double Jocbian[2][2]={GN[0][0]*X[0]+GN[0][1]*X[1]+GN[0][2]*X[2]+GN[0][3]*X[3],GN[0][0]*Y[0]+GN[0][1]*Y[1]+GN[0][2]*Y[2]+GN[0][3]*Y[3],
	GN[1][0]*X[0]+GN[1][1]*X[1]+GN[1][2]*X[2]+GN[1][3]*X[3],GN[1][0]*Y[0]+GN[1][1]*Y[1]+GN[1][2]*Y[2]+GN[1][3]*Y[3]};
	double detJocbian=Jocbian[0][0]*Jocbian[1][1]-Jocbian[1][0]*Jocbian[0][1];
	double Jocbian_inv[2][2]={Jocbian[1][1]/detJocbian,-Jocbian[0][1]/detJocbian,-Jocbian[1][0]/detJocbian,Jocbian[0][0]/detJocbian};
	double B1x=Jocbian_inv[0][0]*GN[0][0]+Jocbian_inv[0][1]*GN[1][0];
	double B2x=Jocbian_inv[0][0]*GN[0][1]+Jocbian_inv[0][1]*GN[1][1];
	double B3x=Jocbian_inv[0][0]*GN[0][2]+Jocbian_inv[0][1]*GN[1][2];
	double B4x=Jocbian_inv[0][0]*GN[0][3]+Jocbian_inv[0][1]*GN[1][3];
	double B1y=Jocbian_inv[1][0]*GN[0][0]+Jocbian_inv[1][1]*GN[1][0];
	double B2y=Jocbian_inv[1][0]*GN[0][1]+Jocbian_inv[1][1]*GN[1][1];
	double B3y=Jocbian_inv[1][0]*GN[0][2]+Jocbian_inv[1][1]*GN[1][2];
	double B4y=Jocbian_inv[1][0]*GN[0][3]+Jocbian_inv[1][1]*GN[1][3];
	double WiWjdetJ=fabs(detJocbian)*weight[i]*weight[j];
	Matrix[0]+=(D1*B1x*B1x + D3*B1y*B1y)*WiWjdetJ;
	Matrix[1]+=(D3*B1x*B1x + D1*B1y*B1y)*WiWjdetJ;
	Matrix[2]+=(B1x*B1y*D2 + B1x*B1y*D3)*WiWjdetJ;
	Matrix[3]+=(D1*B2x*B2x + D3*B2y*B2y)*WiWjdetJ;
	Matrix[4]+=(B2x*B1y*D2 + B1x*B2y*D3)*WiWjdetJ;
	Matrix[5]+=(B1x*B2x*D1 + B1y*B2y*D3)*WiWjdetJ;
	Matrix[6]+=(D3*B2x*B2x + D1*B2y*B2y)*WiWjdetJ;
	Matrix[7]+=(B2x*B2y*D2 + B2x*B2y*D3)*WiWjdetJ;
	Matrix[8]+=(B1x*B2x*D3 + B1y*B2y*D1)*WiWjdetJ;
	Matrix[9]+=(B1x*B2y*D2 + B2x*B1y*D3)*WiWjdetJ;
	Matrix[10]+=(D1*B3x*B3x + D3*B3y*B3y)*WiWjdetJ;
	Matrix[11]+=(B3x*B2y*D2 + B2x*B3y*D3)*WiWjdetJ;
	Matrix[12]+=(B2x*B3x*D1 + B2y*B3y*D3)*WiWjdetJ;
	Matrix[13]+=(B3x*B1y*D2 + B1x*B3y*D3)*WiWjdetJ;
	Matrix[14]+=(B1x*B3x*D1 + B1y*B3y*D3)*WiWjdetJ;
	Matrix[15]+=(D3*B3x*B3x + D1*B3y*B3y)*WiWjdetJ;
	Matrix[16]+=(B3x*B3y*D2 + B3x*B3y*D3)*WiWjdetJ;
	Matrix[17]+=(B2x*B3x*D3 + B2y*B3y*D1)*WiWjdetJ;
	Matrix[18]+=(B2x*B3y*D2 + B3x*B2y*D3)*WiWjdetJ;
	Matrix[19]+=(B1x*B3x*D3 + B1y*B3y*D1)*WiWjdetJ;
	Matrix[20]+=(B1x*B3y*D2 + B3x*B1y*D3)*WiWjdetJ;
	Matrix[21]+=(D1*B4x*B4x + D3*B4y*B4y)*WiWjdetJ;
	Matrix[22]+=(B4x*B3y*D2 + B3x*B4y*D3)*WiWjdetJ;
	Matrix[23]+=(B3x*B4x*D1 + B3y*B4y*D3)*WiWjdetJ;
	Matrix[24]+=(B4x*B2y*D2 + B2x*B4y*D3)*WiWjdetJ;
	Matrix[25]+=(B2x*B4x*D1 + B2y*B4y*D3)*WiWjdetJ;
	Matrix[26]+=(B4x*B1y*D2 + B1x*B4y*D3)*WiWjdetJ;
	Matrix[27]+=(B1x*B4x*D1 + B1y*B4y*D3)*WiWjdetJ;
	Matrix[28]+=(D3*B4x*B4x + D1*B4y*B4y)*WiWjdetJ;
	Matrix[29]+=(B4x*B4y*D2 + B4x*B4y*D3)*WiWjdetJ;
	Matrix[30]+=(B3x*B4x*D3 + B3y*B4y*D1)*WiWjdetJ;
	Matrix[31]+=(B3x*B4y*D2 + B4x*B3y*D3)*WiWjdetJ;
	Matrix[32]+=(B2x*B4x*D3 + B2y*B4y*D1)*WiWjdetJ;
	Matrix[33]+=(B2x*B4y*D2 + B4x*B2y*D3)*WiWjdetJ;
	Matrix[34]+=(B1x*B4x*D3 + B1y*B4y*D1)*WiWjdetJ;
	Matrix[35]+=(B1x*B4y*D2 + B4x*B1y*D3)*WiWjdetJ;
	}
	}
}

void CQ4::ElementStress(double* stress, double* Displacement){
	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);	// Pointer to material of the element
	double E_=material_->E;
	double nu_=material_->nu;
	double D1=E_/(1-nu_*nu_);
	double D2=E_*nu_/(1-nu_*nu_);
	double D3=E_/2/(1+nu_);

	double Disp_element[8];
	for(int i=0;i<8;i++){
		if(!LocationMatrix_[i])
			Disp_element[i]=0;
		else
			Disp_element[i]=Displacement[LocationMatrix_[i]-1];
	}

	double gausspoint[2]={-0.57735027, 0.57735027};	//Gauss point when ngp=2
	unsigned int num=0;
	for(int i=0;i<2;i++){
		double eta=gausspoint[i];
		for(int j=0;j<2;j++){
			double psi=gausspoint[j];
			double X[4]={nodes_[0]->XYZ[0],nodes_[1]->XYZ[0],nodes_[2]->XYZ[0],nodes_[3]->XYZ[0]};
			double Y[4]={nodes_[0]->XYZ[1],nodes_[1]->XYZ[1],nodes_[2]->XYZ[1],nodes_[3]->XYZ[1]};
			double GN[2][4]={(eta-1)/4,(1-eta)/4,(1+eta)/4,(-eta-1)/4,(psi-1)/4,(-psi-1)/4,(1+psi)/4,(1-psi)/4};
			double Jocbian[2][2]={GN[0][0]*X[0]+GN[0][1]*X[1]+GN[0][2]*X[2]+GN[0][3]*X[3],GN[0][0]*Y[0]+GN[0][1]*Y[1]+GN[0][2]*Y[2]+GN[0][3]*Y[3],
			GN[1][0]*X[0]+GN[1][1]*X[1]+GN[1][2]*X[2]+GN[1][3]*X[3],GN[1][0]*Y[0]+GN[1][1]*Y[1]+GN[1][2]*Y[2]+GN[1][3]*Y[3]};
			double detJocbian=Jocbian[0][0]*Jocbian[1][1]-Jocbian[1][0]*Jocbian[0][1];
			double Jocbian_inv[2][2]={Jocbian[1][1]/detJocbian,-Jocbian[0][1]/detJocbian,-Jocbian[1][0]/detJocbian,Jocbian[0][0]/detJocbian};
			double B1x=Jocbian_inv[0][0]*GN[0][0]+Jocbian_inv[0][1]*GN[1][0];
			double B2x=Jocbian_inv[0][0]*GN[0][1]+Jocbian_inv[0][1]*GN[1][1];
			double B3x=Jocbian_inv[0][0]*GN[0][2]+Jocbian_inv[0][1]*GN[1][2];
			double B4x=Jocbian_inv[0][0]*GN[0][3]+Jocbian_inv[0][1]*GN[1][3];
			double B1y=Jocbian_inv[1][0]*GN[0][0]+Jocbian_inv[1][1]*GN[1][0];
			double B2y=Jocbian_inv[1][0]*GN[0][1]+Jocbian_inv[1][1]*GN[1][1];
			double B3y=Jocbian_inv[1][0]*GN[0][2]+Jocbian_inv[1][1]*GN[1][2];
			double B4y=Jocbian_inv[1][0]*GN[0][3]+Jocbian_inv[1][1]*GN[1][3];
			stress[num++]=B1x*D1*Disp_element[0]+B1y*D2*Disp_element[1]+B2x*D1*Disp_element[2]+B2y*D2*Disp_element[3]+B3x*D1*Disp_element[4]+B3y*D2*Disp_element[5]+B4x*D1*Disp_element[6]+B4y*D2*Disp_element[7];
			stress[num++]=B1x*D2*Disp_element[0]+B1y*D1*Disp_element[1]+B2x*D2*Disp_element[2]+B2y*D1*Disp_element[3]+B3x*D2*Disp_element[4]+B3y*D1*Disp_element[5]+B4x*D2*Disp_element[6]+B4y*D1*Disp_element[7];
			stress[num++]=B1y*D3*Disp_element[0]+B1x*D3*Disp_element[1]+B2y*D3*Disp_element[2]+B2x*D3*Disp_element[3]+B3y*D3*Disp_element[4]+B3x*D3*Disp_element[5]+B4y*D3*Disp_element[6]+B4x*D3*Disp_element[7];
		}
	}

}

void CQ4::ElementCoord (double* coord){
	double gausspoint[2]={-0.57735027, 0.57735027};	//Gauss point when ngp=2
	unsigned int num=0;
	for(int i=0;i<2;i++){
		double eta=gausspoint[i];
		for(int j=0;j<2;j++){
			double psi=gausspoint[j];
			coord[num++]=0.25*((1-psi)*(1-eta)*nodes_[0]->XYZ[0]+(1+psi)*(1-eta)*nodes_[1]->XYZ[0]+(1+psi)*(1+eta)*nodes_[2]->XYZ[0]+(1-psi)*(1+eta)*nodes_[3]->XYZ[0]);
			coord[num++]=0.25*((1-psi)*(1-eta)*nodes_[0]->XYZ[1]+(1+psi)*(1-eta)*nodes_[1]->XYZ[1]+(1+psi)*(1+eta)*nodes_[2]->XYZ[1]+(1-psi)*(1+eta)*nodes_[3]->XYZ[1]);
		}
	}
}
