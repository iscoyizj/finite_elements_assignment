#include "Shell.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

CShell::CShell(){
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 24;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CShell::~CShell(){}

bool CShell::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList){
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
	unsigned int N1, N2, N3, N4;	// 4 nodal numbers in counterclock order

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = dynamic_cast<CShellMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

void CShell::Write(COutputter& output, unsigned int Ele){
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber << setw(9) 
		<< nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

void CShell::GenerateLocationMatrix(){
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < 6; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

unsigned int CShell::SizeOfStiffnessMatrix() { return 300; }

void CShell::ElementStiffness(double* Matrix){
	clear(Matrix, SizeOfStiffnessMatrix());
	double X[4]={nodes_[0]->XYZ[0],nodes_[1]->XYZ[0],nodes_[2]->XYZ[0],nodes_[3]->XYZ[0]};
	double Y[4]={nodes_[0]->XYZ[1],nodes_[1]->XYZ[1],nodes_[2]->XYZ[1],nodes_[3]->XYZ[1]};
	double Z[4]={nodes_[0]->XYZ[2],nodes_[1]->XYZ[2],nodes_[2]->XYZ[2],nodes_[3]->XYZ[2]};

	//Calculate the unit normal vector for each node
	double l3[4], m3[4], n3[4];
	l3[0]=(Z[0]*(Y[0]-Y[1])-Z[0]*(Y[0]-Y[3])+Z[1]*(Y[0]-Y[3])-Z[3]*(Y[0]-Y[1]))/4;
	m3[0]=(X[0]*(Z[0]-Z[1])-X[0]*(Z[0]-Z[3])+X[1]*(Z[0]-Z[3])-X[3]*(Z[0]-Z[1]))/4;
	n3[0]=(Y[0]*(X[0]-X[1])-Y[0]*(X[0]-X[3])+Y[1]*(X[0]-X[3])-Y[3]*(X[0]-X[1]))/4;
	l3[1]=(Z[1]*(Y[0]-Y[1])-Z[0]*(Y[1]-Y[2])-Z[2]*(Y[0]-Y[1])+Z[1]*(Y[1]-Y[2]))/4;
	m3[1]=(X[1]*(Z[0]-Z[1])-X[0]*(Z[1]-Z[2])-X[2]*(Z[0]-Z[1])+X[1]*(Z[1]-Z[2]))/4;
	n3[1]=(Y[1]*(X[0]-X[1])-Y[0]*(X[1]-X[2])-Y[2]*(X[0]-X[1])+Y[1]*(X[1]-X[2]))/4;
	l3[2]=(Z[2]*(Y[1]-Y[2])-Z[1]*(Y[2]-Y[3])-Z[3]*(Y[1]-Y[2])+Z[2]*(Y[2]-Y[3]))/4;
	m3[2]=(X[2]*(Z[1]-Z[2])-X[1]*(Z[2]-Z[3])-X[3]*(Z[1]-Z[2])+X[2]*(Z[2]-Z[3]))/4;
	n3[2]=(Y[2]*(X[1]-X[2])-Y[1]*(X[2]-X[3])-Y[3]*(X[1]-X[2])+Y[2]*(X[2]-X[3]))/4;
	l3[3]=(Z[2]*(Y[0]-Y[3])-Z[0]*(Y[2]-Y[3])-Z[3]*(Y[0]-Y[3])+Z[3]*(Y[2]-Y[3]))/4;
	m3[3]=(X[2]*(Z[0]-Z[3])-X[0]*(Z[2]-Z[3])-X[3]*(Z[0]-Z[3])+X[3]*(Z[2]-Z[3]))/4;
	n3[3]=(Y[2]*(X[0]-X[3])-Y[0]*(X[2]-X[3])-Y[3]*(X[0]-X[3])+Y[3]*(X[2]-X[3]))/4;
	for(unsigned int loop=0;loop<4;loop++){
		double lmn3=sqrt(l3[loop]*l3[loop]+m3[loop]*m3[loop]+n3[loop]*n3[loop]);
		l3[loop]/=lmn3;
		m3[loop]/=lmn3;
		n3[loop]/=lmn3;
	}
	double l1[4], m1[4], n1[4];
	double l2[4], m2[4], n2[4];
	for(unsigned int loop=0;loop<4;loop++){
	double lmn1=sqrt(n3[loop]*n3[loop]+l3[loop]*l3[loop]);
	l1[loop]=n3[loop]/lmn1;
	m1[loop]=0;
	n1[loop]=-l3[loop]/lmn1;
	l2[loop]=-m3[loop]*l3[loop]/lmn1;
	m2[loop]=(n3[loop]*n3[loop]+l3[loop]*l3[loop])/lmn1;
	n2[loop]=-m3[loop]*n3[loop]/lmn1;
	}
	//Calculate magnititude of non-zero elements of constitutive matrix
	//D=[D1,D2,0,0,0;D2,D1,0,0,0;0,0,D3,0,0;0,0,0,D3,0;0,0,0,0,D3], dimension=5*5
	CShellMaterial* material_ = dynamic_cast<CShellMaterial*>(ElementMaterial_);
	double E_=material_->E;
	double nu_=material_->nu;
	double D1=E_/(1-nu_*nu_);
	double D2=D1*nu_;
	double D3=E_/2/(1+nu_);

	double thick=material_->thick;

	double gausspoint[2]={-0.57735027, 0.57735027};	//Gauss point when ngp=2
	double weight_gauss[2]={1,1};	//Gauss weight

	for(unsigned int i=0;i<2;i++){
		double xi=gausspoint[i];
		for(unsigned int j=0;j<2;j++){
			double eta=gausspoint[j];
			for(unsigned int k=0;k<2;k++){
				double zeta=gausspoint[k];
	double N[4]={(1-xi)*(1-eta)/4,(1+xi)*(1-eta)/4,(1+xi)*(1+eta)/4,(1-xi)*(1+eta)/4};
	double GNxi[4]={-(1-eta)/4,(1-eta)/4,(1+eta)/4,-(1+eta)/4};		//GNxi=dN/dxi
	double GNeta[4]={-(1-xi)/4,-(1+xi)/4,(1+xi)/4,(1-xi)/4};		//GNeta=dN/deta
	double J[3][3]={0,0,0,0,0,0,0,0,0};			//Calculate the Jocabian Matrix
	for(unsigned int loop=0;loop<4;loop++){
	J[0][0]+=GNxi[loop]*(X[loop]+0.5*zeta*thick*l3[loop]);
	J[0][1]+=GNxi[loop]*(Y[loop]+0.5*zeta*thick*m3[loop]);
	J[0][2]+=GNxi[loop]*(Z[loop]+0.5*zeta*thick*n3[loop]);
	J[1][0]+=GNeta[loop]*(X[loop]+0.5*zeta*thick*l3[loop]);
	J[1][1]+=GNeta[loop]*(Y[loop]+0.5*zeta*thick*m3[loop]);
	J[1][2]+=GNeta[loop]*(Z[loop]+0.5*zeta*thick*n3[loop]);
	J[2][0]+=N[loop]*0.5*thick*l3[loop];
	J[2][1]+=N[loop]*0.5*thick*m3[loop];
	J[2][2]+=N[loop]*0.5*thick*n3[loop];
	}
	double detJ=J[0][0]*J[1][1]*J[2][2]-J[0][0]*J[1][2]*J[2][1]-J[0][1]*J[1][0]*J[2][2]+J[0][1]*J[1][2]*J[2][0]+J[0][2]*J[1][0]*J[2][1]-J[0][2]*J[1][1]*J[2][0];
	double invJ[3][3];		//Calculate the inverse of the Jocabian Matrix
	invJ[0][0]=(J[1][1]*J[2][2]-J[1][2]*J[2][1])/detJ;
	invJ[0][1]=(J[0][2]*J[2][1]-J[0][1]*J[2][2])/detJ;
	invJ[0][2]=(J[0][1]*J[1][2]-J[0][2]*J[1][1])/detJ;
	invJ[1][0]=(J[1][2]*J[2][0]-J[1][0]*J[2][2])/detJ;
	invJ[1][1]=(J[0][0]*J[2][2]-J[0][2]*J[2][0])/detJ;
	invJ[1][2]=(J[0][2]*J[1][0]-J[0][0]*J[1][2])/detJ;
	invJ[2][0]=(J[1][0]*J[2][1]-J[1][1]*J[2][0])/detJ;
	invJ[2][1]=(J[0][1]*J[2][0]-J[0][0]*J[2][1])/detJ;
	invJ[2][2]=(J[0][0]*J[1][1]-J[0][1]*J[1][0])/detJ;
	double A[4][5][5];
	for(unsigned int loop=0;loop<4;loop++){
		A[loop][0][0]=GNeta[loop]*invJ[0][1]+GNxi[loop]*invJ[0][0];
		A[loop][0][1]=0.0;
		A[loop][0][2]=0.0;
		A[loop][0][3]=-N[loop]*invJ[0][2]*l2[loop]*thick/2-GNeta[loop]*invJ[0][1]*l2[loop]*thick/2*zeta-GNxi[loop]*invJ[0][0]*l2[loop]*thick/2*zeta;
		A[loop][0][4]=N[loop]*invJ[0][2]*l1[loop]*thick/2+GNeta[loop]*invJ[0][1]*l1[loop]*thick/2*zeta+GNxi[loop]*invJ[0][0]*l1[loop]*thick/2*zeta;
		A[loop][1][0]=0.0;
		A[loop][1][1]=GNeta[loop]*invJ[1][1]+GNxi[loop]*invJ[1][0];
		A[loop][1][2]=0.0;
		A[loop][1][3]=-N[loop]*invJ[1][2]*m2[loop]*thick/2-GNeta[loop]*invJ[1][1]*m2[loop]*thick/2*zeta-GNxi[loop]*invJ[1][0]*m2[loop]*thick/2*zeta;
		A[loop][1][4]=N[loop]*invJ[1][2]*m1[loop]*thick/2+GNeta[loop]*invJ[1][1]*m1[loop]*thick/2*zeta+GNxi[loop]*invJ[1][0]*m1[loop]*thick/2*zeta;
		A[loop][2][0]=GNeta[loop]*invJ[1][1]+GNxi[loop]*invJ[1][0];
		A[loop][2][1]=GNeta[loop]*invJ[0][1]+GNxi[loop]*invJ[0][0];
		A[loop][2][2]=0.0;
		A[loop][2][3]=-N[loop]*invJ[1][2]*l2[loop]*thick/2-N[loop]*invJ[0][2]*m2[loop]*thick/2-GNeta[loop]*invJ[1][1]*l2[loop]*thick/2*zeta-GNxi[loop]*invJ[1][0]*l2[loop]*thick/2*zeta-GNeta[loop]*invJ[0][1]*m2[loop]*thick/2*zeta-GNxi[loop]*invJ[0][0]*m2[loop]*thick/2*zeta;
		A[loop][2][4]=N[loop]*invJ[1][2]*l1[loop]*thick/2+N[loop]*invJ[0][2]*m1[loop]*thick/2+GNeta[loop]*invJ[1][1]*l1[loop]*thick/2*zeta+GNxi[loop]*invJ[1][0]*l1[loop]*thick/2*zeta+GNeta[loop]*invJ[0][1]*m1[loop]*thick/2*zeta+GNxi[loop]*invJ[0][0]*m1[loop]*thick/2*zeta;
		A[loop][3][0]=0.0;
		A[loop][3][1]=GNeta[loop]*invJ[2][1]+GNxi[loop]*invJ[2][0];
		A[loop][3][2]=GNeta[loop]*invJ[1][1]+GNxi[loop]*invJ[1][0];
		A[loop][3][3]=-N[loop]*invJ[2][2]*m2[loop]*thick/2-N[loop]*invJ[1][2]*n2[loop]*thick/2-GNeta[loop]*invJ[2][1]*m2[loop]*thick/2*zeta-GNxi[loop]*invJ[2][0]*m2[loop]*thick/2*zeta-GNeta[loop]*invJ[1][1]*n2[loop]*thick/2*zeta-GNxi[loop]*invJ[1][0]*n2[loop]*thick/2*zeta;
		A[loop][3][4]=N[loop]*invJ[2][2]*m1[loop]*thick/2+N[loop]*invJ[1][2]*n1[loop]*thick/2+GNeta[loop]*invJ[2][1]*m1[loop]*thick/2*zeta+GNxi[loop]*invJ[2][0]*m1[loop]*thick/2*zeta+GNeta[loop]*invJ[1][1]*n1[loop]*thick/2*zeta+GNxi[loop]*invJ[1][0]*n1[loop]*thick/2*zeta;
		A[loop][4][0]=GNeta[loop]*invJ[2][1]+GNxi[loop]*invJ[2][0];
		A[loop][4][1]=0.0;
		A[loop][4][2]=GNeta[loop]*invJ[0][1]+GNxi[loop]*invJ[0][0];
		A[loop][4][3]=-N[loop]*invJ[2][2]*l2[loop]*thick/2-N[loop]*invJ[0][2]*n2[loop]*thick/2-GNeta[loop]*invJ[2][1]*l2[loop]*thick/2*zeta-GNxi[loop]*invJ[2][0]*l2[loop]*thick/2*zeta-GNeta[loop]*invJ[0][1]*n2[loop]*thick/2*zeta-GNxi[loop]*invJ[0][0]*n2[loop]*thick/2*zeta;
		A[loop][4][4]=N[loop]*invJ[2][2]*l1[loop]*thick/2+N[loop]*invJ[0][2]*n1[loop]*thick/2+GNeta[loop]*invJ[2][1]*l1[loop]*thick/2*zeta+GNxi[loop]*invJ[2][0]*l1[loop]*thick/2*zeta+GNeta[loop]*invJ[0][1]*n1[loop]*thick/2*zeta+GNxi[loop]*invJ[0][0]*n1[loop]*thick/2*zeta;
	}

	double Be[4][5][6];
	for(unsigned int loop=0;loop<4;loop++){
		double A1=A[loop][0][0];
		double A2=A[loop][0][3];
		double A3=A[loop][0][4];
		double A4=A[loop][1][1];
		double A5=A[loop][1][3];
		double A6=A[loop][1][4];
		double A7=A[loop][2][0];
		double A8=A[loop][2][1];
		double A9=A[loop][2][3];
		double A10=A[loop][2][4];
		double A11=A[loop][3][1];
		double A12=A[loop][3][2];
		double A13=A[loop][3][3];
		double A14=A[loop][3][4];
		double A15=A[loop][4][0];
		double A16=A[loop][4][2];
		double A17=A[loop][4][3];
		double A18=A[loop][4][4];
		//Trans=[1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,ts1,ts2,ts3;0,0,0,ts4,ts5,ts6]
		double ts1=m2[loop]*n3[loop]-m3[loop]*n2[loop];
		double ts2=l3[loop]*n2[loop]-l2[loop]*n3[loop];
		double ts3=l2[loop]*m3[loop]-l3[loop]*m2[loop];
		double ts4=m3[loop]*n1[loop]-m1[loop]*n3[loop];
		double ts5=l1[loop]*n3[loop]-l3[loop]*n1[loop];
		double ts6=l3[loop]*m1[loop]-l1[loop]*m3[loop];

		Be[loop][0][0]=A1;
		Be[loop][0][1]=0.0;
		Be[loop][0][2]=0.0;
		Be[loop][0][3]=A2*ts1 + A3*ts4;
		Be[loop][0][4]=A2*ts2 + A3*ts5;
		Be[loop][0][5]=A2*ts3 + A3*ts6;
		Be[loop][1][0]=0.0;
		Be[loop][1][1]=A4;
		Be[loop][1][2]=0.0;
		Be[loop][1][3]=A5*ts1 + A6*ts4;
		Be[loop][1][4]=A5*ts2 + A6*ts5;
		Be[loop][1][5]=A5*ts3 + A6*ts6;
		Be[loop][2][0]=A7;
		Be[loop][2][1]=A8;
		Be[loop][2][2]=0.0;
		Be[loop][2][3]=A9*ts1 + A10*ts4;
		Be[loop][2][4]=A9*ts2 + A10*ts5;
		Be[loop][2][5]=A9*ts3 + A10*ts6;
		Be[loop][3][0]=0.0;
		Be[loop][3][1]=A11;
		Be[loop][3][2]=A12;
		Be[loop][3][3]=A13*ts1 + A14*ts4;
		Be[loop][3][4]=A13*ts2 + A14*ts5;
		Be[loop][3][5]=A13*ts3 + A14*ts6;
		Be[loop][4][0]=A15;
		Be[loop][4][1]=0.0;
		Be[loop][4][2]=A16;
		Be[loop][4][3]=A17*ts1 + A18*ts4;
		Be[loop][4][4]=A17*ts2 + A18*ts5;
		Be[loop][4][5]=A17*ts3 + A18*ts6;
	}

	double K[10][6][6];
	unsigned int num=0;
	for(unsigned int loop1=0;loop1<4;loop1++){
		for(unsigned int loop2=loop1;loop2<4;loop2++){
			//Be(loop1)T*D*Be(loop2), dimension of Be=5*6
			double P1=Be[loop1][0][0];
			double P2=Be[loop1][0][3];
			double P3=Be[loop1][0][4];
			double P4=Be[loop1][0][5];
			double P5=Be[loop1][1][1];
			double P6=Be[loop1][1][3];
			double P7=Be[loop1][1][4];
			double P8=Be[loop1][1][5];
			double P9=Be[loop1][2][0];
			double P10=Be[loop1][2][1];
			double P11=Be[loop1][2][3];
			double P12=Be[loop1][2][4];
			double P13=Be[loop1][2][5];
			double P14=Be[loop1][3][1];
			double P15=Be[loop1][3][2];
			double P16=Be[loop1][3][3];
			double P17=Be[loop1][3][4];
			double P18=Be[loop1][3][5];
			double P19=Be[loop1][4][0];
			double P20=Be[loop1][4][2];
			double P21=Be[loop1][4][3];
			double P22=Be[loop1][4][4];
			double P23=Be[loop1][4][5];

			double Q1=Be[loop2][0][0];
			double Q2=Be[loop2][0][3];
			double Q3=Be[loop2][0][4];
			double Q4=Be[loop2][0][5];
			double Q5=Be[loop2][1][1];
			double Q6=Be[loop2][1][3];
			double Q7=Be[loop2][1][4];
			double Q8=Be[loop2][1][5];
			double Q9=Be[loop2][2][0];
			double Q10=Be[loop2][2][1];
			double Q11=Be[loop2][2][3];
			double Q12=Be[loop2][2][4];
			double Q13=Be[loop2][2][5];
			double Q14=Be[loop2][3][1];
			double Q15=Be[loop2][3][2];
			double Q16=Be[loop2][3][3];
			double Q17=Be[loop2][3][4];
			double Q18=Be[loop2][3][5];
			double Q19=Be[loop2][4][0];
			double Q20=Be[loop2][4][2];
			double Q21=Be[loop2][4][3];
			double Q22=Be[loop2][4][4];
			double Q23=Be[loop2][4][5];

			K[num][0][0]=D1*P1*Q1 + D3*P9*Q9 + D3*P19*Q19;
			K[num][0][1]=D2*P1*Q5 + D3*P9*Q10;
			K[num][0][2]=D3*P19*Q20;
			K[num][0][3]=D1*P1*Q2 + D2*P1*Q6 + D3*P9*Q11 + D3*P19*Q21;
			K[num][0][4]=D1*P1*Q3 + D2*P1*Q7 + D3*P9*Q12 + D3*P19*Q22;
			K[num][0][5]=D1*P1*Q4 + D2*P1*Q8 + D3*P9*Q13 + D3*P19*Q23;
			K[num][1][0]=D2*P5*Q1 + D3*P10*Q9;
			K[num][1][1]=D1*P5*Q5 + D3*P10*Q10 + D3*P14*Q14;
			K[num][1][2]=D3*P14*Q15;
			K[num][1][3]=D2*P5*Q2 + D1*P5*Q6 + D3*P10*Q11 + D3*P14*Q16;
			K[num][1][4]=D2*P5*Q3 + D1*P5*Q7 + D3*P10*Q12 + D3*P14*Q17;
			K[num][1][5]=D2*P5*Q4 + D1*P5*Q8 + D3*P10*Q13 + D3*P14*Q18;
			K[num][2][0]=D3*P20*Q19;
			K[num][2][1]=D3*P15*Q14;
			K[num][2][2]=D3*P15*Q15 + D3*P20*Q20;
			K[num][2][3]=D3*P15*Q16 + D3*P20*Q21;
			K[num][2][4]=D3*P15*Q17 + D3*P20*Q22;
			K[num][2][5]=D3*P15*Q18 + D3*P20*Q23;
			K[num][3][0]=Q1*(D1*P2 + D2*P6) + D3*P11*Q9 + D3*P21*Q19;
			K[num][3][1]=Q5*(D2*P2 + D1*P6) + D3*P11*Q10 + D3*P16*Q14;
			K[num][3][2]=D3*P16*Q15 + D3*P21*Q20;
			K[num][3][3]=Q2*(D1*P2 + D2*P6) + Q6*(D2*P2 + D1*P6) + D3*P11*Q11 + D3*P16*Q16 + D3*P21*Q21;
			K[num][3][4]=Q3*(D1*P2 + D2*P6) + Q7*(D2*P2 + D1*P6) + D3*P11*Q12 + D3*P16*Q17 + D3*P21*Q22;
			K[num][3][5]=Q4*(D1*P2 + D2*P6) + Q8*(D2*P2 + D1*P6) + D3*P11*Q13 + D3*P16*Q18 + D3*P21*Q23;
			K[num][4][0]=Q1*(D1*P3 + D2*P7) + D3*P12*Q9 + D3*P22*Q19;
			K[num][4][1]=Q5*(D2*P3 + D1*P7) + D3*P12*Q10 + D3*P17*Q14, D3*P17*Q15 + D3*P22*Q20;
			K[num][4][2]=D3*P17*Q15 + D3*P22*Q20;
			K[num][4][3]=Q2*(D1*P3 + D2*P7) + Q6*(D2*P3 + D1*P7) + D3*P12*Q11 + D3*P17*Q16 + D3*P22*Q21;
			K[num][4][4]=Q3*(D1*P3 + D2*P7) + Q7*(D2*P3 + D1*P7) + D3*P12*Q12 + D3*P17*Q17 + D3*P22*Q22;
			K[num][4][5]=Q4*(D1*P3 + D2*P7) + Q8*(D2*P3 + D1*P7) + D3*P12*Q13 + D3*P17*Q18 + D3*P22*Q23;
			K[num][5][0]=Q1*(D1*P4 + D2*P8) + D3*P13*Q9 + D3*P23*Q19;
			K[num][5][1]=Q5*(D2*P4 + D1*P8) + D3*P13*Q10 + D3*P18*Q14;
			K[num][5][2]=D3*P18*Q15 + D3*P23*Q20;
			K[num][5][3]=Q2*(D1*P4 + D2*P8) + Q6*(D2*P4 + D1*P8) + D3*P13*Q11 + D3*P18*Q16 + D3*P23*Q21;
			K[num][5][4]=Q3*(D1*P4 + D2*P8) + Q7*(D2*P4 + D1*P8) + D3*P13*Q12 + D3*P18*Q17 + D3*P23*Q22;
			K[num][5][5]=Q4*(D1*P4 + D2*P8) + Q8*(D2*P4 + D1*P8) + D3*P13*Q13 + D3*P18*Q18 + D3*P23*Q23;

			num++;
		}
	}
	double ke[300];
	ke[0]=K[0][0][0];
	ke[1]=K[0][1][1];
	ke[2]=K[0][0][1];
	ke[3]=K[0][2][2];
	ke[4]=K[0][1][2];
	ke[5]=K[0][0][2];
	ke[6]=K[0][3][3];
	ke[7]=K[0][2][3];
	ke[8]=K[0][1][3];
	ke[9]=K[0][0][3];
	ke[10]=K[0][4][4];
	ke[11]=K[0][3][4];
	ke[12]=K[0][2][4];
	ke[13]=K[0][1][4];
	ke[14]=K[0][0][4];
	ke[15]=K[0][5][5];
	ke[16]=K[0][4][5];
	ke[17]=K[0][3][5];
	ke[18]=K[0][2][5];
	ke[19]=K[0][1][5];
	ke[20]=K[0][0][5];
	ke[21]=K[4][0][0];
	ke[22]=K[1][5][0];
	ke[23]=K[1][4][0];
	ke[24]=K[1][3][0];
	ke[25]=K[1][2][0];
	ke[26]=K[1][1][0];
	ke[27]=K[1][0][0];
	ke[28]=K[4][1][1];
	ke[29]=K[4][0][1];
	ke[30]=K[1][5][1];
	ke[31]=K[1][4][1];
	ke[32]=K[1][3][1];
	ke[33]=K[1][2][1];
	ke[34]=K[1][1][1];
	ke[35]=K[1][0][1];
	ke[36]=K[4][2][2];
	ke[37]=K[4][1][2];
	ke[38]=K[4][0][2];
	ke[39]=K[1][5][2];
	ke[40]=K[1][4][2];
	ke[41]=K[1][3][2];
	ke[42]=K[1][2][2];
	ke[43]=K[1][1][2];
	ke[44]=K[1][0][2];
	ke[45]=K[4][3][3];
	ke[46]=K[4][2][3];
	ke[47]=K[4][1][3];
	ke[48]=K[4][0][3];
	ke[49]=K[1][5][3];
	ke[50]=K[1][4][3];
	ke[51]=K[1][3][3];
	ke[52]=K[1][2][3];
	ke[53]=K[1][1][3];
	ke[54]=K[1][0][3];
	ke[55]=K[4][4][4];
	ke[56]=K[4][3][4];
	ke[57]=K[4][2][4];
	ke[58]=K[4][1][4];
	ke[59]=K[4][0][4];
	ke[60]=K[1][5][4];
	ke[61]=K[1][4][4];
	ke[62]=K[1][3][4];
	ke[63]=K[1][2][4];
	ke[64]=K[1][1][4];
	ke[65]=K[1][0][4];
	ke[66]=K[4][5][5];
	ke[67]=K[4][4][5];
	ke[68]=K[4][3][5];
	ke[69]=K[4][2][5];
	ke[70]=K[4][1][5];
	ke[71]=K[4][0][5];
	ke[72]=K[1][5][5];
	ke[73]=K[1][4][5];
	ke[74]=K[1][3][5];
	ke[75]=K[1][2][5];
	ke[76]=K[1][1][5];
	ke[77]=K[1][0][5];
	ke[78]=K[7][0][0];
	ke[79]=K[5][5][0];
	ke[80]=K[5][4][0];
	ke[81]=K[5][3][0];
	ke[82]=K[5][2][0];
	ke[83]=K[5][1][0];
	ke[84]=K[5][0][0];
	ke[85]=K[2][5][0];
	ke[86]=K[2][4][0];
	ke[87]=K[2][3][0];
	ke[88]=K[2][2][0];
	ke[89]=K[2][1][0];
	ke[90]=K[2][0][0];
	ke[91]=K[7][1][1];
	ke[92]=K[7][0][1];
	ke[93]=K[5][5][1];
	ke[94]=K[5][4][1];
	ke[95]=K[5][3][1];
	ke[96]=K[5][2][1];
	ke[97]=K[5][1][1];
	ke[98]=K[5][0][1];
	ke[99]=K[2][5][1];
	ke[100]=K[2][4][1];
	ke[101]=K[2][3][1];
	ke[102]=K[2][2][1];
	ke[103]=K[2][1][1];
	ke[104]=K[2][0][1];
	ke[105]=K[7][2][2];
	ke[106]=K[7][1][2];
	ke[107]=K[7][0][2];
	ke[108]=K[5][5][2];
	ke[109]=K[5][4][2];
	ke[110]=K[5][3][2];
	ke[111]=K[5][2][2];
	ke[112]=K[5][1][2];
	ke[113]=K[5][0][2];
	ke[114]=K[2][5][2];
	ke[115]=K[2][4][2];
	ke[116]=K[2][3][2];
	ke[117]=K[2][2][2];
	ke[118]=K[2][1][2];
	ke[119]=K[2][0][2];
	ke[120]=K[7][3][3];
	ke[121]=K[7][2][3];
	ke[122]=K[7][1][3];
	ke[123]=K[7][0][3];
	ke[124]=K[5][5][3];
	ke[125]=K[5][4][3];
	ke[126]=K[5][3][3];
	ke[127]=K[5][2][3];
	ke[128]=K[5][1][3];
	ke[129]=K[5][0][3];
	ke[130]=K[2][5][3];
	ke[131]=K[2][4][3];
	ke[132]=K[2][3][3];
	ke[133]=K[2][2][3];
	ke[134]=K[2][1][3];
	ke[135]=K[2][0][3];
	ke[136]=K[7][4][4];
	ke[137]=K[7][3][4];
	ke[138]=K[7][2][4];
	ke[139]=K[7][1][4];
	ke[140]=K[7][0][4];
	ke[141]=K[5][5][4];
	ke[142]=K[5][4][4];
	ke[143]=K[5][3][4];
	ke[144]=K[5][2][4];
	ke[145]=K[5][1][4];
	ke[146]=K[5][0][4];
	ke[147]=K[2][5][4];
	ke[148]=K[2][4][4];
	ke[149]=K[2][3][4];
	ke[150]=K[2][2][4];
	ke[151]=K[2][1][4];
	ke[152]=K[2][0][4];
	ke[153]=K[7][5][5];
	ke[154]=K[7][4][5];
	ke[155]=K[7][3][5];
	ke[156]=K[7][2][5];
	ke[157]=K[7][1][5];
	ke[158]=K[7][0][5];
	ke[159]=K[5][5][5];
	ke[160]=K[5][4][5];
	ke[161]=K[5][3][5];
	ke[162]=K[5][2][5];
	ke[163]=K[5][1][5];
	ke[164]=K[5][0][5];
	ke[165]=K[2][5][5];
	ke[166]=K[2][4][5];
	ke[167]=K[2][3][5];
	ke[168]=K[2][2][5];
	ke[169]=K[2][1][5];
	ke[170]=K[2][0][5];
	ke[171]=K[9][0][0];
	ke[172]=K[8][5][0];
	ke[173]=K[8][4][0];
	ke[174]=K[8][3][0];
	ke[175]=K[8][2][0];
	ke[176]=K[8][1][0];
	ke[177]=K[8][0][0];
	ke[178]=K[6][5][0];
	ke[179]=K[6][4][0];
	ke[180]=K[6][3][0];
	ke[181]=K[6][2][0];
	ke[182]=K[6][1][0];
	ke[183]=K[6][0][0];
	ke[184]=K[3][5][0];
	ke[185]=K[3][4][0];
	ke[186]=K[3][3][0];
	ke[187]=K[3][2][0];
	ke[188]=K[3][1][0];
	ke[189]=K[3][0][0];
	ke[190]=K[9][1][1];
	ke[191]=K[9][0][1];
	ke[192]=K[8][5][1];
	ke[193]=K[8][4][1];
	ke[194]=K[8][3][1];
	ke[195]=K[8][2][1];
	ke[196]=K[8][1][1];
	ke[197]=K[8][0][1];
	ke[198]=K[6][5][1];
	ke[199]=K[6][4][1];
	ke[200]=K[6][3][1];
	ke[201]=K[6][2][1];
	ke[202]=K[6][1][1];
	ke[203]=K[6][0][1];
	ke[204]=K[3][5][1];
	ke[205]=K[3][4][1];
	ke[206]=K[3][3][1];
	ke[207]=K[3][2][1];
	ke[208]=K[3][1][1];
	ke[209]=K[3][0][1];
	ke[210]=K[9][2][2];
	ke[211]=K[9][1][2];
	ke[212]=K[9][0][2];
	ke[213]=K[8][5][2];
	ke[214]=K[8][4][2];
	ke[215]=K[8][3][2];
	ke[216]=K[8][2][2];
	ke[217]=K[8][1][2];
	ke[218]=K[8][0][2];
	ke[219]=K[6][5][2];
	ke[220]=K[6][4][2];
	ke[221]=K[6][3][2];
	ke[222]=K[6][2][2];
	ke[223]=K[6][1][2];
	ke[224]=K[6][0][2];
	ke[225]=K[3][5][2];
	ke[226]=K[3][4][2];
	ke[227]=K[3][3][2];
	ke[228]=K[3][2][2];
	ke[229]=K[3][1][2];
	ke[230]=K[3][0][2];
	ke[231]=K[9][3][3];
	ke[232]=K[9][2][3];
	ke[233]=K[9][1][3];
	ke[234]=K[9][0][3];
	ke[235]=K[8][5][3];
	ke[236]=K[8][4][3];
	ke[237]=K[8][3][3];
	ke[238]=K[8][2][3];
	ke[239]=K[8][1][3];
	ke[240]=K[8][0][3];
	ke[241]=K[6][5][3];
	ke[242]=K[6][4][3];
	ke[243]=K[6][3][3];
	ke[244]=K[6][2][3];
	ke[245]=K[6][1][3];
	ke[246]=K[6][0][3];
	ke[247]=K[3][5][3];
	ke[248]=K[3][4][3];
	ke[249]=K[3][3][3];
	ke[250]=K[3][2][3];
	ke[251]=K[3][1][3];
	ke[252]=K[3][0][3];
	ke[253]=K[9][4][4];
	ke[254]=K[9][3][4];
	ke[255]=K[9][2][4];
	ke[256]=K[9][1][4];
	ke[257]=K[9][0][4];
	ke[258]=K[8][5][4];
	ke[259]=K[8][4][4];
	ke[260]=K[8][3][4];
	ke[261]=K[8][2][4];
	ke[262]=K[8][1][4];
	ke[263]=K[8][0][4];
	ke[264]=K[6][5][4];
	ke[265]=K[6][4][4];
	ke[266]=K[6][3][4];
	ke[267]=K[6][2][4];
	ke[268]=K[6][1][4];
	ke[269]=K[6][0][4];
	ke[270]=K[3][5][4];
	ke[271]=K[3][4][4];
	ke[272]=K[3][3][4];
	ke[273]=K[3][2][4];
	ke[274]=K[3][1][4];
	ke[275]=K[3][0][4];
	ke[276]=K[9][5][5];
	ke[277]=K[9][4][5];
	ke[278]=K[9][3][5];
	ke[279]=K[9][2][5];
	ke[280]=K[9][1][5];
	ke[281]=K[9][0][5];
	ke[282]=K[8][5][5];
	ke[283]=K[8][4][5];
	ke[284]=K[8][3][5];
	ke[285]=K[8][2][5];
	ke[286]=K[8][1][5];
	ke[287]=K[8][0][5];
	ke[288]=K[6][5][5];
	ke[289]=K[6][4][5];
	ke[290]=K[6][3][5];
	ke[291]=K[6][2][5];
	ke[292]=K[6][1][5];
	ke[293]=K[6][0][5];
	ke[294]=K[3][5][5];
	ke[295]=K[3][4][5];
	ke[296]=K[3][3][5];
	ke[297]=K[3][2][5];
	ke[298]=K[3][1][5];
	ke[299]=K[3][0][5];
/*	The declaration of various variables shown before can also be woked out by the following code based on loop
	unsigned int num=0;
	for(unsigned int loop1=0;loop1<6;loop1++){
		for(unsigned int loop2=0;loop2<=loop1;loop2++){
			ke[num]=K[0][loop1-loop2][loop1];
			num++;
		}
	}
	for(unsigned int loop1=0;loop1<6;loop1++){
		for(unsigned int loop2=0;loop2<=loop1;loop2++){
			ke[num]=K[4][loop1-loop2][loop1];
			num++;
		}
		for(unsigned int loop2=0;loop2<6;loop2++){
			ke[num]=K[1][5-loop2][loop1];
			num++;
		}
	}
	for(unsigned int loop1=0;loop1<6;loop1++){
		for(unsigned int loop2=0;loop2<=loop1;loop2++){
			ke[num]=K[7][loop1-loop2][loop1];
			num++;
		}
		for(unsigned int loop2=0;loop2<6;loop2++){
			ke[num]=K[5][5-loop2][loop1];
			num++;
		}
		for(unsigned int loop2=0;loop2<6;loop2++){
			ke[num]=K[2][5-loop2][loop1];
			num++;
		}
	}
	for(unsigned int loop1=0;loop1<6;loop1++){
		for(unsigned int loop2=0;loop2<=loop1;loop2++){
			ke[num]=K[9][loop1-loop2][loop1];
			num++;
		}
		for(unsigned int loop2=0;loop2<6;loop2++){
			ke[num]=K[8][5-loop2][loop1];
			num++;
		}
		for(unsigned int loop2=0;loop2<6;loop2++){
			ke[num]=K[6][5-loop2][loop1];
			num++;
		}
		for(unsigned int loop2=0;loop2<6;loop2++){
			ke[num]=K[3][5-loop2][loop1];
			num++;
		}
	}*/
	double WiWjWkdetJ=fabs(detJ)*weight_gauss[i]*weight_gauss[j]*weight_gauss[k];
	for(unsigned int loop=0;loop<300;loop++){
		Matrix[loop]+=ke[loop]*WiWjWkdetJ;
	}

			}
		}
	}
}

void CShell::ElementMass(double* Mass){
	Mass[0]=0;
}

void CShell::ElementStress(double* stress, double* Displacement){
	clear(stress, 56);
	double X[4]={nodes_[0]->XYZ[0],nodes_[1]->XYZ[0],nodes_[2]->XYZ[0],nodes_[3]->XYZ[0]};
	double Y[4]={nodes_[0]->XYZ[1],nodes_[1]->XYZ[1],nodes_[2]->XYZ[1],nodes_[3]->XYZ[1]};
	double Z[4]={nodes_[0]->XYZ[2],nodes_[1]->XYZ[2],nodes_[2]->XYZ[2],nodes_[3]->XYZ[2]};

	//Calculate the unit normal vector for each node
	double l3[4], m3[4], n3[4];
	l3[0]=(Z[0]*(Y[0]-Y[1])-Z[0]*(Y[0]-Y[3])+Z[1]*(Y[0]-Y[3])-Z[3]*(Y[0]-Y[1]))/4;
	m3[0]=(X[0]*(Z[0]-Z[1])-X[0]*(Z[0]-Z[3])+X[1]*(Z[0]-Z[3])-X[3]*(Z[0]-Z[1]))/4;
	n3[0]=(Y[0]*(X[0]-X[1])-Y[0]*(X[0]-X[3])+Y[1]*(X[0]-X[3])-Y[3]*(X[0]-X[1]))/4;
	l3[1]=(Z[1]*(Y[0]-Y[1])-Z[0]*(Y[1]-Y[2])-Z[2]*(Y[0]-Y[1])+Z[1]*(Y[1]-Y[2]))/4;
	m3[1]=(X[1]*(Z[0]-Z[1])-X[0]*(Z[1]-Z[2])-X[2]*(Z[0]-Z[1])+X[1]*(Z[1]-Z[2]))/4;
	n3[1]=(Y[1]*(X[0]-X[1])-Y[0]*(X[1]-X[2])-Y[2]*(X[0]-X[1])+Y[1]*(X[1]-X[2]))/4;
	l3[2]=(Z[2]*(Y[1]-Y[2])-Z[1]*(Y[2]-Y[3])-Z[3]*(Y[1]-Y[2])+Z[2]*(Y[2]-Y[3]))/4;
	m3[2]=(X[2]*(Z[1]-Z[2])-X[1]*(Z[2]-Z[3])-X[3]*(Z[1]-Z[2])+X[2]*(Z[2]-Z[3]))/4;
	n3[2]=(Y[2]*(X[1]-X[2])-Y[1]*(X[2]-X[3])-Y[3]*(X[1]-X[2])+Y[2]*(X[2]-X[3]))/4;
	l3[3]=(Z[2]*(Y[0]-Y[3])-Z[0]*(Y[2]-Y[3])-Z[3]*(Y[0]-Y[3])+Z[3]*(Y[2]-Y[3]))/4;
	m3[3]=(X[2]*(Z[0]-Z[3])-X[0]*(Z[2]-Z[3])-X[3]*(Z[0]-Z[3])+X[3]*(Z[2]-Z[3]))/4;
	n3[3]=(Y[2]*(X[0]-X[3])-Y[0]*(X[2]-X[3])-Y[3]*(X[0]-X[3])+Y[3]*(X[2]-X[3]))/4;
	for(unsigned int loop=0;loop<4;loop++){
		double lmn3=sqrt(l3[loop]*l3[loop]+m3[loop]*m3[loop]+n3[loop]*n3[loop]);
		l3[loop]/=lmn3;
		m3[loop]/=lmn3;
		n3[loop]/=lmn3;
	}
	double l1[4], m1[4], n1[4];
	double l2[4], m2[4], n2[4];
	for(unsigned int loop=0;loop<4;loop++){
		double lmn1=sqrt(n3[loop]*n3[loop]+l3[loop]*l3[loop]);
		l1[loop]=n3[loop]/lmn1;
		m1[loop]=0;
		n1[loop]=-l3[loop]/lmn1;
		l2[loop]=-m3[loop]*l3[loop]/lmn1;
		m2[loop]=(n3[loop]*n3[loop]+l3[loop]*l3[loop])/lmn1;
		n2[loop]=-m3[loop]*n3[loop]/lmn1;
	}

	//Calculate magnititude of non-zero elements of constitutive matrix
	//D=[D1,D2,0,0,0;D2,D1,0,0,0;0,0,D3,0,0;0,0,0,D3,0;0,0,0,0,D3], dimension=5*5
	CShellMaterial* material_ = dynamic_cast<CShellMaterial*>(ElementMaterial_);
	double E_=material_->E;
	double nu_=material_->nu;
	double D1=E_/(1-nu_*nu_);
	double D2=D1*nu_;
	double D3=E_/2/(1+nu_);

	double thick=material_->thick;

	double gausspoint[2]={-0.57735027, 0.57735027};	//Gauss point when ngp=2
	double weight_gauss[2]={1,1};	//Gauss weight

	double Disp_element[24];
	for(int i=0;i<24;i++){
		if(!LocationMatrix_[i])
			Disp_element[i]=0;
		else
			Disp_element[i]=Displacement[LocationMatrix_[i]-1];
	}

	unsigned int num=0;		//stress[num],num=0:7,	8 gauss point, each with 6 stress component

	for(unsigned int i=0;i<2;i++){
		double xi=gausspoint[i];
		for(unsigned int j=0;j<2;j++){
			double eta=gausspoint[j];
			for(unsigned int k=0;k<2;k++){
				double zeta=gausspoint[k];

	double N[4]={(1-xi)*(1-eta)/4,(1+xi)*(1-eta)/4,(1+xi)*(1+eta)/4,(1-xi)*(1+eta)/4};
	double GNxi[4]={-(1-eta)/4,(1-eta)/4,(1+eta)/4,-(1+eta)/4};		//GNxi=dN/dxi
	double GNeta[4]={-(1-xi)/4,-(1+xi)/4,(1+xi)/4,(1-xi)/4};		//GNeta=dN/deta
	double J[3][3]={0,0,0,0,0,0,0,0,0};			//Calculate the Jocabian Matrix
	for(unsigned int loop=0;loop<4;loop++){
	J[0][0]+=GNxi[loop]*(X[loop]+0.5*zeta*thick*l3[loop]);
	J[0][1]+=GNxi[loop]*(Y[loop]+0.5*zeta*thick*m3[loop]);
	J[0][2]+=GNxi[loop]*(Z[loop]+0.5*zeta*thick*n3[loop]);
	J[1][0]+=GNeta[loop]*(X[loop]+0.5*zeta*thick*l3[loop]);
	J[1][1]+=GNeta[loop]*(Y[loop]+0.5*zeta*thick*m3[loop]);
	J[1][2]+=GNeta[loop]*(Z[loop]+0.5*zeta*thick*n3[loop]);
	J[2][0]+=N[loop]*0.5*thick*l3[loop];
	J[2][1]+=N[loop]*0.5*thick*m3[loop];
	J[2][2]+=N[loop]*0.5*thick*n3[loop];
	}
	double detJ=J[0][0]*J[1][1]*J[2][2]-J[0][0]*J[1][2]*J[2][1]-J[0][1]*J[1][0]*J[2][2]+J[0][1]*J[1][2]*J[2][0]+J[0][2]*J[1][0]*J[2][1]-J[0][2]*J[1][1]*J[2][0];
	double invJ[3][3];		//Calculate the inverse of the Jocabian Matrix
	invJ[0][0]=(J[1][1]*J[2][2]-J[1][2]*J[2][1])/detJ;
	invJ[0][1]=(J[0][2]*J[2][1]-J[0][1]*J[2][2])/detJ;
	invJ[0][2]=(J[0][1]*J[1][2]-J[0][2]*J[1][1])/detJ;
	invJ[1][0]=(J[1][2]*J[2][0]-J[1][0]*J[2][2])/detJ;
	invJ[1][1]=(J[0][0]*J[2][2]-J[0][2]*J[2][0])/detJ;
	invJ[1][2]=(J[0][2]*J[1][0]-J[0][0]*J[1][2])/detJ;
	invJ[2][0]=(J[1][0]*J[2][1]-J[1][1]*J[2][0])/detJ;
	invJ[2][1]=(J[0][1]*J[2][0]-J[0][0]*J[2][1])/detJ;
	invJ[2][2]=(J[0][0]*J[1][1]-J[0][1]*J[1][0])/detJ;
	double A[4][5][5];
	for(unsigned int loop=0;loop<4;loop++){
		A[loop][0][0]=GNeta[loop]*invJ[0][1]+GNxi[loop]*invJ[0][0];
		A[loop][0][1]=0.0;
		A[loop][0][2]=0.0;
		A[loop][0][3]=-N[loop]*invJ[0][2]*l2[loop]*thick/2-GNeta[loop]*invJ[0][1]*l2[loop]*thick/2*zeta-GNxi[loop]*invJ[0][0]*l2[loop]*thick/2*zeta;
		A[loop][0][4]=N[loop]*invJ[0][2]*l1[loop]*thick/2+GNeta[loop]*invJ[0][1]*l1[loop]*thick/2*zeta+GNxi[loop]*invJ[0][0]*l1[loop]*thick/2*zeta;
		A[loop][1][0]=0.0;
		A[loop][1][1]=GNeta[loop]*invJ[1][1]+GNxi[loop]*invJ[1][0];
		A[loop][1][2]=0.0;
		A[loop][1][3]=-N[loop]*invJ[1][2]*m2[loop]*thick/2-GNeta[loop]*invJ[1][1]*m2[loop]*thick/2*zeta-GNxi[loop]*invJ[1][0]*m2[loop]*thick/2*zeta;
		A[loop][1][4]=N[loop]*invJ[1][2]*m1[loop]*thick/2+GNeta[loop]*invJ[1][1]*m1[loop]*thick/2*zeta+GNxi[loop]*invJ[1][0]*m1[loop]*thick/2*zeta;
		A[loop][2][0]=GNeta[loop]*invJ[1][1]+GNxi[loop]*invJ[1][0];
		A[loop][2][1]=GNeta[loop]*invJ[0][1]+GNxi[loop]*invJ[0][0];
		A[loop][2][2]=0.0;
		A[loop][2][3]=-N[loop]*invJ[1][2]*l2[loop]*thick/2-N[loop]*invJ[0][2]*m2[loop]*thick/2-GNeta[loop]*invJ[1][1]*l2[loop]*thick/2*zeta-GNxi[loop]*invJ[1][0]*l2[loop]*thick/2*zeta-GNeta[loop]*invJ[0][1]*m2[loop]*thick/2*zeta-GNxi[loop]*invJ[0][0]*m2[loop]*thick/2*zeta;
		A[loop][2][4]=N[loop]*invJ[1][2]*l1[loop]*thick/2+N[loop]*invJ[0][2]*m1[loop]*thick/2+GNeta[loop]*invJ[1][1]*l1[loop]*thick/2*zeta+GNxi[loop]*invJ[1][0]*l1[loop]*thick/2*zeta+GNeta[loop]*invJ[0][1]*m1[loop]*thick/2*zeta+GNxi[loop]*invJ[0][0]*m1[loop]*thick/2*zeta;
		A[loop][3][0]=0.0;
		A[loop][3][1]=GNeta[loop]*invJ[2][1]+GNxi[loop]*invJ[2][0];
		A[loop][3][2]=GNeta[loop]*invJ[1][1]+GNxi[loop]*invJ[1][0];
		A[loop][3][3]=-N[loop]*invJ[2][2]*m2[loop]*thick/2-N[loop]*invJ[1][2]*n2[loop]*thick/2-GNeta[loop]*invJ[2][1]*m2[loop]*thick/2*zeta-GNxi[loop]*invJ[2][0]*m2[loop]*thick/2*zeta-GNeta[loop]*invJ[1][1]*n2[loop]*thick/2*zeta-GNxi[loop]*invJ[1][0]*n2[loop]*thick/2*zeta;
		A[loop][3][4]=N[loop]*invJ[2][2]*m1[loop]*thick/2+N[loop]*invJ[1][2]*n1[loop]*thick/2+GNeta[loop]*invJ[2][1]*m1[loop]*thick/2*zeta+GNxi[loop]*invJ[2][0]*m1[loop]*thick/2*zeta+GNeta[loop]*invJ[1][1]*n1[loop]*thick/2*zeta+GNxi[loop]*invJ[1][0]*n1[loop]*thick/2*zeta;
		A[loop][4][0]=GNeta[loop]*invJ[2][1]+GNxi[loop]*invJ[2][0];
		A[loop][4][1]=0.0;
		A[loop][4][2]=GNeta[loop]*invJ[0][1]+GNxi[loop]*invJ[0][0];
		A[loop][4][3]=-N[loop]*invJ[2][2]*l2[loop]*thick/2-N[loop]*invJ[0][2]*n2[loop]*thick/2-GNeta[loop]*invJ[2][1]*l2[loop]*thick/2*zeta-GNxi[loop]*invJ[2][0]*l2[loop]*thick/2*zeta-GNeta[loop]*invJ[0][1]*n2[loop]*thick/2*zeta-GNxi[loop]*invJ[0][0]*n2[loop]*thick/2*zeta;
		A[loop][4][4]=N[loop]*invJ[2][2]*l1[loop]*thick/2+N[loop]*invJ[0][2]*n1[loop]*thick/2+GNeta[loop]*invJ[2][1]*l1[loop]*thick/2*zeta+GNxi[loop]*invJ[2][0]*l1[loop]*thick/2*zeta+GNeta[loop]*invJ[0][1]*n1[loop]*thick/2*zeta+GNxi[loop]*invJ[0][0]*n1[loop]*thick/2*zeta;
	}

	double Be[4][5][6];
	for(unsigned int loop=0;loop<4;loop++){
		double A1=A[loop][0][0];
		double A2=A[loop][0][3];
		double A3=A[loop][0][4];
		double A4=A[loop][1][1];
		double A5=A[loop][1][3];
		double A6=A[loop][1][4];
		double A7=A[loop][2][0];
		double A8=A[loop][2][1];
		double A9=A[loop][2][3];
		double A10=A[loop][2][4];
		double A11=A[loop][3][1];
		double A12=A[loop][3][2];
		double A13=A[loop][3][3];
		double A14=A[loop][3][4];
		double A15=A[loop][4][0];
		double A16=A[loop][4][2];
		double A17=A[loop][4][3];
		double A18=A[loop][4][4];
		//Trans=[1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,ts1,ts2,ts3;0,0,0,ts4,ts5,ts6]
		double ts1=m2[loop]*n3[loop]-m3[loop]*n2[loop];
		double ts2=l3[loop]*n2[loop]-l2[loop]*n3[loop];
		double ts3=l2[loop]*m3[loop]-l3[loop]*m2[loop];
		double ts4=m3[loop]*n1[loop]-m1[loop]*n3[loop];
		double ts5=l1[loop]*n3[loop]-l3[loop]*n1[loop];
		double ts6=l3[loop]*m1[loop]-l1[loop]*m3[loop];

		Be[loop][0][0]=A1;
		Be[loop][0][1]=0.0;
		Be[loop][0][2]=0.0;
		Be[loop][0][3]=A2*ts1 + A3*ts4;
		Be[loop][0][4]=A2*ts2 + A3*ts5;
		Be[loop][0][5]=A2*ts3 + A3*ts6;
		Be[loop][1][0]=0.0;
		Be[loop][1][1]=A4;
		Be[loop][1][2]=0.0;
		Be[loop][1][3]=A5*ts1 + A6*ts4;
		Be[loop][1][4]=A5*ts2 + A6*ts5;
		Be[loop][1][5]=A5*ts3 + A6*ts6;
		Be[loop][2][0]=A7;
		Be[loop][2][1]=A8;
		Be[loop][2][2]=0.0;
		Be[loop][2][3]=A9*ts1 + A10*ts4;
		Be[loop][2][4]=A9*ts2 + A10*ts5;
		Be[loop][2][5]=A9*ts3 + A10*ts6;
		Be[loop][3][0]=0.0;
		Be[loop][3][1]=A11;
		Be[loop][3][2]=A12;
		Be[loop][3][3]=A13*ts1 + A14*ts4;
		Be[loop][3][4]=A13*ts2 + A14*ts5;
		Be[loop][3][5]=A13*ts3 + A14*ts6;
		Be[loop][4][0]=A15;
		Be[loop][4][1]=0.0;
		Be[loop][4][2]=A16;
		Be[loop][4][3]=A17*ts1 + A18*ts4;
		Be[loop][4][4]=A17*ts2 + A18*ts5;
		Be[loop][4][5]=A17*ts3 + A18*ts6;
	}

	double DBe[4][5][6];
	for(unsigned int loop=0;loop<4;loop++){
			//D*Be(loop), dimension of Be=5*6
			double P1=Be[loop][0][0];
			double P2=Be[loop][0][3];
			double P3=Be[loop][0][4];
			double P4=Be[loop][0][5];
			double P5=Be[loop][1][1];
			double P6=Be[loop][1][3];
			double P7=Be[loop][1][4];
			double P8=Be[loop][1][5];
			double P9=Be[loop][2][0];
			double P10=Be[loop][2][1];
			double P11=Be[loop][2][3];
			double P12=Be[loop][2][4];
			double P13=Be[loop][2][5];
			double P14=Be[loop][3][1];
			double P15=Be[loop][3][2];
			double P16=Be[loop][3][3];
			double P17=Be[loop][3][4];
			double P18=Be[loop][3][5];
			double P19=Be[loop][4][0];
			double P20=Be[loop][4][2];
			double P21=Be[loop][4][3];
			double P22=Be[loop][4][4];
			double P23=Be[loop][4][5];

			DBe[loop][0][0]=D1*P1;
			DBe[loop][0][1]=D2*P5;
			DBe[loop][0][2]=0.0;
			DBe[loop][0][3]=D1*P2 + D2*P6;
			DBe[loop][0][4]=D1*P3 + D2*P7;
			DBe[loop][0][5]=D1*P4 + D2*P8;
			DBe[loop][1][0]=D2*P1;
			DBe[loop][1][1]=D1*P5;
			DBe[loop][1][2]=0.0;
			DBe[loop][1][3]=D2*P2 + D1*P6;
			DBe[loop][1][4]=D2*P3 + D1*P7;
			DBe[loop][1][5]=D2*P4 + D1*P8;
			DBe[loop][2][0]=D3*P9;
			DBe[loop][2][1]=D3*P10;
			DBe[loop][2][2]=0.0;
			DBe[loop][2][3]=D3*P11;
			DBe[loop][2][4]=D3*P12;
			DBe[loop][2][5]=D3*P13;
			DBe[loop][3][0]=0.0;
			DBe[loop][3][1]=D3*P14;
			DBe[loop][3][2]=D3*P15;
			DBe[loop][3][3]=D3*P16;
			DBe[loop][3][4]=D3*P17;
			DBe[loop][3][5]=D3*P18;
			DBe[loop][4][0]=D3*P19;
			DBe[loop][4][1]=0.0;
			DBe[loop][4][2]=D3*P20;
			DBe[loop][4][3]=D3*P21;
			DBe[loop][4][4]=D3*P22;
			DBe[loop][4][5]=D3*P23;
	}

	stress[7*num+2]=0.0;
	for(unsigned int loop2=0;loop2<4;loop2++){
		for(unsigned int loop3=0;loop3<6;loop3++){
			stress[7*num]+=DBe[loop2][0][loop3]*Disp_element[6*loop2+loop3];
			stress[7*num+1]+=DBe[loop2][1][loop3]*Disp_element[6*loop2+loop3];
			stress[7*num+3]+=DBe[loop2][3][loop3]*Disp_element[6*loop2+loop3];
			stress[7*num+4]+=DBe[loop2][4][loop3]*Disp_element[6*loop2+loop3];
			stress[7*num+5]+=DBe[loop2][5][loop3]*Disp_element[6*loop2+loop3];
		}
		stress[7*num+6]=sqrt(stress[7*num]*stress[7*num]-stress[7*num]*stress[7*num+1]+stress[7*num+1]*stress[7*num+1]+3*(stress[7*num+3]*stress[7*num+3]+stress[7*num+4]*stress[7*num+4]+stress[7*num+5]*stress[7*num+5]));
	}
	num++;
			}
		}
	}
}

void CShell::ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* PostPositions){
	;
}

void CShell::GravityCalculation(double* ptr_force){
	double g = 9.8;
	CShellMaterial* material_ = dynamic_cast<CShellMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double density = material_->density;
	double thick=material_->thick;
	double X[4]={nodes_[0]->XYZ[0],nodes_[1]->XYZ[0],nodes_[2]->XYZ[0],nodes_[3]->XYZ[0]};
	double Y[4]={nodes_[0]->XYZ[1],nodes_[1]->XYZ[1],nodes_[2]->XYZ[1],nodes_[3]->XYZ[1]};
	double Z[4]={nodes_[0]->XYZ[2],nodes_[1]->XYZ[2],nodes_[2]->XYZ[2],nodes_[3]->XYZ[2]};

	//Calculate the unit normal vector for each node
	double l3[4], m3[4], n3[4];
	l3[0]=(Z[0]*(Y[0]-Y[1])-Z[0]*(Y[0]-Y[3])+Z[1]*(Y[0]-Y[3])-Z[3]*(Y[0]-Y[1]))/4;
	m3[0]=(X[0]*(Z[0]-Z[1])-X[0]*(Z[0]-Z[3])+X[1]*(Z[0]-Z[3])-X[3]*(Z[0]-Z[1]))/4;
	n3[0]=(Y[0]*(X[0]-X[1])-Y[0]*(X[0]-X[3])+Y[1]*(X[0]-X[3])-Y[3]*(X[0]-X[1]))/4;
	l3[1]=(Z[1]*(Y[0]-Y[1])-Z[0]*(Y[1]-Y[2])-Z[2]*(Y[0]-Y[1])+Z[1]*(Y[1]-Y[2]))/4;
	m3[1]=(X[1]*(Z[0]-Z[1])-X[0]*(Z[1]-Z[2])-X[2]*(Z[0]-Z[1])+X[1]*(Z[1]-Z[2]))/4;
	n3[1]=(Y[1]*(X[0]-X[1])-Y[0]*(X[1]-X[2])-Y[2]*(X[0]-X[1])+Y[1]*(X[1]-X[2]))/4;
	l3[2]=(Z[2]*(Y[1]-Y[2])-Z[1]*(Y[2]-Y[3])-Z[3]*(Y[1]-Y[2])+Z[2]*(Y[2]-Y[3]))/4;
	m3[2]=(X[2]*(Z[1]-Z[2])-X[1]*(Z[2]-Z[3])-X[3]*(Z[1]-Z[2])+X[2]*(Z[2]-Z[3]))/4;
	n3[2]=(Y[2]*(X[1]-X[2])-Y[1]*(X[2]-X[3])-Y[3]*(X[1]-X[2])+Y[2]*(X[2]-X[3]))/4;
	l3[3]=(Z[2]*(Y[0]-Y[3])-Z[0]*(Y[2]-Y[3])-Z[3]*(Y[0]-Y[3])+Z[3]*(Y[2]-Y[3]))/4;
	m3[3]=(X[2]*(Z[0]-Z[3])-X[0]*(Z[2]-Z[3])-X[3]*(Z[0]-Z[3])+X[3]*(Z[2]-Z[3]))/4;
	n3[3]=(Y[2]*(X[0]-X[3])-Y[0]*(X[2]-X[3])-Y[3]*(X[0]-X[3])+Y[3]*(X[2]-X[3]))/4;
	for(unsigned int loop=0;loop<4;loop++){
		double lmn3=sqrt(l3[loop]*l3[loop]+m3[loop]*m3[loop]+n3[loop]*n3[loop]);
		l3[loop]/=lmn3;
		m3[loop]/=lmn3;
		n3[loop]/=lmn3;
	}
	double l1[4], m1[4], n1[4];
	double l2[4], m2[4], n2[4];
	for(unsigned int loop=0;loop<4;loop++){
	double lmn1=sqrt(n3[loop]*n3[loop]+l3[loop]*l3[loop]);
	l1[loop]=n3[loop]/lmn1;
	m1[loop]=0;
	n1[loop]=-l3[loop]/lmn1;
	l2[loop]=-m3[loop]*l3[loop]/lmn1;
	m2[loop]=(n3[loop]*n3[loop]+l3[loop]*l3[loop])/lmn1;
	n2[loop]=-m3[loop]*n3[loop]/lmn1;
	}

	double gausspoint[2]={-0.57735027, 0.57735027};	//Gauss point when ngp=2
	double weight_gauss[2]={1,1};	//Gauss weight

	double volumn=0;
	for(unsigned int i=0;i<2;i++){
		double xi=gausspoint[i];
		for(unsigned int j=0;j<2;j++){
			double eta=gausspoint[j];
			for(unsigned int k=0;k<2;k++){
				double zeta=gausspoint[k];
				double N[4]={(1-xi)*(1-eta)/4,(1+xi)*(1-eta)/4,(1+xi)*(1+eta)/4,(1-xi)*(1+eta)/4};
				double GNxi[4]={-(1-eta)/4,(1-eta)/4,(1+eta)/4,-(1+eta)/4};		//GNxi=dN/dxi
				double GNeta[4]={-(1-xi)/4,-(1+xi)/4,(1+xi)/4,(1-xi)/4};		//GNeta=dN/deta
				double J[3][3]={0,0,0,0,0,0,0,0,0};			//Calculate the Jocabian Matrix
			for(unsigned int loop=0;loop<4;loop++){
				J[0][0]+=GNxi[loop]*(X[loop]+0.5*zeta*thick*l3[loop]);
				J[0][1]+=GNxi[loop]*(Y[loop]+0.5*zeta*thick*m3[loop]);
				J[0][2]+=GNxi[loop]*(Z[loop]+0.5*zeta*thick*n3[loop]);
				J[1][0]+=GNeta[loop]*(X[loop]+0.5*zeta*thick*l3[loop]);
				J[1][1]+=GNeta[loop]*(Y[loop]+0.5*zeta*thick*m3[loop]);
				J[1][2]+=GNeta[loop]*(Z[loop]+0.5*zeta*thick*n3[loop]);
				J[2][0]+=N[loop]*0.5*thick*l3[loop];
				J[2][1]+=N[loop]*0.5*thick*m3[loop];
				J[2][2]+=N[loop]*0.5*thick*n3[loop];
			}
			double detJ=J[0][0]*J[1][1]*J[2][2]-J[0][0]*J[1][2]*J[2][1]-J[0][1]*J[1][0]*J[2][2]+J[0][1]*J[1][2]*J[2][0]+J[0][2]*J[1][0]*J[2][1]-J[0][2]*J[1][1]*J[2][0];
			volumn+=fabs(detJ)*weight_gauss[i]*weight_gauss[j]*weight_gauss[k];


			}
		}
	}
	weight=volumn*g*density;
	double weight_avg=weight/4;
	double Xc=(X[0]+X[1]+X[2]+X[3])/4;
	double Yc=(Y[0]+Y[1]+Y[2]+Y[3])/4;
	double Zc=(Z[0]+Z[1]+Z[2]+Z[3])/4;
	for(unsigned int loop;loop<4;loop++){
		ptr_force[3*loop]=-weight_avg;
		ptr_force[3*loop+1]=2*weight_avg/3*(Yc-Y[loop]);
		ptr_force[3*loop+2]=-2*weight_avg/3*(Xc-X[loop]);
	}

}

void CShell::ElementCoord(double* coord){
	clear(coord,24);
	CShellMaterial* material_ = dynamic_cast<CShellMaterial*>(ElementMaterial_);
	double thick=material_->thick;
	double X[4]={nodes_[0]->XYZ[0],nodes_[1]->XYZ[0],nodes_[2]->XYZ[0],nodes_[3]->XYZ[0]};
	double Y[4]={nodes_[0]->XYZ[1],nodes_[1]->XYZ[1],nodes_[2]->XYZ[1],nodes_[3]->XYZ[1]};
	double Z[4]={nodes_[0]->XYZ[2],nodes_[1]->XYZ[2],nodes_[2]->XYZ[2],nodes_[3]->XYZ[2]};

	//Calculate the unit normal vector for each node
	double l3[4], m3[4], n3[4];
	l3[0]=(Z[0]*(Y[0]-Y[1])-Z[0]*(Y[0]-Y[3])+Z[1]*(Y[0]-Y[3])-Z[3]*(Y[0]-Y[1]))/4;
	m3[0]=(X[0]*(Z[0]-Z[1])-X[0]*(Z[0]-Z[3])+X[1]*(Z[0]-Z[3])-X[3]*(Z[0]-Z[1]))/4;
	n3[0]=(Y[0]*(X[0]-X[1])-Y[0]*(X[0]-X[3])+Y[1]*(X[0]-X[3])-Y[3]*(X[0]-X[1]))/4;
	l3[1]=(Z[1]*(Y[0]-Y[1])-Z[0]*(Y[1]-Y[2])-Z[2]*(Y[0]-Y[1])+Z[1]*(Y[1]-Y[2]))/4;
	m3[1]=(X[1]*(Z[0]-Z[1])-X[0]*(Z[1]-Z[2])-X[2]*(Z[0]-Z[1])+X[1]*(Z[1]-Z[2]))/4;
	n3[1]=(Y[1]*(X[0]-X[1])-Y[0]*(X[1]-X[2])-Y[2]*(X[0]-X[1])+Y[1]*(X[1]-X[2]))/4;
	l3[2]=(Z[2]*(Y[1]-Y[2])-Z[1]*(Y[2]-Y[3])-Z[3]*(Y[1]-Y[2])+Z[2]*(Y[2]-Y[3]))/4;
	m3[2]=(X[2]*(Z[1]-Z[2])-X[1]*(Z[2]-Z[3])-X[3]*(Z[1]-Z[2])+X[2]*(Z[2]-Z[3]))/4;
	n3[2]=(Y[2]*(X[1]-X[2])-Y[1]*(X[2]-X[3])-Y[3]*(X[1]-X[2])+Y[2]*(X[2]-X[3]))/4;
	l3[3]=(Z[2]*(Y[0]-Y[3])-Z[0]*(Y[2]-Y[3])-Z[3]*(Y[0]-Y[3])+Z[3]*(Y[2]-Y[3]))/4;
	m3[3]=(X[2]*(Z[0]-Z[3])-X[0]*(Z[2]-Z[3])-X[3]*(Z[0]-Z[3])+X[3]*(Z[2]-Z[3]))/4;
	n3[3]=(Y[2]*(X[0]-X[3])-Y[0]*(X[2]-X[3])-Y[3]*(X[0]-X[3])+Y[3]*(X[2]-X[3]))/4;
	for(unsigned int loop=0;loop<4;loop++){
		double lmn3=sqrt(l3[loop]*l3[loop]+m3[loop]*m3[loop]+n3[loop]*n3[loop]);
		l3[loop]/=lmn3;
		m3[loop]/=lmn3;
		n3[loop]/=lmn3;
	}
	double l1[4], m1[4], n1[4];
	double l2[4], m2[4], n2[4];
	for(unsigned int loop=0;loop<4;loop++){
		double lmn1=sqrt(n3[loop]*n3[loop]+l3[loop]*l3[loop]);
		l1[loop]=n3[loop]/lmn1;
		m1[loop]=0;
		n1[loop]=-l3[loop]/lmn1;
		l2[loop]=-m3[loop]*l3[loop]/lmn1;
		m2[loop]=(n3[loop]*n3[loop]+l3[loop]*l3[loop])/lmn1;
		n2[loop]=-m3[loop]*n3[loop]/lmn1;
	}
	double gausspoint[2]={-0.57735027, 0.57735027};	//Gauss point when ngp=2
	unsigned int num=0;
	for(unsigned int i=0;i<2;i++){
		double xi=gausspoint[i];
		for(unsigned int j=0;j<2;j++){
			double eta=gausspoint[j];
			for(unsigned int k=0;k<2;k++){
				double zeta=gausspoint[k];
				double N[4]={(1-xi)*(1-eta)/4,(1+xi)*(1-eta)/4,(1+xi)*(1+eta)/4,(1-xi)*(1+eta)/4};
				coord[num++]=N[0]*(X[0]+zeta*0.5*thick*l3[0])+N[1]*(X[1]+zeta*0.5*thick*l3[1])+N[2]*(X[2]+zeta*0.5*thick*l3[2])+N[3]*(X[3]+zeta*0.5*thick*l3[3]);
				coord[num++]=N[0]*(Y[0]+zeta*0.5*thick*m3[0])+N[1]*(Y[1]+zeta*0.5*thick*m3[1])+N[2]*(Y[2]+zeta*0.5*thick*m3[2])+N[3]*(Y[3]+zeta*0.5*thick*m3[3]);
				coord[num++]=N[0]*(Z[0]+zeta*0.5*thick*n3[0])+N[1]*(Z[1]+zeta*0.5*thick*n3[1])+N[2]*(Z[2]+zeta*0.5*thick*n3[2])+N[3]*(Z[3]+zeta*0.5*thick*n3[3]);
			}
		}
	}


}


