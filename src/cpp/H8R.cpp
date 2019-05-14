/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "H8R.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CH8R::CH8R()
{
	NEN_ = 8;	// Each element has 8 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 24;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CH8R::~CH8R()
{
}

//	Read element data from stream Input
bool CH8R::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8 ;	// four node numbers

	Input >> N1 >> N2 >> N3 >> N4>> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial_ = dynamic_cast<CH8Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[4] = &NodeList[N5 - 1];
	nodes_[5] = &NodeList[N6 - 1];
	nodes_[6] = &NodeList[N7 - 1];
	nodes_[7] = &NodeList[N8 - 1];

	return true;
}

//	Write element data to stream
void CH8R::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber 
		   << setw(9) << nodes_[3]->NodeNumber << setw(9) << nodes_[4]->NodeNumber << 
		    setw(9) << nodes_[5]->NodeNumber << setw(9) << nodes_[6]->NodeNumber << 
		   setw(9) << nodes_[7]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

void CH8R::WritePlot(COutPlot& output, unsigned int Ele)
{
	output << 8 << setw(11) << nodes_[0]->NodeNumber-1 << setw(9) << nodes_[1]->NodeNumber-1 << setw(9) << nodes_[3]->NodeNumber-1 
		   << setw(9) << nodes_[2]->NodeNumber-1 << setw(9) << nodes_[4]->NodeNumber-1 << 
		    setw(9) << nodes_[5]->NodeNumber-1 << setw(9) << nodes_[7]->NodeNumber-1 << 
		   setw(9) << nodes_[6]->NodeNumber-1 << endl;
}

void CH8R::WritePlotPost(COutPlotPost& output, unsigned int Ele)
{
	output << 8 << setw(11) << nodes_[0]->NodeNumber-1 << setw(9) << nodes_[1]->NodeNumber-1 << setw(9) << nodes_[3]->NodeNumber-1 
		   << setw(9) << nodes_[2]->NodeNumber-1 << setw(9) << nodes_[4]->NodeNumber-1 << 
		    setw(9) << nodes_[5]->NodeNumber-1 << setw(9) << nodes_[7]->NodeNumber-1 << 
		   setw(9) << nodes_[6]->NodeNumber-1 << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CH8R::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
	
/*	printf("LM\n");
	for (int ii=0;ii<8;ii++)
		printf("%d\n",LocationMatrix_[ii]);
*/
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 8H element, element stiffness is a 24x24 matrix, whose upper triangular part
//	has 300 elements
unsigned int CH8R::SizeOfStiffnessMatrix() { return 300; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element

void CH8R::CalculateBe (double* Be, double psi, double eta, double zet)
{
	
//GN matrix
    double GN[12];
	
	GN[0] = (1-eta)*(1-zet)/8;
	GN[1] = (1-eta)*(1+zet)/8;
	GN[2] = (1+eta)*(1-zet)/8;
	GN[3] = (1+eta)*(1+zet)/8;
	GN[4] = (1-psi)*(1-zet)/8;
	GN[5] = (1-psi)*(1+zet)/8;
	GN[6] = (1+psi)*(1-zet)/8;
	GN[7] = (1+psi)*(1+zet)/8;
	GN[8] = (1-eta)*(1-psi)/8;
	GN[9] = (1-eta)*(1+psi)/8;
	GN[10] = (1+eta)*(1-psi)/8;
	GN[11] = (1+eta)*(1+psi)/8;
	
	
//Jaccobi
// [0	3	6]
// [1	4	7]
// [2	5	8]
	double Ja[9];
	Ja[0] = -GN[0]*nodes_[0]->XYZ[0] +GN[0]*nodes_[1]->XYZ[0] +GN[2]*nodes_[2]->XYZ[0] -GN[2]*nodes_[3]->XYZ[0]
		-GN[1]*nodes_[4]->XYZ[0] +GN[1]*nodes_[5]->XYZ[0] +GN[3]*nodes_[6]->XYZ[0] -GN[3]*nodes_[7]->XYZ[0];
	Ja[1] = -GN[4]*nodes_[0]->XYZ[0] -GN[6]*nodes_[1]->XYZ[0] +GN[6]*nodes_[2]->XYZ[0] +GN[4]*nodes_[3]->XYZ[0]
		-GN[5]*nodes_[4]->XYZ[0] -GN[7]*nodes_[5]->XYZ[0] +GN[7]*nodes_[6]->XYZ[0] +GN[5]*nodes_[7]->XYZ[0];
	Ja[2] = -GN[8]*nodes_[0]->XYZ[0] -GN[9]*nodes_[1]->XYZ[0] -GN[11]*nodes_[2]->XYZ[0] -GN[10]*nodes_[3]->XYZ[0]
		+GN[8]*nodes_[4]->XYZ[0] +GN[9]*nodes_[5]->XYZ[0] +GN[11]*nodes_[6]->XYZ[0] +GN[10]*nodes_[7]->XYZ[0];
	Ja[3] = -GN[0]*nodes_[0]->XYZ[1] +GN[0]*nodes_[1]->XYZ[1] +GN[2]*nodes_[2]->XYZ[1] -GN[2]*nodes_[3]->XYZ[1]
		-GN[1]*nodes_[4]->XYZ[1] +GN[1]*nodes_[5]->XYZ[1] +GN[3]*nodes_[6]->XYZ[1] -GN[3]*nodes_[7]->XYZ[1];
	Ja[4] = -GN[4]*nodes_[0]->XYZ[1] -GN[6]*nodes_[1]->XYZ[1] +GN[6]*nodes_[2]->XYZ[1] +GN[4]*nodes_[3]->XYZ[1]
		-GN[5]*nodes_[4]->XYZ[1] -GN[7]*nodes_[5]->XYZ[1] +GN[7]*nodes_[6]->XYZ[1] +GN[5]*nodes_[7]->XYZ[1];
	Ja[5] = -GN[8]*nodes_[0]->XYZ[1] -GN[9]*nodes_[1]->XYZ[1] -GN[11]*nodes_[2]->XYZ[1] -GN[10]*nodes_[3]->XYZ[1]
		+GN[8]*nodes_[4]->XYZ[1] +GN[9]*nodes_[5]->XYZ[1] +GN[11]*nodes_[6]->XYZ[1] +GN[10]*nodes_[7]->XYZ[1];
	Ja[6] = -GN[2]*nodes_[0]->XYZ[2] +GN[2]*nodes_[1]->XYZ[2] +GN[2]*nodes_[2]->XYZ[2] -GN[2]*nodes_[3]->XYZ[2]
		-GN[1]*nodes_[4]->XYZ[2] +GN[1]*nodes_[5]->XYZ[2] +GN[3]*nodes_[6]->XYZ[2] -GN[3]*nodes_[7]->XYZ[2];
	Ja[7] = -GN[4]*nodes_[0]->XYZ[2] -GN[6]*nodes_[1]->XYZ[2] +GN[6]*nodes_[2]->XYZ[2] +GN[4]*nodes_[3]->XYZ[2]
		-GN[5]*nodes_[4]->XYZ[2] -GN[7]*nodes_[5]->XYZ[2] +GN[7]*nodes_[6]->XYZ[2] +GN[5]*nodes_[7]->XYZ[2];
	Ja[8] = -GN[8]*nodes_[0]->XYZ[2] -GN[9]*nodes_[1]->XYZ[2] -GN[11]*nodes_[2]->XYZ[2] -GN[10]*nodes_[3]->XYZ[2]
		+GN[8]*nodes_[4]->XYZ[2] +GN[9]*nodes_[5]->XYZ[2] +GN[11]*nodes_[6]->XYZ[2] +GN[10]*nodes_[7]->XYZ[2];
	
	// printf("jac: \n");
	// for (int i=0;i<9;i++)
		// printf("%lf\n",Ja[i]);
	// printf("\n");
	
	double JaCof[9];
	JaCof[0] = Ja[4]*Ja[8]-Ja[5]*Ja[7];
	JaCof[1] = -Ja[3]*Ja[8]+Ja[5]*Ja[6];
	JaCof[2] = Ja[3]*Ja[7]-Ja[4]*Ja[6];
	JaCof[3] = -Ja[1]*Ja[8]+Ja[2]*Ja[7];
	JaCof[4] = Ja[0]*Ja[8]-Ja[2]*Ja[6];
	JaCof[5] = -Ja[0]*Ja[7]+Ja[1]*Ja[6];
	JaCof[6] = Ja[1]*Ja[5]-Ja[2]*Ja[4];
	JaCof[7] = -Ja[0]*Ja[5]+Ja[2]*Ja[3];
	JaCof[8] = Ja[0]*Ja[4]-Ja[1]*Ja[3];
	
	double JaDet;
	JaDet = Ja[0]*JaCof[0] + Ja[1]*JaCof[1] + Ja[2]*JaCof[2];
	
	// printf("jacdet: %lf\n\n",JaDet);
//inversion
	double JaInv[9];
	JaInv[0] = JaCof[0]/JaDet;
	JaInv[1] = JaCof[3]/JaDet;
	JaInv[2] = JaCof[6]/JaDet;
	JaInv[3] = JaCof[1]/JaDet;
	JaInv[4] = JaCof[4]/JaDet;
	JaInv[5] = JaCof[7]/JaDet;
	JaInv[6] = JaCof[2]/JaDet;
	JaInv[7] = JaCof[5]/JaDet;
	JaInv[8] = JaCof[8]/JaDet;
	

//Be [1x 1y 1z ... 8z JaDet] 25 menbers
	
	Be[0] = -JaInv[0]*GN[0] -JaInv[3]*GN[4] -JaInv[6]*GN[8];
	Be[3] =  JaInv[0]*GN[0] -JaInv[3]*GN[6] -JaInv[6]*GN[9];
	Be[6] =  JaInv[0]*GN[2] +JaInv[3]*GN[6] -JaInv[6]*GN[11];
	Be[9] = -JaInv[0]*GN[2] +JaInv[3]*GN[4] -JaInv[6]*GN[10];
	Be[12] = -JaInv[0]*GN[1] -JaInv[3]*GN[5] +JaInv[6]*GN[8];
	Be[15] =  JaInv[0]*GN[1] -JaInv[3]*GN[7] +JaInv[6]*GN[9];
	Be[18] =  JaInv[0]*GN[3] +JaInv[3]*GN[7] +JaInv[6]*GN[11];
	Be[21] = -JaInv[0]*GN[3] +JaInv[3]*GN[5] +JaInv[6]*GN[10];
	
	Be[1] = -JaInv[1]*GN[0] -JaInv[4]*GN[4] -JaInv[7]*GN[8];
	Be[4] =  JaInv[1]*GN[0] -JaInv[4]*GN[6] -JaInv[7]*GN[9];
	Be[7] =  JaInv[1]*GN[2] +JaInv[4]*GN[6] -JaInv[7]*GN[11];
	Be[10] = -JaInv[1]*GN[2] +JaInv[4]*GN[4] -JaInv[7]*GN[10];
	Be[13] = -JaInv[1]*GN[1] -JaInv[4]*GN[5] +JaInv[7]*GN[8];
	Be[16] =  JaInv[1]*GN[1] -JaInv[4]*GN[7] +JaInv[7]*GN[9];
	Be[19] =  JaInv[1]*GN[3] +JaInv[4]*GN[7] +JaInv[7]*GN[11];
	Be[22] = -JaInv[1]*GN[3] +JaInv[4]*GN[5] +JaInv[7]*GN[10];
	
	Be[2] = -JaInv[2]*GN[0] -JaInv[5]*GN[4] -JaInv[8]*GN[8];
	Be[5] =  JaInv[2]*GN[0] -JaInv[5]*GN[6] -JaInv[8]*GN[9];
	Be[8] =  JaInv[2]*GN[2] +JaInv[5]*GN[6] -JaInv[8]*GN[11];
	Be[11] = -JaInv[2]*GN[2] +JaInv[5]*GN[4] -JaInv[8]*GN[10];
	Be[14] = -JaInv[2]*GN[1] -JaInv[5]*GN[5] +JaInv[8]*GN[8];
	Be[17] =  JaInv[2]*GN[1] -JaInv[5]*GN[7] +JaInv[8]*GN[9];
	Be[20] =  JaInv[2]*GN[3] +JaInv[5]*GN[7] +JaInv[8]*GN[11];
	Be[23] = -JaInv[2]*GN[3] +JaInv[5]*GN[5] +JaInv[8]*GN[10];
	
	Be[24] = JaDet;

/*	
	 printf("Be\n");
	 for (int i=0;i<24;i++)
		 printf("%lf\n",Be[i]);
	 printf("\n");
*/
}

	
void CH8R::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	double Be[25];

	CalculateBe(Be,0,0,0);
	CalculateMatrix(Be, Matrix);
	
	for (unsigned int i=0;i<SizeOfStiffnessMatrix();i++)
		Matrix[i] = Matrix[i]*8;

}

//	Calculate element stress 
void CH8R::ElementStress(double* stress, double* Displacement)
{
	double Be[25];
	CalculateBe(Be,0,0,0);
	CalculateStress(Be, stress, Displacement);

}


void CH8R::ElementCoord (double* coord)
{
	CalculateCoord (coord, 0, 0, 0);
	
}


void CH8R::CalculateCoord (double* coord, double psi, double eta, double zet)
{
	double N[8];
	N[0] = (1-psi)*(1-eta)*(1-zet)/8;
	N[1] = (1+psi)*(1-eta)*(1-zet)/8;
	N[2] = (1+psi)*(1+eta)*(1-zet)/8;
	N[3] = (1-psi)*(1+eta)*(1-zet)/8;
	N[4] = (1-psi)*(1-eta)*(1+zet)/8;
	N[5] = (1+psi)*(1-eta)*(1+zet)/8;
	N[6] = (1+psi)*(1+eta)*(1+zet)/8;
	N[7] = (1-psi)*(1+eta)*(1+zet)/8;
	
	coord[0] = nodes_[0]->XYZ[0]*N[0] +nodes_[1]->XYZ[0]*N[1] +nodes_[2]->XYZ[0]*N[2] +nodes_[3]->XYZ[0]*N[3] +nodes_[4]->XYZ[0]*N[4] +nodes_[5]->XYZ[0]*N[5] +nodes_[6]->XYZ[0]*N[6] +nodes_[7]->XYZ[0]*N[7];
	coord[1] = nodes_[0]->XYZ[1]*N[0] +nodes_[1]->XYZ[1]*N[1] +nodes_[2]->XYZ[1]*N[2] +nodes_[3]->XYZ[1]*N[3] +nodes_[4]->XYZ[1]*N[4] +nodes_[5]->XYZ[1]*N[5] +nodes_[6]->XYZ[1]*N[6] +nodes_[7]->XYZ[1]*N[7];
	coord[2] = nodes_[0]->XYZ[2]*N[0] +nodes_[1]->XYZ[2]*N[1] +nodes_[2]->XYZ[2]*N[2] +nodes_[3]->XYZ[2]*N[3] +nodes_[4]->XYZ[2]*N[4] +nodes_[5]->XYZ[2]*N[5] +nodes_[6]->XYZ[2]*N[6] +nodes_[7]->XYZ[2]*N[7];
	
}
	
void CH8R::CalculateStress (double* Be, double* stress, double* Displacement)
{
	CH8Material* material_ = dynamic_cast<CH8Material*>(ElementMaterial_);
	double CM[3];
	CM[1]= material_->Lam;
	CM[2]= material_->G;
	CM[0]= 2*CM[2] + CM[1];
	double strain[6];

	strain[0] = Be[0]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) +Be[3]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[6]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0) +Be[9]*Displacement[LocationMatrix_[9]-1]*(double)(!LocationMatrix_[9]==0) +Be[12]*Displacement[LocationMatrix_[12]-1]*(double)(!LocationMatrix_[12]==0) +Be[15]*Displacement[LocationMatrix_[15]-1]*(double)(!LocationMatrix_[15]==0) +Be[18]*Displacement[LocationMatrix_[18]-1]*(double)(!LocationMatrix_[18]==0) +Be[21]*Displacement[LocationMatrix_[21]-1]*(double)(!LocationMatrix_[21]==0);
	strain[1] = Be[1]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) +Be[4]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) +Be[7]*Displacement[LocationMatrix_[7]-1]*(double)(!LocationMatrix_[7]==0) +Be[10]*Displacement[LocationMatrix_[10]-1]*(double)(!LocationMatrix_[10]==0) +Be[13]*Displacement[LocationMatrix_[13]-1]*(double)(!LocationMatrix_[13]==0) +Be[16]*Displacement[LocationMatrix_[16]-1]*(double)(!LocationMatrix_[16]==0) +Be[19]*Displacement[LocationMatrix_[19]-1]*(double)(!LocationMatrix_[19]==0) +Be[22]*Displacement[LocationMatrix_[22]-1]*(double)(!LocationMatrix_[22]==0);
	strain[2] = Be[2]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[5]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) +Be[8]*Displacement[LocationMatrix_[8]-1]*(double)(!LocationMatrix_[8]==0) +Be[11]*Displacement[LocationMatrix_[11]-1]*(double)(!LocationMatrix_[11]==0) +Be[14]*Displacement[LocationMatrix_[14]-1]*(double)(!LocationMatrix_[14]==0) +Be[17]*Displacement[LocationMatrix_[17]-1]*(double)(!LocationMatrix_[17]==0) +Be[20]*Displacement[LocationMatrix_[20]-1]*(double)(!LocationMatrix_[20]==0) +Be[23]*Displacement[LocationMatrix_[23]-1]*(double)(!LocationMatrix_[23]==0);
	strain[3] = Be[2]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) +Be[1]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[5]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) +Be[4]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) +Be[8]*Displacement[LocationMatrix_[7]-1]*(double)(!LocationMatrix_[7]==0) +Be[7]*Displacement[LocationMatrix_[8]-1]*(double)(!LocationMatrix_[8]==0) +Be[11]*Displacement[LocationMatrix_[10]-1]*(double)(!LocationMatrix_[10]==0) +Be[10]*Displacement[LocationMatrix_[11]-1]*(double)(!LocationMatrix_[11]==0) +Be[14]*Displacement[LocationMatrix_[13]-1]*(double)(!LocationMatrix_[13]==0) +Be[13]*Displacement[LocationMatrix_[14]-1]*(double)(!LocationMatrix_[14]==0) +Be[17]*Displacement[LocationMatrix_[16]-1]*(double)(!LocationMatrix_[16]==0) +Be[16]*Displacement[LocationMatrix_[17]-1]*(double)(!LocationMatrix_[17]==0) +Be[20]*Displacement[LocationMatrix_[19]-1]*(double)(!LocationMatrix_[19]==0) +Be[19]*Displacement[LocationMatrix_[20]-1]*(double)(!LocationMatrix_[20]==0) +Be[23]*Displacement[LocationMatrix_[22]-1]*(double)(!LocationMatrix_[22]==0) +Be[22]*Displacement[LocationMatrix_[23]-1]*(double)(!LocationMatrix_[23]==0);
	strain[4] = Be[2]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) +Be[0]*Displacement[LocationMatrix_[2]-1]*(double)(!LocationMatrix_[2]==0) +Be[5]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[3]*Displacement[LocationMatrix_[5]-1]*(double)(!LocationMatrix_[5]==0) +Be[8]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0) +Be[6]*Displacement[LocationMatrix_[8]-1]*(double)(!LocationMatrix_[8]==0) +Be[11]*Displacement[LocationMatrix_[9]-1]*(double)(!LocationMatrix_[9]==0) +Be[9]*Displacement[LocationMatrix_[11]-1]*(double)(!LocationMatrix_[11]==0) +Be[14]*Displacement[LocationMatrix_[12]-1]*(double)(!LocationMatrix_[12]==0) +Be[12]*Displacement[LocationMatrix_[14]-1]*(double)(!LocationMatrix_[14]==0) +Be[17]*Displacement[LocationMatrix_[15]-1]*(double)(!LocationMatrix_[15]==0) +Be[15]*Displacement[LocationMatrix_[17]-1]*(double)(!LocationMatrix_[17]==0) +Be[20]*Displacement[LocationMatrix_[18]-1]*(double)(!LocationMatrix_[18]==0) +Be[18]*Displacement[LocationMatrix_[20]-1]*(double)(!LocationMatrix_[20]==0) +Be[23]*Displacement[LocationMatrix_[21]-1]*(double)(!LocationMatrix_[21]==0) +Be[21]*Displacement[LocationMatrix_[23]-1]*(double)(!LocationMatrix_[23]==0);
	strain[5] = Be[1]*Displacement[LocationMatrix_[0]-1]*(double)(!LocationMatrix_[0]==0) +Be[0]*Displacement[LocationMatrix_[1]-1]*(double)(!LocationMatrix_[1]==0) +Be[4]*Displacement[LocationMatrix_[3]-1]*(double)(!LocationMatrix_[3]==0) +Be[3]*Displacement[LocationMatrix_[4]-1]*(double)(!LocationMatrix_[4]==0) +Be[7]*Displacement[LocationMatrix_[6]-1]*(double)(!LocationMatrix_[6]==0) +Be[6]*Displacement[LocationMatrix_[7]-1]*(double)(!LocationMatrix_[7]==0) +Be[10]*Displacement[LocationMatrix_[9]-1]*(double)(!LocationMatrix_[9]==0) +Be[9]*Displacement[LocationMatrix_[10]-1]*(double)(!LocationMatrix_[10]==0) +Be[13]*Displacement[LocationMatrix_[12]-1]*(double)(!LocationMatrix_[12]==0) +Be[12]*Displacement[LocationMatrix_[13]-1]*(double)(!LocationMatrix_[13]==0) +Be[16]*Displacement[LocationMatrix_[15]-1]*(double)(!LocationMatrix_[15]==0) +Be[15]*Displacement[LocationMatrix_[16]-1]*(double)(!LocationMatrix_[16]==0) +Be[19]*Displacement[LocationMatrix_[18]-1]*(double)(!LocationMatrix_[18]==0) +Be[18]*Displacement[LocationMatrix_[19]-1]*(double)(!LocationMatrix_[19]==0) +Be[22]*Displacement[LocationMatrix_[21]-1]*(double)(!LocationMatrix_[21]==0) +Be[21]*Displacement[LocationMatrix_[22]-1]*(double)(!LocationMatrix_[22]==0);

	stress[0] = CM[0]*strain[0] +CM[1]*strain[1] +CM[1]*strain[2];
	stress[1] = CM[1]*strain[0] +CM[0]*strain[1] +CM[1]*strain[2];
	stress[2] = CM[1]*strain[0] +CM[1]*strain[1] +CM[0]*strain[2];
	stress[3] = CM[2]*strain[3];
	stress[4] = CM[2]*strain[4];
	stress[5] = CM[2]*strain[5];

	
}

void CH8R::CalculateVolume (double* EG, double psi, double eta, double zet)
{
//GN matrix
    double GN[12];
	
	GN[0] = (1-eta)*(1-zet)/8;
	GN[1] = (1-eta)*(1+zet)/8;
	GN[2] = (1+eta)*(1-zet)/8;
	GN[3] = (1+eta)*(1+zet)/8;
	GN[4] = (1-psi)*(1-zet)/8;
	GN[5] = (1-psi)*(1+zet)/8;
	GN[6] = (1+psi)*(1-zet)/8;
	GN[7] = (1+psi)*(1+zet)/8;
	GN[8] = (1-eta)*(1-psi)/8;
	GN[9] = (1-eta)*(1+psi)/8;
	GN[10] = (1+eta)*(1-psi)/8;
	GN[11] = (1+eta)*(1+psi)/8;
	
	
//Jaccobi
// [0	3	6]
// [1	4	7]
// [2	5	8]
	double Ja[9];
	Ja[0] = -GN[0]*nodes_[0]->XYZ[0] +GN[0]*nodes_[1]->XYZ[0] +GN[2]*nodes_[2]->XYZ[0] -GN[2]*nodes_[3]->XYZ[0]
		-GN[1]*nodes_[4]->XYZ[0] +GN[1]*nodes_[5]->XYZ[0] +GN[3]*nodes_[6]->XYZ[0] -GN[3]*nodes_[7]->XYZ[0];
	Ja[1] = -GN[4]*nodes_[0]->XYZ[0] -GN[6]*nodes_[1]->XYZ[0] +GN[6]*nodes_[2]->XYZ[0] +GN[4]*nodes_[3]->XYZ[0]
		-GN[5]*nodes_[4]->XYZ[0] -GN[7]*nodes_[5]->XYZ[0] +GN[7]*nodes_[6]->XYZ[0] +GN[5]*nodes_[7]->XYZ[0];
	Ja[2] = -GN[8]*nodes_[0]->XYZ[0] -GN[9]*nodes_[1]->XYZ[0] -GN[11]*nodes_[2]->XYZ[0] -GN[10]*nodes_[3]->XYZ[0]
		+GN[8]*nodes_[4]->XYZ[0] +GN[9]*nodes_[5]->XYZ[0] +GN[11]*nodes_[6]->XYZ[0] +GN[10]*nodes_[7]->XYZ[0];
	Ja[3] = -GN[0]*nodes_[0]->XYZ[1] +GN[0]*nodes_[1]->XYZ[1] +GN[2]*nodes_[2]->XYZ[1] -GN[2]*nodes_[3]->XYZ[1]
		-GN[1]*nodes_[4]->XYZ[1] +GN[1]*nodes_[5]->XYZ[1] +GN[3]*nodes_[6]->XYZ[1] -GN[3]*nodes_[7]->XYZ[1];
	Ja[4] = -GN[4]*nodes_[0]->XYZ[1] -GN[6]*nodes_[1]->XYZ[1] +GN[6]*nodes_[2]->XYZ[1] +GN[4]*nodes_[3]->XYZ[1]
		-GN[5]*nodes_[4]->XYZ[1] -GN[7]*nodes_[5]->XYZ[1] +GN[7]*nodes_[6]->XYZ[1] +GN[5]*nodes_[7]->XYZ[1];
	Ja[5] = -GN[8]*nodes_[0]->XYZ[1] -GN[9]*nodes_[1]->XYZ[1] -GN[11]*nodes_[2]->XYZ[1] -GN[10]*nodes_[3]->XYZ[1]
		+GN[8]*nodes_[4]->XYZ[1] +GN[9]*nodes_[5]->XYZ[1] +GN[11]*nodes_[6]->XYZ[1] +GN[10]*nodes_[7]->XYZ[1];
	Ja[6] = -GN[2]*nodes_[0]->XYZ[2] +GN[2]*nodes_[1]->XYZ[2] +GN[2]*nodes_[2]->XYZ[2] -GN[2]*nodes_[3]->XYZ[2]
		-GN[1]*nodes_[4]->XYZ[2] +GN[1]*nodes_[5]->XYZ[2] +GN[3]*nodes_[6]->XYZ[2] -GN[3]*nodes_[7]->XYZ[2];
	Ja[7] = -GN[4]*nodes_[0]->XYZ[2] -GN[6]*nodes_[1]->XYZ[2] +GN[6]*nodes_[2]->XYZ[2] +GN[4]*nodes_[3]->XYZ[2]
		-GN[5]*nodes_[4]->XYZ[2] -GN[7]*nodes_[5]->XYZ[2] +GN[7]*nodes_[6]->XYZ[2] +GN[5]*nodes_[7]->XYZ[2];
	Ja[8] = -GN[8]*nodes_[0]->XYZ[2] -GN[9]*nodes_[1]->XYZ[2] -GN[11]*nodes_[2]->XYZ[2] -GN[10]*nodes_[3]->XYZ[2]
		+GN[8]*nodes_[4]->XYZ[2] +GN[9]*nodes_[5]->XYZ[2] +GN[11]*nodes_[6]->XYZ[2] +GN[10]*nodes_[7]->XYZ[2];
	
	
	double JaDet;
	JaDet = Ja[0]*(Ja[4]*Ja[8]-Ja[5]*Ja[7]) + Ja[1]*(-Ja[3]*Ja[8]+Ja[5]*Ja[6]) + Ja[2]*(Ja[3]*Ja[7]-Ja[4]*Ja[6]);
	
	double N[8];
	N[0] = (1-psi)*(1-eta)*(1-zet)/8;
	N[1] = (1+psi)*(1-eta)*(1-zet)/8;
	N[2] = (1+psi)*(1+eta)*(1-zet)/8;
	N[3] = (1-psi)*(1+eta)*(1-zet)/8;
	N[4] = (1-psi)*(1-eta)*(1+zet)/8;
	N[5] = (1+psi)*(1-eta)*(1+zet)/8;
	N[6] = (1+psi)*(1+eta)*(1+zet)/8;
	N[7] = (1-psi)*(1+eta)*(1+zet)/8;
	
	for (unsigned int i=0;i<NEN_;i++)
		EG[i] = EG[i] + N[i]*JaDet;
	
}
	
void CH8R::GravityCalculation (double* EG)
{
	clear(EG, NEN_);
	double Ga=0.57735;
	
	CalculateVolume(EG, -Ga, -Ga, -Ga);
	CalculateVolume(EG, Ga, -Ga, -Ga);
	CalculateVolume(EG, Ga, Ga, -Ga);
	CalculateVolume(EG, -Ga, Ga, -Ga);
	CalculateVolume(EG, -Ga, -Ga, Ga);
	CalculateVolume(EG, Ga, -Ga, Ga);
	CalculateVolume(EG, Ga, Ga, Ga);
	CalculateVolume(EG, -Ga, Ga, Ga);
	
	CH8Material* material_ = dynamic_cast<CH8Material*>(ElementMaterial_);
	
	for (unsigned int  i=0;i<NEN_;i++)
//		EG[i]= 0;
		EG[i]= EG[i]*material_->Rou*9.8;
}

/*
void CH8R::ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* PostPositions)
{
	for (unsigned int i=0;i<NEN_;i++)
	{
		for (unsigned int j=0;j<3;j++)
		{
			PrePositions[3*i+j] = nodes_[i] ->XYZ[0];
			if (LocationMatrix_[3*i+j])
				PostPositions[3*i+j] = PrePositions[3*i+j] + Displacement[LocationMatrix_[3*i+j]-1];
			else
				PostPositions[3*i+j] = PrePositions[3*i+j];
		}	
	}
	ElementStress(stress, Displacement);
}
*/

void CH8R::CalculateMatrix (double* Be, double* Matrix)
{
	CH8Material* material_ = dynamic_cast<CH8Material*>(ElementMaterial_);
	double CM[3];
	CM[1]= material_->Lam;
	CM[2]= material_->G;
	CM[0]= 2*CM[2] + CM[1];
Matrix[  0] = Matrix[0] + ( Be[0]*CM[0]* Be[0]  + Be[1]*CM[2]* Be[1]  + Be[2]*CM[2]* Be[2] )*Be[24];
Matrix[  1] = Matrix[1] + ( Be[0]*CM[2]* Be[0]  + Be[1]*CM[0]* Be[1]  + Be[2]*CM[2]* Be[2] )*Be[24];
Matrix[  2] = Matrix[2] + ( Be[1]*CM[1]* Be[0]  + Be[0]*CM[2]* Be[1] )*Be[24];
Matrix[  3] = Matrix[3] + ( Be[0]*CM[2]* Be[0]  + Be[1]*CM[2]* Be[1]  + Be[2]*CM[0]* Be[2] )*Be[24];
Matrix[  4] = Matrix[4] + ( Be[2]*CM[1]* Be[1]  + Be[1]*CM[2]* Be[2] )*Be[24];
Matrix[  5] = Matrix[5] + ( Be[2]*CM[1]* Be[0]  + Be[0]*CM[2]* Be[2] )*Be[24];
Matrix[  6] = Matrix[6] + ( Be[3]*CM[0]* Be[3]  + Be[4]*CM[2]* Be[4]  + Be[5]*CM[2]* Be[5] )*Be[24];
Matrix[  7] = Matrix[7] + ( Be[3]*CM[1]* Be[2]  + Be[5]*CM[2]* Be[0] )*Be[24];
Matrix[  8] = Matrix[8] + ( Be[3]*CM[1]* Be[1]  + Be[4]*CM[2]* Be[0] )*Be[24];
Matrix[  9] = Matrix[9] + ( Be[3]*CM[0]* Be[0]  + Be[4]*CM[2]* Be[1]  + Be[5]*CM[2]* Be[2] )*Be[24];
Matrix[ 10] = Matrix[10] + ( Be[3]*CM[2]* Be[3]  + Be[4]*CM[0]* Be[4]  + Be[5]*CM[2]* Be[5] )*Be[24];
Matrix[ 11] = Matrix[11] + ( Be[4]*CM[1]* Be[3]  + Be[3]*CM[2]* Be[4] )*Be[24];
Matrix[ 12] = Matrix[12] + ( Be[4]*CM[1]* Be[2]  + Be[5]*CM[2]* Be[1] )*Be[24];
Matrix[ 13] = Matrix[13] + ( Be[3]*CM[2]* Be[0]  + Be[4]*CM[0]* Be[1]  + Be[5]*CM[2]* Be[2] )*Be[24];
Matrix[ 14] = Matrix[14] + ( Be[4]*CM[1]* Be[0]  + Be[3]*CM[2]* Be[1] )*Be[24];
Matrix[ 15] = Matrix[15] + ( Be[3]*CM[2]* Be[3]  + Be[4]*CM[2]* Be[4]  + Be[5]*CM[0]* Be[5] )*Be[24];
Matrix[ 16] = Matrix[16] + ( Be[5]*CM[1]* Be[4]  + Be[4]*CM[2]* Be[5] )*Be[24];
Matrix[ 17] = Matrix[17] + ( Be[5]*CM[1]* Be[3]  + Be[3]*CM[2]* Be[5] )*Be[24];
Matrix[ 18] = Matrix[18] + ( Be[3]*CM[2]* Be[0]  + Be[4]*CM[2]* Be[1]  + Be[5]*CM[0]* Be[2] )*Be[24];
Matrix[ 19] = Matrix[19] + ( Be[5]*CM[1]* Be[1]  + Be[4]*CM[2]* Be[2] )*Be[24];
Matrix[ 20] = Matrix[20] + ( Be[5]*CM[1]* Be[0]  + Be[3]*CM[2]* Be[2] )*Be[24];
Matrix[ 21] = Matrix[21] + ( Be[6]*CM[0]* Be[6]  + Be[7]*CM[2]* Be[7]  + Be[8]*CM[2]* Be[8] )*Be[24];
Matrix[ 22] = Matrix[22] + ( Be[6]*CM[1]* Be[5]  + Be[8]*CM[2]* Be[3] )*Be[24];
Matrix[ 23] = Matrix[23] + ( Be[6]*CM[1]* Be[4]  + Be[7]*CM[2]* Be[3] )*Be[24];
Matrix[ 24] = Matrix[24] + ( Be[6]*CM[0]* Be[3]  + Be[7]*CM[2]* Be[4]  + Be[8]*CM[2]* Be[5] )*Be[24];
Matrix[ 25] = Matrix[25] + ( Be[6]*CM[1]* Be[2]  + Be[8]*CM[2]* Be[0] )*Be[24];
Matrix[ 26] = Matrix[26] + ( Be[6]*CM[1]* Be[1]  + Be[7]*CM[2]* Be[0] )*Be[24];
Matrix[ 27] = Matrix[27] + ( Be[6]*CM[0]* Be[0]  + Be[7]*CM[2]* Be[1]  + Be[8]*CM[2]* Be[2] )*Be[24];
Matrix[ 28] = Matrix[28] + ( Be[6]*CM[2]* Be[6]  + Be[7]*CM[0]* Be[7]  + Be[8]*CM[2]* Be[8] )*Be[24];
Matrix[ 29] = Matrix[29] + ( Be[7]*CM[1]* Be[6]  + Be[6]*CM[2]* Be[7] )*Be[24];
Matrix[ 30] = Matrix[30] + ( Be[7]*CM[1]* Be[5]  + Be[8]*CM[2]* Be[4] )*Be[24];
Matrix[ 31] = Matrix[31] + ( Be[6]*CM[2]* Be[3]  + Be[7]*CM[0]* Be[4]  + Be[8]*CM[2]* Be[5] )*Be[24];
Matrix[ 32] = Matrix[32] + ( Be[7]*CM[1]* Be[3]  + Be[6]*CM[2]* Be[4] )*Be[24];
Matrix[ 33] = Matrix[33] + ( Be[7]*CM[1]* Be[2]  + Be[8]*CM[2]* Be[1] )*Be[24];
Matrix[ 34] = Matrix[34] + ( Be[6]*CM[2]* Be[0]  + Be[7]*CM[0]* Be[1]  + Be[8]*CM[2]* Be[2] )*Be[24];
Matrix[ 35] = Matrix[35] + ( Be[7]*CM[1]* Be[0]  + Be[6]*CM[2]* Be[1] )*Be[24];
Matrix[ 36] = Matrix[36] + ( Be[6]*CM[2]* Be[6]  + Be[7]*CM[2]* Be[7]  + Be[8]*CM[0]* Be[8] )*Be[24];
Matrix[ 37] = Matrix[37] + ( Be[8]*CM[1]* Be[7]  + Be[7]*CM[2]* Be[8] )*Be[24];
Matrix[ 38] = Matrix[38] + ( Be[8]*CM[1]* Be[6]  + Be[6]*CM[2]* Be[8] )*Be[24];
Matrix[ 39] = Matrix[39] + ( Be[6]*CM[2]* Be[3]  + Be[7]*CM[2]* Be[4]  + Be[8]*CM[0]* Be[5] )*Be[24];
Matrix[ 40] = Matrix[40] + ( Be[8]*CM[1]* Be[4]  + Be[7]*CM[2]* Be[5] )*Be[24];
Matrix[ 41] = Matrix[41] + ( Be[8]*CM[1]* Be[3]  + Be[6]*CM[2]* Be[5] )*Be[24];
Matrix[ 42] = Matrix[42] + ( Be[6]*CM[2]* Be[0]  + Be[7]*CM[2]* Be[1]  + Be[8]*CM[0]* Be[2] )*Be[24];
Matrix[ 43] = Matrix[43] + ( Be[8]*CM[1]* Be[1]  + Be[7]*CM[2]* Be[2] )*Be[24];
Matrix[ 44] = Matrix[44] + ( Be[8]*CM[1]* Be[0]  + Be[6]*CM[2]* Be[2] )*Be[24];
Matrix[ 45] = Matrix[45] + ( Be[9]*CM[0]* Be[9]  + Be[10]*CM[2]* Be[10]  + Be[11]*CM[2]* Be[11] )*Be[24];
Matrix[ 46] = Matrix[46] + ( Be[9]*CM[1]* Be[8]  + Be[11]*CM[2]* Be[6] )*Be[24];
Matrix[ 47] = Matrix[47] + ( Be[9]*CM[1]* Be[7]  + Be[10]*CM[2]* Be[6] )*Be[24];
Matrix[ 48] = Matrix[48] + ( Be[9]*CM[0]* Be[6]  + Be[10]*CM[2]* Be[7]  + Be[11]*CM[2]* Be[8] )*Be[24];
Matrix[ 49] = Matrix[49] + ( Be[9]*CM[1]* Be[5]  + Be[11]*CM[2]* Be[3] )*Be[24];
Matrix[ 50] = Matrix[50] + ( Be[9]*CM[1]* Be[4]  + Be[10]*CM[2]* Be[3] )*Be[24];
Matrix[ 51] = Matrix[51] + ( Be[9]*CM[0]* Be[3]  + Be[10]*CM[2]* Be[4]  + Be[11]*CM[2]* Be[5] )*Be[24];
Matrix[ 52] = Matrix[52] + ( Be[9]*CM[1]* Be[2]  + Be[11]*CM[2]* Be[0] )*Be[24];
Matrix[ 53] = Matrix[53] + ( Be[9]*CM[1]* Be[1]  + Be[10]*CM[2]* Be[0] )*Be[24];
Matrix[ 54] = Matrix[54] + ( Be[9]*CM[0]* Be[0]  + Be[10]*CM[2]* Be[1]  + Be[11]*CM[2]* Be[2] )*Be[24];
Matrix[ 55] = Matrix[55] + ( Be[9]*CM[2]* Be[9]  + Be[10]*CM[0]* Be[10]  + Be[11]*CM[2]* Be[11] )*Be[24];
Matrix[ 56] = Matrix[56] + ( Be[10]*CM[1]* Be[9]  + Be[9]*CM[2]* Be[10] )*Be[24];
Matrix[ 57] = Matrix[57] + ( Be[10]*CM[1]* Be[8]  + Be[11]*CM[2]* Be[7] )*Be[24];
Matrix[ 58] = Matrix[58] + ( Be[9]*CM[2]* Be[6]  + Be[10]*CM[0]* Be[7]  + Be[11]*CM[2]* Be[8] )*Be[24];
Matrix[ 59] = Matrix[59] + ( Be[10]*CM[1]* Be[6]  + Be[9]*CM[2]* Be[7] )*Be[24];
Matrix[ 60] = Matrix[60] + ( Be[10]*CM[1]* Be[5]  + Be[11]*CM[2]* Be[4] )*Be[24];
Matrix[ 61] = Matrix[61] + ( Be[9]*CM[2]* Be[3]  + Be[10]*CM[0]* Be[4]  + Be[11]*CM[2]* Be[5] )*Be[24];
Matrix[ 62] = Matrix[62] + ( Be[10]*CM[1]* Be[3]  + Be[9]*CM[2]* Be[4] )*Be[24];
Matrix[ 63] = Matrix[63] + ( Be[10]*CM[1]* Be[2]  + Be[11]*CM[2]* Be[1] )*Be[24];
Matrix[ 64] = Matrix[64] + ( Be[9]*CM[2]* Be[0]  + Be[10]*CM[0]* Be[1]  + Be[11]*CM[2]* Be[2] )*Be[24];
Matrix[ 65] = Matrix[65] + ( Be[10]*CM[1]* Be[0]  + Be[9]*CM[2]* Be[1] )*Be[24];
Matrix[ 66] = Matrix[66] + ( Be[9]*CM[2]* Be[9]  + Be[10]*CM[2]* Be[10]  + Be[11]*CM[0]* Be[11] )*Be[24];
Matrix[ 67] = Matrix[67] + ( Be[11]*CM[1]* Be[10]  + Be[10]*CM[2]* Be[11] )*Be[24];
Matrix[ 68] = Matrix[68] + ( Be[11]*CM[1]* Be[9]  + Be[9]*CM[2]* Be[11] )*Be[24];
Matrix[ 69] = Matrix[69] + ( Be[9]*CM[2]* Be[6]  + Be[10]*CM[2]* Be[7]  + Be[11]*CM[0]* Be[8] )*Be[24];
Matrix[ 70] = Matrix[70] + ( Be[11]*CM[1]* Be[7]  + Be[10]*CM[2]* Be[8] )*Be[24];
Matrix[ 71] = Matrix[71] + ( Be[11]*CM[1]* Be[6]  + Be[9]*CM[2]* Be[8] )*Be[24];
Matrix[ 72] = Matrix[72] + ( Be[9]*CM[2]* Be[3]  + Be[10]*CM[2]* Be[4]  + Be[11]*CM[0]* Be[5] )*Be[24];
Matrix[ 73] = Matrix[73] + ( Be[11]*CM[1]* Be[4]  + Be[10]*CM[2]* Be[5] )*Be[24];
Matrix[ 74] = Matrix[74] + ( Be[11]*CM[1]* Be[3]  + Be[9]*CM[2]* Be[5] )*Be[24];
Matrix[ 75] = Matrix[75] + ( Be[9]*CM[2]* Be[0]  + Be[10]*CM[2]* Be[1]  + Be[11]*CM[0]* Be[2] )*Be[24];
Matrix[ 76] = Matrix[76] + ( Be[11]*CM[1]* Be[1]  + Be[10]*CM[2]* Be[2] )*Be[24];
Matrix[ 77] = Matrix[77] + ( Be[11]*CM[1]* Be[0]  + Be[9]*CM[2]* Be[2] )*Be[24];
Matrix[ 78] = Matrix[78] + ( Be[12]*CM[0]* Be[12]  + Be[13]*CM[2]* Be[13]  + Be[14]*CM[2]* Be[14] )*Be[24];
Matrix[ 79] = Matrix[79] + ( Be[12]*CM[1]* Be[11]  + Be[14]*CM[2]* Be[9] )*Be[24];
Matrix[ 80] = Matrix[80] + ( Be[12]*CM[1]* Be[10]  + Be[13]*CM[2]* Be[9] )*Be[24];
Matrix[ 81] = Matrix[81] + ( Be[12]*CM[0]* Be[9]  + Be[13]*CM[2]* Be[10]  + Be[14]*CM[2]* Be[11] )*Be[24];
Matrix[ 82] = Matrix[82] + ( Be[12]*CM[1]* Be[8]  + Be[14]*CM[2]* Be[6] )*Be[24];
Matrix[ 83] = Matrix[83] + ( Be[12]*CM[1]* Be[7]  + Be[13]*CM[2]* Be[6] )*Be[24];
Matrix[ 84] = Matrix[84] + ( Be[12]*CM[0]* Be[6]  + Be[13]*CM[2]* Be[7]  + Be[14]*CM[2]* Be[8] )*Be[24];
Matrix[ 85] = Matrix[85] + ( Be[12]*CM[1]* Be[5]  + Be[14]*CM[2]* Be[3] )*Be[24];
Matrix[ 86] = Matrix[86] + ( Be[12]*CM[1]* Be[4]  + Be[13]*CM[2]* Be[3] )*Be[24];
Matrix[ 87] = Matrix[87] + ( Be[12]*CM[0]* Be[3]  + Be[13]*CM[2]* Be[4]  + Be[14]*CM[2]* Be[5] )*Be[24];
Matrix[ 88] = Matrix[88] + ( Be[12]*CM[1]* Be[2]  + Be[14]*CM[2]* Be[0] )*Be[24];
Matrix[ 89] = Matrix[89] + ( Be[12]*CM[1]* Be[1]  + Be[13]*CM[2]* Be[0] )*Be[24];
Matrix[ 90] = Matrix[90] + ( Be[12]*CM[0]* Be[0]  + Be[13]*CM[2]* Be[1]  + Be[14]*CM[2]* Be[2] )*Be[24];
Matrix[ 91] = Matrix[91] + ( Be[12]*CM[2]* Be[12]  + Be[13]*CM[0]* Be[13]  + Be[14]*CM[2]* Be[14] )*Be[24];
Matrix[ 92] = Matrix[92] + ( Be[13]*CM[1]* Be[12]  + Be[12]*CM[2]* Be[13] )*Be[24];
Matrix[ 93] = Matrix[93] + ( Be[13]*CM[1]* Be[11]  + Be[14]*CM[2]* Be[10] )*Be[24];
Matrix[ 94] = Matrix[94] + ( Be[12]*CM[2]* Be[9]  + Be[13]*CM[0]* Be[10]  + Be[14]*CM[2]* Be[11] )*Be[24];
Matrix[ 95] = Matrix[95] + ( Be[13]*CM[1]* Be[9]  + Be[12]*CM[2]* Be[10] )*Be[24];
Matrix[ 96] = Matrix[96] + ( Be[13]*CM[1]* Be[8]  + Be[14]*CM[2]* Be[7] )*Be[24];
Matrix[ 97] = Matrix[97] + ( Be[12]*CM[2]* Be[6]  + Be[13]*CM[0]* Be[7]  + Be[14]*CM[2]* Be[8] )*Be[24];
Matrix[ 98] = Matrix[98] + ( Be[13]*CM[1]* Be[6]  + Be[12]*CM[2]* Be[7] )*Be[24];
Matrix[ 99] = Matrix[99] + ( Be[13]*CM[1]* Be[5]  + Be[14]*CM[2]* Be[4] )*Be[24];
Matrix[100] = Matrix[100] + ( Be[12]*CM[2]* Be[3]  + Be[13]*CM[0]* Be[4]  + Be[14]*CM[2]* Be[5] )*Be[24];
Matrix[101] = Matrix[101] + ( Be[13]*CM[1]* Be[3]  + Be[12]*CM[2]* Be[4] )*Be[24];
Matrix[102] = Matrix[102] + ( Be[13]*CM[1]* Be[2]  + Be[14]*CM[2]* Be[1] )*Be[24];
Matrix[103] = Matrix[103] + ( Be[12]*CM[2]* Be[0]  + Be[13]*CM[0]* Be[1]  + Be[14]*CM[2]* Be[2] )*Be[24];
Matrix[104] = Matrix[104] + ( Be[13]*CM[1]* Be[0]  + Be[12]*CM[2]* Be[1] )*Be[24];
Matrix[105] = Matrix[105] + ( Be[12]*CM[2]* Be[12]  + Be[13]*CM[2]* Be[13]  + Be[14]*CM[0]* Be[14] )*Be[24];
Matrix[106] = Matrix[106] + ( Be[14]*CM[1]* Be[13]  + Be[13]*CM[2]* Be[14] )*Be[24];
Matrix[107] = Matrix[107] + ( Be[14]*CM[1]* Be[12]  + Be[12]*CM[2]* Be[14] )*Be[24];
Matrix[108] = Matrix[108] + ( Be[12]*CM[2]* Be[9]  + Be[13]*CM[2]* Be[10]  + Be[14]*CM[0]* Be[11] )*Be[24];
Matrix[109] = Matrix[109] + ( Be[14]*CM[1]* Be[10]  + Be[13]*CM[2]* Be[11] )*Be[24];
Matrix[110] = Matrix[110] + ( Be[14]*CM[1]* Be[9]  + Be[12]*CM[2]* Be[11] )*Be[24];
Matrix[111] = Matrix[111] + ( Be[12]*CM[2]* Be[6]  + Be[13]*CM[2]* Be[7]  + Be[14]*CM[0]* Be[8] )*Be[24];
Matrix[112] = Matrix[112] + ( Be[14]*CM[1]* Be[7]  + Be[13]*CM[2]* Be[8] )*Be[24];
Matrix[113] = Matrix[113] + ( Be[14]*CM[1]* Be[6]  + Be[12]*CM[2]* Be[8] )*Be[24];
Matrix[114] = Matrix[114] + ( Be[12]*CM[2]* Be[3]  + Be[13]*CM[2]* Be[4]  + Be[14]*CM[0]* Be[5] )*Be[24];
Matrix[115] = Matrix[115] + ( Be[14]*CM[1]* Be[4]  + Be[13]*CM[2]* Be[5] )*Be[24];
Matrix[116] = Matrix[116] + ( Be[14]*CM[1]* Be[3]  + Be[12]*CM[2]* Be[5] )*Be[24];
Matrix[117] = Matrix[117] + ( Be[12]*CM[2]* Be[0]  + Be[13]*CM[2]* Be[1]  + Be[14]*CM[0]* Be[2] )*Be[24];
Matrix[118] = Matrix[118] + ( Be[14]*CM[1]* Be[1]  + Be[13]*CM[2]* Be[2] )*Be[24];
Matrix[119] = Matrix[119] + ( Be[14]*CM[1]* Be[0]  + Be[12]*CM[2]* Be[2] )*Be[24];
Matrix[120] = Matrix[120] + ( Be[15]*CM[0]* Be[15]  + Be[16]*CM[2]* Be[16]  + Be[17]*CM[2]* Be[17] )*Be[24];
Matrix[121] = Matrix[121] + ( Be[15]*CM[1]* Be[14]  + Be[17]*CM[2]* Be[12] )*Be[24];
Matrix[122] = Matrix[122] + ( Be[15]*CM[1]* Be[13]  + Be[16]*CM[2]* Be[12] )*Be[24];
Matrix[123] = Matrix[123] + ( Be[15]*CM[0]* Be[12]  + Be[16]*CM[2]* Be[13]  + Be[17]*CM[2]* Be[14] )*Be[24];
Matrix[124] = Matrix[124] + ( Be[15]*CM[1]* Be[11]  + Be[17]*CM[2]* Be[9] )*Be[24];
Matrix[125] = Matrix[125] + ( Be[15]*CM[1]* Be[10]  + Be[16]*CM[2]* Be[9] )*Be[24];
Matrix[126] = Matrix[126] + ( Be[15]*CM[0]* Be[9]  + Be[16]*CM[2]* Be[10]  + Be[17]*CM[2]* Be[11] )*Be[24];
Matrix[127] = Matrix[127] + ( Be[15]*CM[1]* Be[8]  + Be[17]*CM[2]* Be[6] )*Be[24];
Matrix[128] = Matrix[128] + ( Be[15]*CM[1]* Be[7]  + Be[16]*CM[2]* Be[6] )*Be[24];
Matrix[129] = Matrix[129] + ( Be[15]*CM[0]* Be[6]  + Be[16]*CM[2]* Be[7]  + Be[17]*CM[2]* Be[8] )*Be[24];
Matrix[130] = Matrix[130] + ( Be[15]*CM[1]* Be[5]  + Be[17]*CM[2]* Be[3] )*Be[24];
Matrix[131] = Matrix[131] + ( Be[15]*CM[1]* Be[4]  + Be[16]*CM[2]* Be[3] )*Be[24];
Matrix[132] = Matrix[132] + ( Be[15]*CM[0]* Be[3]  + Be[16]*CM[2]* Be[4]  + Be[17]*CM[2]* Be[5] )*Be[24];
Matrix[133] = Matrix[133] + ( Be[15]*CM[1]* Be[2]  + Be[17]*CM[2]* Be[0] )*Be[24];
Matrix[134] = Matrix[134] + ( Be[15]*CM[1]* Be[1]  + Be[16]*CM[2]* Be[0] )*Be[24];
Matrix[135] = Matrix[135] + ( Be[15]*CM[0]* Be[0]  + Be[16]*CM[2]* Be[1]  + Be[17]*CM[2]* Be[2] )*Be[24];
Matrix[136] = Matrix[136] + ( Be[15]*CM[2]* Be[15]  + Be[16]*CM[0]* Be[16]  + Be[17]*CM[2]* Be[17] )*Be[24];
Matrix[137] = Matrix[137] + ( Be[16]*CM[1]* Be[15]  + Be[15]*CM[2]* Be[16] )*Be[24];
Matrix[138] = Matrix[138] + ( Be[16]*CM[1]* Be[14]  + Be[17]*CM[2]* Be[13] )*Be[24];
Matrix[139] = Matrix[139] + ( Be[15]*CM[2]* Be[12]  + Be[16]*CM[0]* Be[13]  + Be[17]*CM[2]* Be[14] )*Be[24];
Matrix[140] = Matrix[140] + ( Be[16]*CM[1]* Be[12]  + Be[15]*CM[2]* Be[13] )*Be[24];
Matrix[141] = Matrix[141] + ( Be[16]*CM[1]* Be[11]  + Be[17]*CM[2]* Be[10] )*Be[24];
Matrix[142] = Matrix[142] + ( Be[15]*CM[2]* Be[9]  + Be[16]*CM[0]* Be[10]  + Be[17]*CM[2]* Be[11] )*Be[24];
Matrix[143] = Matrix[143] + ( Be[16]*CM[1]* Be[9]  + Be[15]*CM[2]* Be[10] )*Be[24];
Matrix[144] = Matrix[144] + ( Be[16]*CM[1]* Be[8]  + Be[17]*CM[2]* Be[7] )*Be[24];
Matrix[145] = Matrix[145] + ( Be[15]*CM[2]* Be[6]  + Be[16]*CM[0]* Be[7]  + Be[17]*CM[2]* Be[8] )*Be[24];
Matrix[146] = Matrix[146] + ( Be[16]*CM[1]* Be[6]  + Be[15]*CM[2]* Be[7] )*Be[24];
Matrix[147] = Matrix[147] + ( Be[16]*CM[1]* Be[5]  + Be[17]*CM[2]* Be[4] )*Be[24];
Matrix[148] = Matrix[148] + ( Be[15]*CM[2]* Be[3]  + Be[16]*CM[0]* Be[4]  + Be[17]*CM[2]* Be[5] )*Be[24];
Matrix[149] = Matrix[149] + ( Be[16]*CM[1]* Be[3]  + Be[15]*CM[2]* Be[4] )*Be[24];
Matrix[150] = Matrix[150] + ( Be[16]*CM[1]* Be[2]  + Be[17]*CM[2]* Be[1] )*Be[24];
Matrix[151] = Matrix[151] + ( Be[15]*CM[2]* Be[0]  + Be[16]*CM[0]* Be[1]  + Be[17]*CM[2]* Be[2] )*Be[24];
Matrix[152] = Matrix[152] + ( Be[16]*CM[1]* Be[0]  + Be[15]*CM[2]* Be[1] )*Be[24];
Matrix[153] = Matrix[153] + ( Be[15]*CM[2]* Be[15]  + Be[16]*CM[2]* Be[16]  + Be[17]*CM[0]* Be[17] )*Be[24];
Matrix[154] = Matrix[154] + ( Be[17]*CM[1]* Be[16]  + Be[16]*CM[2]* Be[17] )*Be[24];
Matrix[155] = Matrix[155] + ( Be[17]*CM[1]* Be[15]  + Be[15]*CM[2]* Be[17] )*Be[24];
Matrix[156] = Matrix[156] + ( Be[15]*CM[2]* Be[12]  + Be[16]*CM[2]* Be[13]  + Be[17]*CM[0]* Be[14] )*Be[24];
Matrix[157] = Matrix[157] + ( Be[17]*CM[1]* Be[13]  + Be[16]*CM[2]* Be[14] )*Be[24];
Matrix[158] = Matrix[158] + ( Be[17]*CM[1]* Be[12]  + Be[15]*CM[2]* Be[14] )*Be[24];
Matrix[159] = Matrix[159] + ( Be[15]*CM[2]* Be[9]  + Be[16]*CM[2]* Be[10]  + Be[17]*CM[0]* Be[11] )*Be[24];
Matrix[160] = Matrix[160] + ( Be[17]*CM[1]* Be[10]  + Be[16]*CM[2]* Be[11] )*Be[24];
Matrix[161] = Matrix[161] + ( Be[17]*CM[1]* Be[9]  + Be[15]*CM[2]* Be[11] )*Be[24];
Matrix[162] = Matrix[162] + ( Be[15]*CM[2]* Be[6]  + Be[16]*CM[2]* Be[7]  + Be[17]*CM[0]* Be[8] )*Be[24];
Matrix[163] = Matrix[163] + ( Be[17]*CM[1]* Be[7]  + Be[16]*CM[2]* Be[8] )*Be[24];
Matrix[164] = Matrix[164] + ( Be[17]*CM[1]* Be[6]  + Be[15]*CM[2]* Be[8] )*Be[24];
Matrix[165] = Matrix[165] + ( Be[15]*CM[2]* Be[3]  + Be[16]*CM[2]* Be[4]  + Be[17]*CM[0]* Be[5] )*Be[24];
Matrix[166] = Matrix[166] + ( Be[17]*CM[1]* Be[4]  + Be[16]*CM[2]* Be[5] )*Be[24];
Matrix[167] = Matrix[167] + ( Be[17]*CM[1]* Be[3]  + Be[15]*CM[2]* Be[5] )*Be[24];
Matrix[168] = Matrix[168] + ( Be[15]*CM[2]* Be[0]  + Be[16]*CM[2]* Be[1]  + Be[17]*CM[0]* Be[2] )*Be[24];
Matrix[169] = Matrix[169] + ( Be[17]*CM[1]* Be[1]  + Be[16]*CM[2]* Be[2] )*Be[24];
Matrix[170] = Matrix[170] + ( Be[17]*CM[1]* Be[0]  + Be[15]*CM[2]* Be[2] )*Be[24];
Matrix[171] = Matrix[171] + ( Be[18]*CM[0]* Be[18]  + Be[19]*CM[2]* Be[19]  + Be[20]*CM[2]* Be[20] )*Be[24];
Matrix[172] = Matrix[172] + ( Be[18]*CM[1]* Be[17]  + Be[20]*CM[2]* Be[15] )*Be[24];
Matrix[173] = Matrix[173] + ( Be[18]*CM[1]* Be[16]  + Be[19]*CM[2]* Be[15] )*Be[24];
Matrix[174] = Matrix[174] + ( Be[18]*CM[0]* Be[15]  + Be[19]*CM[2]* Be[16]  + Be[20]*CM[2]* Be[17] )*Be[24];
Matrix[175] = Matrix[175] + ( Be[18]*CM[1]* Be[14]  + Be[20]*CM[2]* Be[12] )*Be[24];
Matrix[176] = Matrix[176] + ( Be[18]*CM[1]* Be[13]  + Be[19]*CM[2]* Be[12] )*Be[24];
Matrix[177] = Matrix[177] + ( Be[18]*CM[0]* Be[12]  + Be[19]*CM[2]* Be[13]  + Be[20]*CM[2]* Be[14] )*Be[24];
Matrix[178] = Matrix[178] + ( Be[18]*CM[1]* Be[11]  + Be[20]*CM[2]* Be[9] )*Be[24];
Matrix[179] = Matrix[179] + ( Be[18]*CM[1]* Be[10]  + Be[19]*CM[2]* Be[9] )*Be[24];
Matrix[180] = Matrix[180] + ( Be[18]*CM[0]* Be[9]  + Be[19]*CM[2]* Be[10]  + Be[20]*CM[2]* Be[11] )*Be[24];
Matrix[181] = Matrix[181] + ( Be[18]*CM[1]* Be[8]  + Be[20]*CM[2]* Be[6] )*Be[24];
Matrix[182] = Matrix[182] + ( Be[18]*CM[1]* Be[7]  + Be[19]*CM[2]* Be[6] )*Be[24];
Matrix[183] = Matrix[183] + ( Be[18]*CM[0]* Be[6]  + Be[19]*CM[2]* Be[7]  + Be[20]*CM[2]* Be[8] )*Be[24];
Matrix[184] = Matrix[184] + ( Be[18]*CM[1]* Be[5]  + Be[20]*CM[2]* Be[3] )*Be[24];
Matrix[185] = Matrix[185] + ( Be[18]*CM[1]* Be[4]  + Be[19]*CM[2]* Be[3] )*Be[24];
Matrix[186] = Matrix[186] + ( Be[18]*CM[0]* Be[3]  + Be[19]*CM[2]* Be[4]  + Be[20]*CM[2]* Be[5] )*Be[24];
Matrix[187] = Matrix[187] + ( Be[18]*CM[1]* Be[2]  + Be[20]*CM[2]* Be[0] )*Be[24];
Matrix[188] = Matrix[188] + ( Be[18]*CM[1]* Be[1]  + Be[19]*CM[2]* Be[0] )*Be[24];
Matrix[189] = Matrix[189] + ( Be[18]*CM[0]* Be[0]  + Be[19]*CM[2]* Be[1]  + Be[20]*CM[2]* Be[2] )*Be[24];
Matrix[190] = Matrix[190] + ( Be[18]*CM[2]* Be[18]  + Be[19]*CM[0]* Be[19]  + Be[20]*CM[2]* Be[20] )*Be[24];
Matrix[191] = Matrix[191] + ( Be[19]*CM[1]* Be[18]  + Be[18]*CM[2]* Be[19] )*Be[24];
Matrix[192] = Matrix[192] + ( Be[19]*CM[1]* Be[17]  + Be[20]*CM[2]* Be[16] )*Be[24];
Matrix[193] = Matrix[193] + ( Be[18]*CM[2]* Be[15]  + Be[19]*CM[0]* Be[16]  + Be[20]*CM[2]* Be[17] )*Be[24];
Matrix[194] = Matrix[194] + ( Be[19]*CM[1]* Be[15]  + Be[18]*CM[2]* Be[16] )*Be[24];
Matrix[195] = Matrix[195] + ( Be[19]*CM[1]* Be[14]  + Be[20]*CM[2]* Be[13] )*Be[24];
Matrix[196] = Matrix[196] + ( Be[18]*CM[2]* Be[12]  + Be[19]*CM[0]* Be[13]  + Be[20]*CM[2]* Be[14] )*Be[24];
Matrix[197] = Matrix[197] + ( Be[19]*CM[1]* Be[12]  + Be[18]*CM[2]* Be[13] )*Be[24];
Matrix[198] = Matrix[198] + ( Be[19]*CM[1]* Be[11]  + Be[20]*CM[2]* Be[10] )*Be[24];
Matrix[199] = Matrix[199] + ( Be[18]*CM[2]* Be[9]  + Be[19]*CM[0]* Be[10]  + Be[20]*CM[2]* Be[11] )*Be[24];
Matrix[200] = Matrix[200] + ( Be[19]*CM[1]* Be[9]  + Be[18]*CM[2]* Be[10] )*Be[24];
Matrix[201] = Matrix[201] + ( Be[19]*CM[1]* Be[8]  + Be[20]*CM[2]* Be[7] )*Be[24];
Matrix[202] = Matrix[202] + ( Be[18]*CM[2]* Be[6]  + Be[19]*CM[0]* Be[7]  + Be[20]*CM[2]* Be[8] )*Be[24];
Matrix[203] = Matrix[203] + ( Be[19]*CM[1]* Be[6]  + Be[18]*CM[2]* Be[7] )*Be[24];
Matrix[204] = Matrix[204] + ( Be[19]*CM[1]* Be[5]  + Be[20]*CM[2]* Be[4] )*Be[24];
Matrix[205] = Matrix[205] + ( Be[18]*CM[2]* Be[3]  + Be[19]*CM[0]* Be[4]  + Be[20]*CM[2]* Be[5] )*Be[24];
Matrix[206] = Matrix[206] + ( Be[19]*CM[1]* Be[3]  + Be[18]*CM[2]* Be[4] )*Be[24];
Matrix[207] = Matrix[207] + ( Be[19]*CM[1]* Be[2]  + Be[20]*CM[2]* Be[1] )*Be[24];
Matrix[208] = Matrix[208] + ( Be[18]*CM[2]* Be[0]  + Be[19]*CM[0]* Be[1]  + Be[20]*CM[2]* Be[2] )*Be[24];
Matrix[209] = Matrix[209] + ( Be[19]*CM[1]* Be[0]  + Be[18]*CM[2]* Be[1] )*Be[24];
Matrix[210] = Matrix[210] + ( Be[18]*CM[2]* Be[18]  + Be[19]*CM[2]* Be[19]  + Be[20]*CM[0]* Be[20] )*Be[24];
Matrix[211] = Matrix[211] + ( Be[20]*CM[1]* Be[19]  + Be[19]*CM[2]* Be[20] )*Be[24];
Matrix[212] = Matrix[212] + ( Be[20]*CM[1]* Be[18]  + Be[18]*CM[2]* Be[20] )*Be[24];
Matrix[213] = Matrix[213] + ( Be[18]*CM[2]* Be[15]  + Be[19]*CM[2]* Be[16]  + Be[20]*CM[0]* Be[17] )*Be[24];
Matrix[214] = Matrix[214] + ( Be[20]*CM[1]* Be[16]  + Be[19]*CM[2]* Be[17] )*Be[24];
Matrix[215] = Matrix[215] + ( Be[20]*CM[1]* Be[15]  + Be[18]*CM[2]* Be[17] )*Be[24];
Matrix[216] = Matrix[216] + ( Be[18]*CM[2]* Be[12]  + Be[19]*CM[2]* Be[13]  + Be[20]*CM[0]* Be[14] )*Be[24];
Matrix[217] = Matrix[217] + ( Be[20]*CM[1]* Be[13]  + Be[19]*CM[2]* Be[14] )*Be[24];
Matrix[218] = Matrix[218] + ( Be[20]*CM[1]* Be[12]  + Be[18]*CM[2]* Be[14] )*Be[24];
Matrix[219] = Matrix[219] + ( Be[18]*CM[2]* Be[9]  + Be[19]*CM[2]* Be[10]  + Be[20]*CM[0]* Be[11] )*Be[24];
Matrix[220] = Matrix[220] + ( Be[20]*CM[1]* Be[10]  + Be[19]*CM[2]* Be[11] )*Be[24];
Matrix[221] = Matrix[221] + ( Be[20]*CM[1]* Be[9]  + Be[18]*CM[2]* Be[11] )*Be[24];
Matrix[222] = Matrix[222] + ( Be[18]*CM[2]* Be[6]  + Be[19]*CM[2]* Be[7]  + Be[20]*CM[0]* Be[8] )*Be[24];
Matrix[223] = Matrix[223] + ( Be[20]*CM[1]* Be[7]  + Be[19]*CM[2]* Be[8] )*Be[24];
Matrix[224] = Matrix[224] + ( Be[20]*CM[1]* Be[6]  + Be[18]*CM[2]* Be[8] )*Be[24];
Matrix[225] = Matrix[225] + ( Be[18]*CM[2]* Be[3]  + Be[19]*CM[2]* Be[4]  + Be[20]*CM[0]* Be[5] )*Be[24];
Matrix[226] = Matrix[226] + ( Be[20]*CM[1]* Be[4]  + Be[19]*CM[2]* Be[5] )*Be[24];
Matrix[227] = Matrix[227] + ( Be[20]*CM[1]* Be[3]  + Be[18]*CM[2]* Be[5] )*Be[24];
Matrix[228] = Matrix[228] + ( Be[18]*CM[2]* Be[0]  + Be[19]*CM[2]* Be[1]  + Be[20]*CM[0]* Be[2] )*Be[24];
Matrix[229] = Matrix[229] + ( Be[20]*CM[1]* Be[1]  + Be[19]*CM[2]* Be[2] )*Be[24];
Matrix[230] = Matrix[230] + ( Be[20]*CM[1]* Be[0]  + Be[18]*CM[2]* Be[2] )*Be[24];
Matrix[231] = Matrix[231] + ( Be[21]*CM[0]* Be[21]  + Be[22]*CM[2]* Be[22]  + Be[23]*CM[2]* Be[23] )*Be[24];
Matrix[232] = Matrix[232] + ( Be[21]*CM[1]* Be[20]  + Be[23]*CM[2]* Be[18] )*Be[24];
Matrix[233] = Matrix[233] + ( Be[21]*CM[1]* Be[19]  + Be[22]*CM[2]* Be[18] )*Be[24];
Matrix[234] = Matrix[234] + ( Be[21]*CM[0]* Be[18]  + Be[22]*CM[2]* Be[19]  + Be[23]*CM[2]* Be[20] )*Be[24];
Matrix[235] = Matrix[235] + ( Be[21]*CM[1]* Be[17]  + Be[23]*CM[2]* Be[15] )*Be[24];
Matrix[236] = Matrix[236] + ( Be[21]*CM[1]* Be[16]  + Be[22]*CM[2]* Be[15] )*Be[24];
Matrix[237] = Matrix[237] + ( Be[21]*CM[0]* Be[15]  + Be[22]*CM[2]* Be[16]  + Be[23]*CM[2]* Be[17] )*Be[24];
Matrix[238] = Matrix[238] + ( Be[21]*CM[1]* Be[14]  + Be[23]*CM[2]* Be[12] )*Be[24];
Matrix[239] = Matrix[239] + ( Be[21]*CM[1]* Be[13]  + Be[22]*CM[2]* Be[12] )*Be[24];
Matrix[240] = Matrix[240] + ( Be[21]*CM[0]* Be[12]  + Be[22]*CM[2]* Be[13]  + Be[23]*CM[2]* Be[14] )*Be[24];
Matrix[241] = Matrix[241] + ( Be[21]*CM[1]* Be[11]  + Be[23]*CM[2]* Be[9] )*Be[24];
Matrix[242] = Matrix[242] + ( Be[21]*CM[1]* Be[10]  + Be[22]*CM[2]* Be[9] )*Be[24];
Matrix[243] = Matrix[243] + ( Be[21]*CM[0]* Be[9]  + Be[22]*CM[2]* Be[10]  + Be[23]*CM[2]* Be[11] )*Be[24];
Matrix[244] = Matrix[244] + ( Be[21]*CM[1]* Be[8]  + Be[23]*CM[2]* Be[6] )*Be[24];
Matrix[245] = Matrix[245] + ( Be[21]*CM[1]* Be[7]  + Be[22]*CM[2]* Be[6] )*Be[24];
Matrix[246] = Matrix[246] + ( Be[21]*CM[0]* Be[6]  + Be[22]*CM[2]* Be[7]  + Be[23]*CM[2]* Be[8] )*Be[24];
Matrix[247] = Matrix[247] + ( Be[21]*CM[1]* Be[5]  + Be[23]*CM[2]* Be[3] )*Be[24];
Matrix[248] = Matrix[248] + ( Be[21]*CM[1]* Be[4]  + Be[22]*CM[2]* Be[3] )*Be[24];
Matrix[249] = Matrix[249] + ( Be[21]*CM[0]* Be[3]  + Be[22]*CM[2]* Be[4]  + Be[23]*CM[2]* Be[5] )*Be[24];
Matrix[250] = Matrix[250] + ( Be[21]*CM[1]* Be[2]  + Be[23]*CM[2]* Be[0] )*Be[24];
Matrix[251] = Matrix[251] + ( Be[21]*CM[1]* Be[1]  + Be[22]*CM[2]* Be[0] )*Be[24];
Matrix[252] = Matrix[252] + ( Be[21]*CM[0]* Be[0]  + Be[22]*CM[2]* Be[1]  + Be[23]*CM[2]* Be[2] )*Be[24];
Matrix[253] = Matrix[253] + ( Be[21]*CM[2]* Be[21]  + Be[22]*CM[0]* Be[22]  + Be[23]*CM[2]* Be[23] )*Be[24];
Matrix[254] = Matrix[254] + ( Be[22]*CM[1]* Be[21]  + Be[21]*CM[2]* Be[22] )*Be[24];
Matrix[255] = Matrix[255] + ( Be[22]*CM[1]* Be[20]  + Be[23]*CM[2]* Be[19] )*Be[24];
Matrix[256] = Matrix[256] + ( Be[21]*CM[2]* Be[18]  + Be[22]*CM[0]* Be[19]  + Be[23]*CM[2]* Be[20] )*Be[24];
Matrix[257] = Matrix[257] + ( Be[22]*CM[1]* Be[18]  + Be[21]*CM[2]* Be[19] )*Be[24];
Matrix[258] = Matrix[258] + ( Be[22]*CM[1]* Be[17]  + Be[23]*CM[2]* Be[16] )*Be[24];
Matrix[259] = Matrix[259] + ( Be[21]*CM[2]* Be[15]  + Be[22]*CM[0]* Be[16]  + Be[23]*CM[2]* Be[17] )*Be[24];
Matrix[260] = Matrix[260] + ( Be[22]*CM[1]* Be[15]  + Be[21]*CM[2]* Be[16] )*Be[24];
Matrix[261] = Matrix[261] + ( Be[22]*CM[1]* Be[14]  + Be[23]*CM[2]* Be[13] )*Be[24];
Matrix[262] = Matrix[262] + ( Be[21]*CM[2]* Be[12]  + Be[22]*CM[0]* Be[13]  + Be[23]*CM[2]* Be[14] )*Be[24];
Matrix[263] = Matrix[263] + ( Be[22]*CM[1]* Be[12]  + Be[21]*CM[2]* Be[13] )*Be[24];
Matrix[264] = Matrix[264] + ( Be[22]*CM[1]* Be[11]  + Be[23]*CM[2]* Be[10] )*Be[24];
Matrix[265] = Matrix[265] + ( Be[21]*CM[2]* Be[9]  + Be[22]*CM[0]* Be[10]  + Be[23]*CM[2]* Be[11] )*Be[24];
Matrix[266] = Matrix[266] + ( Be[22]*CM[1]* Be[9]  + Be[21]*CM[2]* Be[10] )*Be[24];
Matrix[267] = Matrix[267] + ( Be[22]*CM[1]* Be[8]  + Be[23]*CM[2]* Be[7] )*Be[24];
Matrix[268] = Matrix[268] + ( Be[21]*CM[2]* Be[6]  + Be[22]*CM[0]* Be[7]  + Be[23]*CM[2]* Be[8] )*Be[24];
Matrix[269] = Matrix[269] + ( Be[22]*CM[1]* Be[6]  + Be[21]*CM[2]* Be[7] )*Be[24];
Matrix[270] = Matrix[270] + ( Be[22]*CM[1]* Be[5]  + Be[23]*CM[2]* Be[4] )*Be[24];
Matrix[271] = Matrix[271] + ( Be[21]*CM[2]* Be[3]  + Be[22]*CM[0]* Be[4]  + Be[23]*CM[2]* Be[5] )*Be[24];
Matrix[272] = Matrix[272] + ( Be[22]*CM[1]* Be[3]  + Be[21]*CM[2]* Be[4] )*Be[24];
Matrix[273] = Matrix[273] + ( Be[22]*CM[1]* Be[2]  + Be[23]*CM[2]* Be[1] )*Be[24];
Matrix[274] = Matrix[274] + ( Be[21]*CM[2]* Be[0]  + Be[22]*CM[0]* Be[1]  + Be[23]*CM[2]* Be[2] )*Be[24];
Matrix[275] = Matrix[275] + ( Be[22]*CM[1]* Be[0]  + Be[21]*CM[2]* Be[1] )*Be[24];
Matrix[276] = Matrix[276] + ( Be[21]*CM[2]* Be[21]  + Be[22]*CM[2]* Be[22]  + Be[23]*CM[0]* Be[23] )*Be[24];
Matrix[277] = Matrix[277] + ( Be[23]*CM[1]* Be[22]  + Be[22]*CM[2]* Be[23] )*Be[24];
Matrix[278] = Matrix[278] + ( Be[23]*CM[1]* Be[21]  + Be[21]*CM[2]* Be[23] )*Be[24];
Matrix[279] = Matrix[279] + ( Be[21]*CM[2]* Be[18]  + Be[22]*CM[2]* Be[19]  + Be[23]*CM[0]* Be[20] )*Be[24];
Matrix[280] = Matrix[280] + ( Be[23]*CM[1]* Be[19]  + Be[22]*CM[2]* Be[20] )*Be[24];
Matrix[281] = Matrix[281] + ( Be[23]*CM[1]* Be[18]  + Be[21]*CM[2]* Be[20] )*Be[24];
Matrix[282] = Matrix[282] + ( Be[21]*CM[2]* Be[15]  + Be[22]*CM[2]* Be[16]  + Be[23]*CM[0]* Be[17] )*Be[24];
Matrix[283] = Matrix[283] + ( Be[23]*CM[1]* Be[16]  + Be[22]*CM[2]* Be[17] )*Be[24];
Matrix[284] = Matrix[284] + ( Be[23]*CM[1]* Be[15]  + Be[21]*CM[2]* Be[17] )*Be[24];
Matrix[285] = Matrix[285] + ( Be[21]*CM[2]* Be[12]  + Be[22]*CM[2]* Be[13]  + Be[23]*CM[0]* Be[14] )*Be[24];
Matrix[286] = Matrix[286] + ( Be[23]*CM[1]* Be[13]  + Be[22]*CM[2]* Be[14] )*Be[24];
Matrix[287] = Matrix[287] + ( Be[23]*CM[1]* Be[12]  + Be[21]*CM[2]* Be[14] )*Be[24];
Matrix[288] = Matrix[288] + ( Be[21]*CM[2]* Be[9]  + Be[22]*CM[2]* Be[10]  + Be[23]*CM[0]* Be[11] )*Be[24];
Matrix[289] = Matrix[289] + ( Be[23]*CM[1]* Be[10]  + Be[22]*CM[2]* Be[11] )*Be[24];
Matrix[290] = Matrix[290] + ( Be[23]*CM[1]* Be[9]  + Be[21]*CM[2]* Be[11] )*Be[24];
Matrix[291] = Matrix[291] + ( Be[21]*CM[2]* Be[6]  + Be[22]*CM[2]* Be[7]  + Be[23]*CM[0]* Be[8] )*Be[24];
Matrix[292] = Matrix[292] + ( Be[23]*CM[1]* Be[7]  + Be[22]*CM[2]* Be[8] )*Be[24];
Matrix[293] = Matrix[293] + ( Be[23]*CM[1]* Be[6]  + Be[21]*CM[2]* Be[8] )*Be[24];
Matrix[294] = Matrix[294] + ( Be[21]*CM[2]* Be[3]  + Be[22]*CM[2]* Be[4]  + Be[23]*CM[0]* Be[5] )*Be[24];
Matrix[295] = Matrix[295] + ( Be[23]*CM[1]* Be[4]  + Be[22]*CM[2]* Be[5] )*Be[24];
Matrix[296] = Matrix[296] + ( Be[23]*CM[1]* Be[3]  + Be[21]*CM[2]* Be[5] )*Be[24];
Matrix[297] = Matrix[297] + ( Be[21]*CM[2]* Be[0]  + Be[22]*CM[2]* Be[1]  + Be[23]*CM[0]* Be[2] )*Be[24];
Matrix[298] = Matrix[298] + ( Be[23]*CM[1]* Be[1]  + Be[22]*CM[2]* Be[2] )*Be[24];
Matrix[299] = Matrix[299] + ( Be[23]*CM[1]* Be[0]  + Be[21]*CM[2]* Be[2] )*Be[24];

}