/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! Bar element class
class CSubpara : public CElement
{
public:

//!	Constructor
	CSubpara();

//!	Desconstructor
	~CSubpara();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output, unsigned int Ele);

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

// Caculate Gravity of Elements
	virtual void GravityCalculation(double* ptr_force) = 0;

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();
	
//!	Calculate element stiffness matrix 
	virtual void ElementMass(double* Mass){Mass = 0;};
	
//!	Calculate the values required in the POSTPROCESS 
	virtual void ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* PostPositions);

//! Recover element stress
	virtual void RecoverElementStress(double* Displacement, double* A);

private:

//! calculate shape function N at (xi,eta)
	void calN(double N[], double xi, double eta);

//! calculate [GB] at (xi,eta),dimension of GB is 2*NEN_
	void calGB(double GB[], double xi, double eta);

//! calculate [B] at (xi,eta),dimension of B is 2*NEN_
	void calB(double B[], double xi, double eta);

//! calculate invJ and determint J at (xi,eta)
	double detinvJ(double invJ[][2], double GB[], double xi, double eta);

//! detJ at (xi,eta)
	double detJ;

//! gaussian points
	double gspos[3]; double w[3];
//! consititutive matrix
	double D[3];
};
