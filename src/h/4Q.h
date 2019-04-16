#pragma once
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
class C4Q : public CElement
{
public:

	//!	Constructor
	C4Q();

	//!	Desconstructor
	~C4Q();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output, unsigned int Ele);

	//! Generate location matrix: the global equation number that corresponding to each DOF of the element
	//	Caution:  Equation number is numbered from 1 !
	virtual void GenerateLocationMatrix();

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

	//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();
	
	virtual double GravityofElement();

	//£¡ Return the N matrix for 4Q element
	void Nmat4Q(double eta, double psi, double * Nmat);
	
	//!  Return the B matrix for 4Q element
	double Bmat4Q(double eta, double psi, double * Bmat);

	//!  Gauss Integral at gauss node,consists of element stiffness matrix
	void Integral_gauss_node(const double & eta, const double & psi, const double & weight, double * k_e);

	void ElementStress_gauss_node(double eta, double psi, double * stress, double * Displacement);

	void Gauss_node_coordinate(double eta, double psi, double * coordinte);

	void Stress_coordinate(double * coordinate);
    
};
