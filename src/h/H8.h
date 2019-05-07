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
class CH8 : public CElement
{
public:

//!	Constructor
	CH8();

//!	Desconstructor
	~CH8();

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

//	virtual void ElementGravity (double* EG) ;

	//!	Calculate element stiffness matrix 
	virtual void ElementMass(double* Mass){;};

	//! Recover element stress
	virtual void RecoverElementStress(double* Displacement, double* A){;};

	// Caculate Gravity of Elements
	virtual void GravityCalculation(double* ptr_force);

	//!	Calculate the values required in the POSTPROCESS 
	virtual void ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* PostPositions);


	
	void ElementCoord (double* coord);

	void CalculateBe (double* Be, double psi, double eta, double zet);
	
	void CalculateMatrix (double* Be, double* Matrix);
	
	void CalculateStress (double* Be, double* stress, double* Displacement);
	
	void CalculateCoord (double* coord, double psi, double eta, double zet);
	
	double GetJaccobi (double psi, double eta, double zet);
	
	void CalculateVolume (double* EG, double psi, double eta, double zet);
};
