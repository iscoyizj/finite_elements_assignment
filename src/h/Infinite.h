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

//! Infinite element class
class CInfi : public CElement
{
public:

//!	Constructor
	CInfi();

//!	Desconstructor
	~CInfi();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output, unsigned int Ele);

	virtual void WritePlot(COutPlot& output, unsigned int Ele);

	virtual void WritePlotPost(COutPlotPost& output, unsigned int Ele);

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stiffness matrix 
	virtual void ElementMass(double* Mass){Mass = 0;};

//!	Calculate element stress
	virtual void ElementStress( double* stress, double* Displacement);

//!	Calculate the values required in the POSTPROCESS 
	virtual void ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* PostPositions);

//! Recover element stress
	virtual void RecoverElementStress(double* Displacement, double* A);
	
//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();
// Caculate Gravity of Elements
	virtual void GravityCalculation(double* ptr_force);

};
