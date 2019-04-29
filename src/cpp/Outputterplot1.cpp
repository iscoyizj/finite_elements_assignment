/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Outputterplot1.h"
#include "SkylineMatrix.h"

#include <iostream>
#include <iomanip>
#include <ctime>

using namespace std;

COutputterplot1* COutputterplot1::_instanceplot1 = nullptr;

//	Constructor
COutputterplot1::COutputterplot1(string FileName)
{
	OutputFileplot1.open(FileName);

	if (!OutputFileplot1)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputterplot1* COutputterplot1::Instanceplot1(string FileName)
{
	if (!_instanceplot1)
		_instanceplot1 = new COutputterplot1(FileName);
	return _instanceplot1;
}

//	Print program logo
void COutputterplot1::OutputHeading()
{
	CDomain* FEMData = CDomain::Instance();

	*this << "TITLE = " << "\""<< "The Stress Nephogram of" << FEMData->GetTitle() <<"\""<<endl;
	
}

//	Calculate stresses
void COutputterplot1::OutputElementStress()
{
	CDomain* FEMData = CDomain::Instance();

	double* Displacement = FEMData->GetDisplacement();

	unsigned int NUMEG = FEMData->GetNUMEG();

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{

		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();
	
		switch (ElementType)
		{
			case ElementTypes::Q4: // 4Q element
			*this << "  ELEMENT  LOCATION    X       Y       Z        SXX            SYY            TXY" << endl
				<< "  NUMBER   POINT" << endl;

			for (unsigned int Ele = 0; Ele < NUME; Ele++)
			{
				double *newlocation = new double[36];
				for (unsigned int m = 0; m < 36; m++)
				{
					newlocation[m] = 0;
				}
				CElement& Element = EleGrp[Ele];
				Element.ElementStressplot1(newlocation, Displacement);

				for (unsigned int i = 0; i < 4; i++)
				{
					*this << setw(5) << Ele + 1 << setw(9) << i+1
						<< setw(15) << newlocation[0 + 9 * i] << setw(15) << newlocation[1 + 9 * i] << setw(15) << newlocation[2 + 9 * i]
						<< setw(15) << newlocation[3 + 9 * i] << setw(15) << newlocation[4 + 9 * i] << setw(15) << newlocation[6 + 9 * i] << endl;
					
				}
				delete[] newlocation;
			}
				

				*this << endl;
				break;
		}
	}
}