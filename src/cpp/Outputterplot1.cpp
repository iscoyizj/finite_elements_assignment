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

#define coeff 1
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
				double *newlocation = new double[12];
				clear(newlocation, 12);
				double *stress_4Q = new double[24];
				double *origin_loc = new double[12];
				clear(origin_loc, 12);
				clear(stress_4Q, 24);
				CElement& Element = EleGrp[Ele];
				Element.ElementPostInfo(stress_4Q,Displacement,origin_loc,newlocation);

				for (unsigned int i = 0; i < 4; i++)
				{
					*this << setw(5) << Ele + 1 << setw(9) << i + 1
						<< setw(15) << newlocation[0 + 3 * i] << setw(15) << newlocation[1 + 3 * i] << setw(15) << newlocation[2 + 3 * i]
						<< setw(15) << stress_4Q[0 + 6 * i] << setw(15) << stress_4Q[1 + 6 * i] << setw(15) << stress_4Q[3 + 6 * i] << endl;
					
				}
				delete[] newlocation;
			}
			*this << endl;
			break;

			case ElementTypes::Beam: // Beam element
				*this << "  ELEMENT  LOCATION    X       Y       Z        SXX            SYY            SZZ            TXY            TYZ            TXZ" << endl
					<< "  NUMBER   POINT" << endl;

				double beamstress[48];
				double prePositionBeam[24];
				double postPositionBeam[24];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementPostInfo(beamstress, Displacement, prePositionBeam,postPositionBeam);

					CBeamMaterial& material =
						*dynamic_cast<CBeamMaterial*>(Element.GetElementMaterial());
				
					for (unsigned i = 0; i < 8; i++)
					{
						*this << setw(5) << Ele + 1 << setw(12) << i + 1;
						for (unsigned DegOF = 0; DegOF < 3; DegOF++)
						{
							*this << setw(15)
								<< (1 - coeff) * prePositionBeam[3 * i + DegOF] +
								coeff * postPositionBeam[3 * i + DegOF];
						}
						for (unsigned DegOF = 0; DegOF < 6; DegOF++)
						{
							*this << setw(15) << beamstress[6 * i + DegOF];
						}
						*this << endl;
					}
				}

				*this << endl;

				break;
				

				
		}
	}
}