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
#include "OutPlotPost.h"
#include "SkylineMatrix.h"

#include <iostream>
#include <iomanip>

using namespace std;


COutPlotPost* COutPlotPost::_instance = nullptr;

//	Constructor
COutPlotPost::COutPlotPost(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutPlotPost* COutPlotPost::Instance(string FileName)
{
	if (!_instance)
		_instance = new COutPlotPost(FileName);
	return _instance;
}

void COutPlotPost::OutputHeading()
{
	CDomain* FEMData = CDomain::Instance();
	*this << "# vtk DataFile Version 3.0" << endl;
	*this << FEMData->GetTitle() << endl;
	*this << "ASCII" << endl;
	*this << "DATASET UNSTRUCTURED_GRID" << endl;
	*this << endl;
}

//	Print nodal data
void COutPlotPost::OutNode()
{
	CDomain* FEMData = CDomain::Instance();

	CNode* NodeList = FEMData->GetNodeList();

	unsigned int NUMNP = FEMData->GetNUMNP();

	double* Displacement = FEMData->GetDisplacement();

	*this << "POINTS " << setw(9) << NUMNP << setw(9) << "double" << endl;
	
	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		{
			for (unsigned int i=0; i<3; i++)
			{
				if (NodeList[np].bcode[i] == 0)
				{
					*this << NodeList[np].XYZ[i] <<setw(20) ;
				}
				else
				{
					*this << NodeList[np].XYZ[i] + Displacement[NodeList[np].bcode[i] - 1]*MAGN <<setw(20) ;
				}
			}
			*this << endl;
		}

	*this << endl;
}


//	Output element data
void COutPlotPost::OutputElementInfo(unsigned int n, unsigned int sum)
{
	//	Print element group control line
	*this << "CELLS" << setw(12) << sum << setw(12) << n << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
		unsigned int NUME = ElementGroup.GetNUME();

		//	Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
			ElementGroup[Ele].WritePlotPost(*this, Ele);
	}
	*this << endl;
}

void COutPlotPost::OutputEleType(double n, double nnd)
{
	CDomain* FEMData = CDomain::Instance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << "CELL_TYPES"  << setw(12) << n << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
		unsigned int NUME = ElementGroup.GetNUME();
		ElementTypes ElementType = ElementGroup.GetElementType();
		unsigned int type_;
		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				type_=3;
				break;
			case ElementTypes::Q4: // 4Q element
				type_=7;
				break;
			case ElementTypes::T3: // 4Q element
				type_=5;
				break;
			case ElementTypes::H8: // 4Q element
				type_=11;
				break;
			case ElementTypes::Plate: // Plate element
				type_=7;
				break;
			case ElementTypes::Beam: // Beam element
				type_=3;
				break;
			case ElementTypes::Shell: // Shell element
				type_ = 7;
				break;
		}
		for (unsigned int i=0;i<NUME;i++)
			*this << type_ << endl;
	}
	*this << endl;
	*this << "POINT_DATA" << setw(12) << nnd << endl;
}

void COutPlotPost::OutputNodalDisplacement(unsigned int lcase)
{
	CDomain* FEMData = CDomain::Instance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();
	string label[6] = { "x", "y", "z", "theta_1", "theta_2", "theta_3"};

	for (unsigned int i=0;i<6;i++)
	{
		*this << "SCALARS" << setw(5) << "load" << lcase+1 << "_" << label[i] << setw(5) << "   double  1" << endl;
		*this << "LOOKUP_TABLE" << setw(5) << "load" << lcase+1 << "_" << label[i] << endl;
		for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		{
			if (NodeList[np].bcode[i] == 0)
			{
				*this << 0.0 << endl;
			}
			else
			{
				*this << Displacement[NodeList[np].bcode[i] - 1] << endl;
			}
		}
	}

	*this << endl;
}

void COutPlotPost::StressHead(unsigned int lcase, unsigned int n)
{
	*this << "CELL_DATA" << setw(12) << n << endl;
	*this << "SCALARS" << setw(5) << "load" << lcase+1 << "_Mises" << setw(5) << "   double  1" << endl;
	*this << "LOOKUP_TABLE" << setw(5) << "load" << lcase+1 << "_" << "_Mises" << endl;
}

void COutPlotPost::ElementStress(double stress)
{
	*this << stress <<endl;
}





