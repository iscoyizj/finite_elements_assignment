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
#include "Material.h"
#include "CSRMatrix.h"
using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	
	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;
	NUMELE = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::Instance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile, string PlotFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::Instance(OutFile);
	COutPlot* Outplot = COutPlot::Instance(PlotFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();
	Outplot->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	if (ReadNodalPoints())
	{
		Output->OutputNodeInfo();
		Outplot->OutNode();
	}
    else
        return false;

//	Read load data
	if (ReadLoadCases())
		Output->OutputLoadInfo();
    else
        return false;

	unsigned int n = 0;
	//	Read element data
	if (ReadElements(&n, &NUMELE))
	{
		Output->OutputElementInfo();
		Outplot->OutputElementInfo(n, NUMELE);
		Outplot->OutputEleType(NUMELE, NUMNP);
	}
	else
		return false;

	//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{
//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
		if (!NodeList[np].Read(Input, np))
			return false;

	return true;
}



//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
		if (!LoadCases[lcase].Read(Input, lcase))
			return false;

	return true;
}


// Read element data
bool CDomain::ReadElements(unsigned int* n, unsigned int* sum)
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
        if (!EleGrpList[EleGrp].Read(Input,n,sum))
            return false;
    
    return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
	{
		CElementGroup& ElementGrp = EleGrpList[EleGrp];
		ElementTypes ElementType_;
		ElementType_ = ElementGrp.GetElementType();
		switch (ElementType_)
		{
			//		 case Beam: // Bar element
			//		{
			//			unsigned int NUME = ElementGrp.GetNUME();
			//			for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
			//			{
			//				CElement& Element = ElementGrp[Ele];
			//				CNode** node_ = Element.CElement::GetNodes();
			//				unsigned int node_left = node_[0]->NodeNumber;
			//				unsigned int node_right = node_[1]->NodeNumber;
			//				NodeList[node_left - 1].bcode[3] = 0;
			//				NodeList[node_left - 1].bcode[4] = 0;
			//				NodeList[node_left - 1].bcode[5] = 0;
			//				NodeList[node_right - 1].bcode[3] = 0;
			//				NodeList[node_right - 1].bcode[4] = 0;
			//			    NodeList[node_right - 1].bcode[5] = 0;
			//			}
			//		}
			//		break;

			//		case Shell: // Shell element
			//		{
			//			unsigned int NUME = ElementGrp.GetNUME();
			//			for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
			//			{
			//				CElement& Element = ElementGrp[Ele];
			//				CNode** node_ = Element.CElement::GetNodes();
			//				unsigned int node_one = node_[0]->NodeNumber;
			//				unsigned int node_two = node_[1]->NodeNumber;
			//				unsigned int node_three = node_[2]->NodeNumber;
			//				unsigned int node_four = node_[3]->NodeNumber;
			//				NodeList[node_one - 1].bcode[3] = 0;
			//				NodeList[node_one - 1].bcode[4] = 0;
			//				NodeList[node_one - 1].bcode[5] = 0;
			//				NodeList[node_two - 1].bcode[3] = 0;
				//			NodeList[node_two - 1].bcode[4] = 0;
			//				NodeList[node_two - 1].bcode[5] = 0;
			//				NodeList[node_three - 1].bcode[3] = 0;
			//				NodeList[node_three - 1].bcode[4] = 0;
			//				NodeList[node_three - 1].bcode[5] = 0;
			//				NodeList[node_four - 1].bcode[3] = 0;
			//				NodeList[node_four - 1].bcode[4] = 0;
			//				NodeList[node_four - 1].bcode[5] = 0;
			//			}
			//		}
			//		break;
		case Bar: // Bar element
		{
			unsigned int NUME = ElementGrp.GetNUME();
			for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
			{
				CElement& Element = ElementGrp[Ele];
				CNode** node_ = Element.CElement::GetNodes();
				unsigned int node_left = node_[0]->NodeNumber;
				unsigned int node_right = node_[1]->NodeNumber;
				NodeList[node_left - 1].bcode[3] = 1;
				NodeList[node_left - 1].bcode[4] = 1;

				NodeList[node_right - 1].bcode[3] = 1;
				NodeList[node_right - 1].bcode[4] = 1;

			}
		}
		break;
		case H8: // 8H element
		{
			unsigned int NUME = ElementGrp.GetNUME();
			for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
			{
				CElement& Element = ElementGrp[Ele];
				CNode** node_ = Element.CElement::GetNodes();
				unsigned int node_1 = node_[0]->NodeNumber;
				unsigned int node_2 = node_[1]->NodeNumber;
				unsigned int node_3 = node_[2]->NodeNumber;
				unsigned int node_4 = node_[3]->NodeNumber;
				unsigned int node_5 = node_[4]->NodeNumber;
				unsigned int node_6 = node_[5]->NodeNumber;
				unsigned int node_7 = node_[6]->NodeNumber;
				unsigned int node_8 = node_[7]->NodeNumber;
				NodeList[node_1 - 1].bcode[3] = 1;
				NodeList[node_1 - 1].bcode[4] = 1;
				NodeList[node_2 - 1].bcode[3] = 1;
				NodeList[node_2 - 1].bcode[4] = 1;
				NodeList[node_3 - 1].bcode[3] = 1;
				NodeList[node_3 - 1].bcode[4] = 1;
				NodeList[node_4 - 1].bcode[3] = 1;
				NodeList[node_4 - 1].bcode[4] = 1;
				NodeList[node_5 - 1].bcode[3] = 1;
				NodeList[node_5 - 1].bcode[4] = 1;
				NodeList[node_6 - 1].bcode[3] = 1;
				NodeList[node_6 - 1].bcode[4] = 1;
				NodeList[node_7 - 1].bcode[3] = 1;
				NodeList[node_7 - 1].bcode[4] = 1;
				NodeList[node_8 - 1].bcode[3] = 1;
				NodeList[node_8 - 1].bcode[4] = 1;
			}
		}
		break;
		}
	}
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof])
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // Generate location matrix
            Element.GenerateLocationMatrix();
            
            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }
    
    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	Output->PrintColumnHeights();
#endif

}

void CDomain::AssembleStiffnessMatrix()
{
	//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		CElementGroup& ElementGrp = EleGrpList[EleGrp];
		unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

		//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = ElementGrp[Ele];
			Element.ElementStiffness(Matrix);

#ifdef _PARDISO_
			CSRStiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
#else
			StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
#ifdef _DEBUG_
			COutputter* Output = COutputter::Instance();
			Output->PrintStiffnessMatrix();
#endif
#endif // _PARDISO_
		}

		delete[] Matrix;
		Matrix = nullptr;
	}
}


//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
	}
	return true;
}
#ifdef _PARDISO_
void CDomain::CalculateCSRColumns()
{
	CSparseMatrix<double>* matrix = CSRStiffnessMatrix;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		CElementGroup& ElementGrp = EleGrpList[EleGrp];
		unsigned int NUME = ElementGrp.GetNUME();

		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = ElementGrp[Ele];
			Element.GenerateLocationMatrix();
			unsigned LMSize = Element.GetND();
			unsigned* LM = Element.GetLocationMatrix();
			for (unsigned i = 0; i < LMSize; ++i)
			{
				unsigned index1 = LM[i];
				if (!index1) continue;
				for (unsigned j = i; j < LMSize; ++j)
				{
					unsigned index2 = LM[j];
					if (!index2) continue;
					if (index1 < index2) {
						matrix->markPosition(index1, index2);
					}
					else
					{
						matrix->markPosition(index2, index1);
					}
				}
			}
		}
	}
}
#endif // _PARDISO_
void CDomain::Gravity()
{

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
	{
		CElementGroup& ElementGrp = EleGrpList[EleGrp];
		ElementTypes ElementType_;
		ElementType_ = ElementGrp.GetElementType();
		double* ptr_force = nullptr;
		switch (ElementType_ )
		{
		 case Bar: //Bar element
		     {
			    unsigned int NUME = ElementGrp.GetNUME();

				for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
				{
					CElement& Element = ElementGrp[Ele];
					Element.GravityCalculation(ptr_force);
					CNode** node_ = Element.CElement::GetNodes();
					unsigned int node_left = node_[0]->NodeNumber;
					unsigned int node_right = node_[1]->NodeNumber;
					unsigned int dof_left = NodeList[node_left - 1].bcode[2];
					unsigned int dof_right = NodeList[node_right - 1].bcode[2];
					if (dof_left)
						Force[dof_left - 1] += -Element.GetGravity() / 2;
					if (dof_right)
						Force[dof_right - 1] += -Element.GetGravity() / 2;
				}
		     }
			 break;
		 case Q4: //4Q element
		 {
			 unsigned int NUME = ElementGrp.GetNUME();
			 double* ptr_force = nullptr;
			 for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
			 {
				 CElement& Element = ElementGrp[Ele];
				 Element.GravityCalculation(ptr_force);
				 CNode** node_ = Element.CElement::GetNodes();
				 unsigned int node_one = node_[0]->NodeNumber;
				 unsigned int node_two = node_[1]->NodeNumber;
				 unsigned int node_three = node_[2]->NodeNumber;
				 unsigned int node_four = node_[3]->NodeNumber;
				 unsigned int dof_one = NodeList[node_one - 1].bcode[2];
				 unsigned int dof_two = NodeList[node_two - 1].bcode[2];
				 unsigned int dof_three = NodeList[node_three - 1].bcode[2];
				 unsigned int dof_four = NodeList[node_four - 1].bcode[2];
				 if (dof_one)
					 Force[dof_one - 1] += -Element.GetGravity() / 4;
				 if (dof_two)
					 Force[dof_two - 1] += -Element.GetGravity() / 4;
				 if (dof_three)
					 Force[dof_three - 1] += -Element.GetGravity() / 4;
				 if (dof_four)
					 Force[dof_four - 1] += -Element.GetGravity() / 4;
			 }
		 }
		 	break;
		 case Beam:
			 {
				 unsigned int NUME = ElementGrp.GetNUME();
				 double* ptr_force = new double[6];
				 clear(ptr_force, 6);
				 for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
				 {
					 CElement& Element = ElementGrp[Ele];
					 Element.GravityCalculation(ptr_force);
					 CNode** node_ = Element.CElement::GetNodes();
					 double gravity = Element.GetGravity();
					 for (int i = 0; i < 2; i++)
					 {
						 for  (int j = 2; j < 5; j++)
						 {
							 if (NodeList[node_[i]->NodeNumber - 1].bcode[j]) 
							 {
								 Force[NodeList[node_[i]->NodeNumber - 1].bcode[j] - 1] += ptr_force[i * 3 + j - 2];
							 }
						 }
					 }
				 }
			 }
			 break;

		 case Plate: //Plate element
		 {
			 unsigned int NUME = ElementGrp.GetNUME();
			 double* ptr_force = new double[12];
			 clear(ptr_force, 12);
			 for (unsigned int Ele = 0; Ele < NUME; Ele++) {
				 CElement& Element = ElementGrp[Ele];
				 Element.GravityCalculation(ptr_force);
				 CNode** node_ = Element.CElement::GetNodes();
				 for (unsigned int i = 0; i < 4; i++) {
					 for (unsigned int j = 2; j < 5; j++)
						 if (NodeList[node_[i]->NodeNumber - 1].bcode[j])
							 Force[NodeList[node_[i]->NodeNumber - 1].bcode[j] - 1] += ptr_force[i * 3 + j - 2];
				 }
			 }
		 }
		 break;

		case Shell:
		{
			unsigned int NUME = ElementGrp.GetNUME();
			double* ptr_force=new double[12];
			clear(ptr_force,12);
			for(unsigned int Ele=0;Ele<NUME;Ele++){
				CElement& Element = ElementGrp[Ele];
				Element.GravityCalculation(ptr_force);
				CNode** node_ = Element.CElement::GetNodes();
				for(unsigned int i=0;i<4;i++){
					for(unsigned int j=2;j<5;j++)
						if (NodeList[node_[i]->NodeNumber - 1].bcode[j]) 
							 Force[NodeList[node_[i]->NodeNumber - 1].bcode[j] - 1] += ptr_force[i * 3 + j - 2];
				}
			}
		}
			 break;
		case T3:
			;
			break;
		case H8:
			{
				unsigned int NUME = ElementGrp.GetNUME();
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = ElementGrp[Ele];
					unsigned int nen=Element.GetNEN();
					double* EG = new double[nen];
					Element.GravityCalculation(EG);
					for (unsigned int nn = 0; nn < nen; nn++) 
					{
						unsigned int dof = Element.GetNodes()[nn] -> bcode[2];
						if (dof)
						{
							Force[dof - 1] -= EG[nn]; 
						}
					}
					delete [] EG;
				}
			}
			break;
		case H8R:
			{
				unsigned int NUME = ElementGrp.GetNUME();
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = ElementGrp[Ele];
					unsigned int nen=Element.GetNEN();
					double* EG = new double[nen];
					Element.GravityCalculation(EG);
					for (unsigned int nn = 0; nn < nen; nn++) 
					{
						unsigned int dof = Element.GetNodes()[nn] -> bcode[2];
						if (dof)
						{
							Force[dof - 1] -= EG[nn]; 
						}
					}
					delete [] EG;
				}
			}
			break;
		}
	}
}

//	Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//	and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
//	Allocate for global force/displacement vector
	Force = new double[NEQ];
    clear(Force, NEQ);

//  Create the banded stiffness matrix
#ifdef _PARDISO_
	CSRStiffnessMatrix = new CSparseMatrix<double>(NEQ);
	CalculateCSRColumns();
	CSRStiffnessMatrix->Allocate();
#else
	StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
//	Calculate column heights
	CalculateColumnHeights();
//	Calculate address of diagonal elements in banded matrix
	StiffnessMatrix->CalculateDiagnoalAddress();

//	Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();
	COutputter* Output = COutputter::Instance();
	Output->OutputTotalSystemData();
#endif // _PARDISO_

    




	
}
