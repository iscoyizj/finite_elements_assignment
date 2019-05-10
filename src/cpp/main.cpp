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
#include "Bar.h"
#include "Outputter.h"
#include "Clock.h"
#include "Outputterplot1.h"

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName\n";
		exit(1);
	}

	string filename(argv[1]);
    size_t found = filename.find_last_of('.');

    // If the input file name is provided with an extension
    if (found != std::string::npos) 
	{
        if (filename.substr(found) == ".dat")
            filename = filename.substr(0, found);
        else 
		{
            // The input file name must has an extension of 'dat'
            cout << "*** Error *** Invalid file extension: "
                 << filename.substr(found+1) << endl;
            exit(1);
        }
    }

    string InFile = filename + ".dat";
	string OutFile = filename + ".out";
	string PostFile = filename + "_post.out";

	CDomain* FEMData = CDomain::Instance();

    Clock timer;
    timer.Start();

//  Read data and define the problem domain
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();

//  Allocate global vectors and matrices, such as the Force, ColumnHeights,
//  DiagonalAddress and StiffnessMatrix, and calculate the column heights
//  and address of diagonal elements
	FEMData->AllocateMatrices();
    
//  Assemble the banded gloabl stiffness matrix
	FEMData->AssembleStiffnessMatrix();

// Gravity
	FEMData->Gravity();
    
    double time_assemble = timer.ElapsedTime();

//  Solve the linear equilibrium equations for displacements
#ifdef _PARDISO_
	CSRSolver* SRSolver = new CSRSolver(FEMData->GetCSRStiffnessMatrix());
#else
	CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
	Solver->LDLT();
#endif // _PARDISO_
//  Perform L*D*L(T) factorization of stiffness matrix
    

    COutputter* Output = COutputter::Instance();

#ifdef _DEBUG_
#ifndef _PARDISO_
 Output->PrintStiffnessMatrix();
#endif // !_PARDISO_
#endif
        
//  Loop over for all load cases
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
//      Assemble righ-hand-side vector (force vector)
        
		FEMData->AssembleForce(lcase + 1);
		
//      Reduce right-hand-side force vector and back substitute
        
#ifdef  _PARDISO_
		
		SRSolver->solve(FEMData->GetForce(), 1);  
#else
		//      Reduce right-hand-side force vector and back substitute
		Solver->BackSubstitution(FEMData->GetForce());
#endif //  _PARDISO_

        
#ifdef _DEBUG_
        Output->PrintDisplacement(lcase);
		
#endif
            
        Output->OutputNodalDisplacement(lcase);
    }

    double time_solution = timer.ElapsedTime();

//  Calculate and output stresses of all elements
	Output->OutputElementStress();
	
	COutputterplot1* Outputplot1 = COutputterplot1::Instanceplot1(PostFile);
	Outputplot1->OutputElementStress();

    double time_stress = timer.ElapsedTime();
    
    timer.Stop();
    
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     TIME FOR INPUT PHASE = " << time_input << endl
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
            << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl
            << "     T O T A L   S O L U T I O N   T I M E = " << time_stress << endl;

	
	return 0;
}
