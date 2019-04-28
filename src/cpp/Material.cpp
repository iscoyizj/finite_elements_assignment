/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> Area>> density;	// Young's modulus and section area and density

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << Area << density<< endl;
}





//	Read material data from stream Input
bool C4QMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> poisson >> etype >> density >> thick;	// Young's modulus,Poisson ratio and element type

	return true;
}


//	Write material data to Stream
void C4QMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << poisson << setw(16) << density << setw(16)<< thick << endl;
}

//	Read material data from stream Input
bool CPlateMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> poisson >> thick>> density;	// Young's modulus,Poisson ratio and element type

	return true;
}

void C4QMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << poisson << setw(16) << thick<< setw(16) << density << endl;
}
