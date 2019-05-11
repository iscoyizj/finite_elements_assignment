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
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << Area << density << endl;
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

bool CH8Material::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> Nu >> Rou;	// Young's modulus and section area
	G = E / 2 / (1 + Nu);
	Lam = Nu * E / (1+Nu) / (1-2*Nu);

	return true;
}

void CH8Material::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << Nu << setw(16) << G << setw(16) << Lam << setw(16) << Rou << endl;
}

bool CBeamMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> mu >> density                        // material properties
		  >> width >> height                           // geometry properties
		  >> t_side >> t_uplow                         // thickness of the beam's side 
		  >> n_x >> n_y >> n_z;	                       // the cirection of y axis


	return true;
}

//	Write material data to Stream
void CBeamMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << mu << setw(16) << density
	<< setw(16) << width << setw(16) << height << setw(16) << t_side << setw(16) << t_uplow << endl;
}

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

	Input >> E >> poisson >> thick >> density;	// Young's modulus,Poisson ratio and element type

	return true;
}

void CPlateMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << poisson << setw(16) << thick << setw(16) << density << endl;
}


bool CInfiMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set
  Input >> E >> poisson >> etype;	// Young's modulus,Poisson ratio and element type

	return true;
}

bool CShellMaterial::Read(ifstream& Input, unsigned int mset){
	Input >> nset;	// Number of property set;
  if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input>>E>>nu>>density>>thick;
	return true;
}

void CShellMaterial::Write(COutputter& output, unsigned int mset){
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << nu << setw(16) << density << setw(16)<< thick << endl;
}


//	Write material data to Stream
void CInfiMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << poisson << endl;
}


bool CSubparaMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cout << "*** Error *** Material sets must be inputted in order !" << endl
			<< "   Expected set : " << mset + 1 << endl
			<< "   Provided set : " << nset << endl;
		return false;
	}

	Input >> E >> poisson;	// Young's modulus and Poisson's ratio
	return true;
}

//	Write Subparametric material data to Stream OutputFile
void CSubparaMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << poisson << endl;
}
