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

#include "Outputter.h"

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

	double density; //!< density 

public:

//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output, unsigned int mset) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};



//!	Material class for 4Q element
class C4QMaterial : public CMaterial
{
public:

	double poisson;	//!Poisson ratio of a 4Q element
	double etype;//!element type of plane strain (2) or plane stress(1) or The column symmetry (3)
	double thick; // the thickness of element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};


class CH8Material : public CMaterial
{
public:

	double Nu, G, Lam, Rou;


public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};


class CBeamMaterial : public CMaterial
{
public:

	double mu;        // Poisson ratio
	double width;     // width of rectangle
	double height;    // height of rectangle
	double t_side;    // flank thickness
	double t_uplow;   // upper and lowwer surface thickness
	double n_x;       // x component of y' axis
	double n_y;       // y component of y' axis
	double n_z;       // z component of y' axis


public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};


class CPlateMaterial : public CMaterial
{
public:

	double poisson;	//!Poisson ratio of a 4Q element
	double thick; // the thickness of element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};


class CInfiMaterial : public CMaterial
{
public:

	double poisson;	//!Poisson ratio of a 4Q element
	double etype;//!element type of strain and stress

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};
