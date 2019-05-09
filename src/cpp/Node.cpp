/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>
#include<string>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
	Rotation_Flag = 0;
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;
	XYZ[3] = 0;
	XYZ[4] = 0;
	XYZ[5] = 0;
	NodeNumber = 0;
    
    bcode[0] = 0;	// Boundary codes
    bcode[1] = 0;
    bcode[2] = 0;
	bcode[3] = 1;
	bcode[4] = 1;
	bcode[5] = 1;

	stress_node[0] = 0.0;
	stress_node[1] = 0.0;
	stress_node[2] = 0.0;
	stress_node[3] = 0.0;
	stress_node[4] = 0.0;
	stress_node[5] = 0.0;

	
};

// return total count of non-blank args in string
int Count_for_rotation_arg(std::string getstr)
{
	int count = 0;
	bool on_num = false;
	for (unsigned i = 0; i < getstr.length(); ++i)
	{
		if (getstr[i] != ' ' && getstr[i] != '\t')
		{
			if (!on_num)
			{
				on_num = true;
				count++;
			}
		}
		else
		{
			on_num = false;
		}
	}
	return count;
}

//	Read element data from stream Input
bool CNode::Read(ifstream& Input, unsigned int np)
{
	unsigned int N;

	Input >> N;	// node number
	if (N != np + 1) 
	{
		cerr << "*** Error *** Nodes must be inputted in order !" << endl 
			 << "   Expected node number : " << np + 1 << endl
			 << "   Provided node number : " << N << endl;

		return false;
	}

	NodeNumber = N;
	// Read the dataline
	string NodeInfo, ModiNodeInfo;
	getline(Input, NodeInfo);

	int input_count = Count_for_rotation_arg(NodeInfo);

	// Determine the input format:
	//     if allow rotation for plate and beam , input_count is 9;
	//     for elements with no rotation at nodes, input_count is 6.
	

	
	if (input_count == 9)
	{
		Rotation_Flag = 1;
		sscanf(NodeInfo.c_str(), "%d%d%d%d%d%d%lf%lf%lf",
			bcode, bcode + 1, bcode + 2,
			bcode + 3, bcode + 4, bcode + 5,
			XYZ, XYZ + 1, XYZ + 2);
	}
	else if (input_count == 6)
	{
		Rotation_Flag = 0;
		sscanf(NodeInfo.c_str(), "%d%d%d%lf%lf%lf",
			bcode, bcode + 1, bcode + 2,
			XYZ, XYZ + 1, XYZ + 2);
	}
	else
	{
		cerr << "*** Error *** NodeInfos must be inputted in the correct format! " << endl
			<< "  Present Number of Nodeinfos: " << input_count << endl
			<< "  Please input 6 or 9" << endl;
		return false;
	}

	return true;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output, unsigned int np)
{
	output << setw(9) << np + 1 << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
		   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output, unsigned int np)
{
	output << setw(9) << np+1 << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, unsigned int np, double* Displacement)
{
	output << setw(5) << np + 1 << "        ";

	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}
