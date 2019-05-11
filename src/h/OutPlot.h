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

#include <fstream>
#include <iostream>

using namespace std;

//! Outputer class is used to output results
class COutPlot
{
private:

//!	File stream for output
	ofstream OutputFile;

protected:

//!	Constructor
	COutPlot(string FileName);

//!	Designed as a single instance class
	static COutPlot* _instance;

public:

//!	Return pointer to the output file stream
	inline ofstream* GetOutputFile() { return &OutputFile; }

//!	Return the single instance of the class
	static COutPlot* Instance(string FileName = " ");

//!	Output current time and date
	void PrintTime(const struct tm * ptm, COutPlot& output);

//!	Output logo and heading 
	void OutputHeading();

//!	Output nodal point data
	void OutNode();
	
	void OutputElementInfo(unsigned int n, unsigned int sum);

	void OutputEleType(double n, double nnd);

	void OutputNodalDisplacement(unsigned int lcase);

	void StressHead(unsigned int lcase, unsigned int n);

	void ElementStress(double stress);

	//! Overload the operator <<
	template <typename T>
	COutPlot& operator<<(const T& item) 
	{
//		std::cout << item;
		OutputFile << item;
		return *this;
	}

	typedef std::basic_ostream<char, std::char_traits<char> > CharOstream;
	COutPlot& operator<<(CharOstream& (*op)(CharOstream&)) 
	{
//		op(std::cout);
		op(OutputFile);
		return *this;
	}

};