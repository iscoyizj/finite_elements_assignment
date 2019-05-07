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
#include <string>
using namespace std;

//! Outputer class is used to output results
class COutputterplot1
{
private:

//!	File stream for output
	ofstream OutputFileplot1;

protected:

//!	Constructor
	COutputterplot1(string FileName);

//!	Designed as a single instance class
	static COutputterplot1* _instanceplot1;

public:

//!	Return pointer to the output file stream
	inline ofstream* GetOutputFileplot1() { return &OutputFileplot1; }

//!	Return the single instance of the class
	static COutputterplot1* Instanceplot1(string FileName = " ");

//!	Output logo and heading 
	void OutputHeading();

//!	Output element stresses 
	void OutputElementStress();

//! Overload the operator <<
	template <typename T>
	COutputterplot1& operator<<(const T& item) 
	{
		std::cout << item;
		OutputFileplot1 << item;
		return *this;
	}

	typedef std::basic_ostream<char, std::char_traits<char> > CharOstream;
	COutputterplot1& operator<<(CharOstream& (*op)(CharOstream&)) 
	{
		op(std::cout);
		op(OutputFileplot1);
		return *this;
	}


};
