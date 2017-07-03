/**************************************************************************
* Copyright 2017 Rutgers University and CUNY Queens College.
* Code by Pengxiang Wu based on the implementations of Hubert Wagner and Chao Chen.
* 
* License: GPLv3. (https://www.gnu.org/licenses/gpl.html)
*
* The SOFTWARE PACKAGE provided in this page is provided "as is", without any guarantee
* made as to its suitability or fitness for any particular use. It may contain bugs, so use of 
* this tool is at your own risk. We take no responsibility for any damage of any sort that may
* unintentionally be caused through its use.
*
* If you have any questions, please contact Pengxiang Wu (pxiangwu@gmail.com)
**************************************************************************/
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <deque>
#include <cstring>
#include <ctime>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
using namespace std;

#include "PersistenceIO.h"
#include "Debugging.h"
#include "Filtration/CubicalFiltration.h"

#include "InputRunner.h"
#include "InputFileInfo.h"
#include "PersistentPair.h"

#include "PersistenceCalculator.h"
#include "PersistenceCalcRunner.h"
#include "Globals.h"
#include "ParseCommandLine.h"

void runPersistenceHomology(const InputFileInfo &);
void debugStart();
void debugEnd();

int main(int argc, const char* argv[])
{
	time_t startTime, endTime;		

	parseCommandLine(argc, argv);

	InputFileInfo input_file_info(Globals::inputFileName);
	
	debugStart(); // debug settings
	time(&startTime);	

	// Run persistence homology algorithm
	runPersistenceHomology(input_file_info);

	time(&endTime);

	double ellapsed1 = difftime(endTime, startTime);
	cout << "All calculations computed in " << setprecision(3) << ellapsed1 / 60.0 << " Min" << endl;
	
	debugEnd();

	return 0;
}


void debugStart()
{
	string logFile = "log.txt";
	string errorFile = "error.txt";
	DebuggerClass::init(false, logFile, errorFile);

	int max_pers_pts = Globals::BIG_INT; // the maximal number of persistence pairs allowed

	fstream filestr;
	filestr.open(logFile.c_str(), fstream::out | fstream::trunc);

	filestr << "################################################" << endl;
	filestr << "Start computing persistence" << endl;
	filestr << "Input file = \'" << Globals::inputFileName << "\'" << endl;
	filestr << "Maximal number of persistence points allowed = " << max_pers_pts << endl;
	filestr.close();
}

void debugEnd()
{
	DebuggerClass::finish();
}