#ifndef GLOBALS_H
#define GLOBALS_H

#include <climits>

/**************************************************
* This file contains all the global variable definitions and marcos, 
* along with several type definitions.
***************************************************/

typedef int CellNrType; // it could be long for really large inputs...
typedef vector<int> MatrixListType;


namespace Globals 
{
	const int BIG_INT = INT_MAX;			// be careful, could be too small compared to # of cubes

	double reduction_threshold = 0.0;		// if the persistence exceeds this threshold, then a further reduction
															// is performed tocompute the shortest cycle.

	bool use_optimal_alg = false;				// whether to perform the cycle optimization algorithm

	int which_alg = 0;								// 0: proposed heuristic-based annotation algorithm
															// 1: classical annotation algorithm (exhaustive search)

	int max_dim = 2;								// the maximum dimension to be computed

	string inputFileName;							// input data file name

	string memoryFileName_HeuristicAlg = "Memory_Footprint_HeuristicAlg.txt";
	string memoryFileName_ClassicalAlg = "Memory_Footprint_ClassicalAlg.txt";

	enum Algorithm
	{
		HEURISTIC_BASED_ALG = 0,
		CLASSICAL_ALG = 1,
	};

	enum FileType
	{
		IMAGE_DATA = 0,
		DENSE_DISTANCE_MATRIX = 1,
		GENERAL_SIMPLICIAL_COMPLEX = 2
	};
}


#endif // !GLOBALS_H

