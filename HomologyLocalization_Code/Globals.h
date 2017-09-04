#ifndef GLOBALS_H
#define GLOBALS_H

#include <climits>
#include "BitSet.h"

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


// structure of covering graph node
struct cgNode
{
	cgNode() {}

	cgNode(int BettiNum)
	{
		sumAnnotation = BitSet(BettiNum);
	}

	bool operator < (const cgNode & rhs) const // min-heap, with tie-breaking
	{
		if (fScore > rhs.fScore)
			return true;
		else if (fScore < rhs.fScore)
			return false;
		else if (fScore == rhs.fScore)
		{
			if (gScore != rhs.gScore)
				return gScore < rhs.gScore;
			else
				return true;
		}
	}

	int vertex = -1;								// the vertex index
	double fScore = 0.0;						// fScore = the length of walked path + the length of estimated remaining path
	double gScore = 0.0;						// the length of walked path;
	BitSet sumAnnotation;					// sum of annotations of encountered edges when searching along some path
	vector<pair<int, int>> previous;		// used for backtracing the shortest path
};


// Hasher for pair<int, BitSet>
struct KeyHasher
{
	std::size_t operator () (const std::pair<int, BitSet> & key) const
	{
		size_t res = 17;
		int first = key.first;
		int second = key.second.convertBitsToInt();

		res = res * 31 + hash<int>()(first);
		res = res * 31 + hash<int>()(second);

		return res;
	}
};


#endif // !GLOBALS_H

