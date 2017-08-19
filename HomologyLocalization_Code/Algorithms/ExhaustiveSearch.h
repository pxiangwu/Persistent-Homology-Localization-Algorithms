#ifndef EXHAUSTIVE_SEARCH_H
#define EXHAUSTIVE_SEARCH_H

#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <ctime>

#include "../Filtration/CubicalFiltration.h"
#include "../PersistenceIO.h"
#include "../Algorithms/DijkstraShortestPath.h"
#include "AnnotatingEdges.h"
#include "../BitSet.h"
#include "AStar.h"
#include "../Globals.h"

#include "../External/Mem_usage.h"


// convert annotation to integer
int convertAnnotation2Int(const BitSet & vertexAnnotation)
{
	double pos = 0;
	for (int i = 0; i < vertexAnnotation.getBitSize(); ++i)
	{
		if (vertexAnnotation.checkBit(i) == true)
			pos += pow(2.0, (double)i);
	}
	return (int)pos;
}


// convert integer to annotation
BitSet convertInt2Annotation(int pos, int BettiNum)
{
	BitSet resAnnotation(BettiNum);
	int idx = 0;
	int bit;

	while (pos != 0)
	{
		bit = pos % 2;
		if (bit == 1)
			resAnnotation.set(idx);

		pos = pos >> 1;
		idx++;
	}
	return resAnnotation;
}

// construct the whole covering graph
void constructCoveringGraph(const map<pair<int, int>, BitSet> & edgeAnnotations, int low,
	int vertexNum, const vector<MatrixListType> & cell2v_list, adjacency_list_t & resCoveringGraph)
{
	int BettiNum = (edgeAnnotations.cbegin()->second).getBitSize();
	int numCopy = pow(2.0, double(BettiNum));

	resCoveringGraph.clear();
	resCoveringGraph.resize(vertexNum * numCopy);

	MatrixListType edge;
	pair<int, int> key;
	int ptr_1, ptr_2;
	int vCopyPos;
	BitSet tempAnnotation(BettiNum);
	map<pair<int, int>, BitSet>::const_iterator it_annotation;

	for (int i = 0; i < low; ++i)
	{
		edge = cell2v_list[i];

		ptr_1 = edge[0];
		ptr_2 = edge[1];

		if (ptr_1 > ptr_2)
			SWAP(ptr_1, ptr_2);

		key.first = ptr_1;
		key.second = ptr_2;

		// check if it is sentinel edge
		it_annotation = edgeAnnotations.find(key);

		if (it_annotation == edgeAnnotations.end()) // Not a sentinel edge
		{
			for (int j = 0; j < numCopy; ++j)
			{
				resCoveringGraph[ptr_1 + j * vertexNum].push_back(neighbor(ptr_2 + j * vertexNum));
				resCoveringGraph[ptr_2 + j * vertexNum].push_back(neighbor(ptr_1 + j * vertexNum));
			}
		}
		else // A sentinel edge
		{
			for (int j = 0; j < numCopy; ++j)
			{
				tempAnnotation = convertInt2Annotation(j, BettiNum);
				tempAnnotation ^= it_annotation->second;
				vCopyPos = convertAnnotation2Int(tempAnnotation);

				resCoveringGraph[ptr_1 + j * vertexNum].push_back(neighbor(ptr_2 + vCopyPos * vertexNum));
				resCoveringGraph[ptr_2 + vCopyPos * vertexNum].push_back(neighbor(ptr_1 + j * vertexNum));
			}
		}
	}
}

// the algorithm of exhautive search
void ExhaustiveSearch(const adjacency_list_t & coveringGraph, const MatrixListType & inputCycle, const vector<MatrixListType> & cell2v_list,
	const map<pair<int, int>, BitSet> & edgeAnnotations, int vertexNum, const std::map<std::pair<int, int>, int> & edgeMap,
	MatrixListType & resShortestCycle)
{
	resShortestCycle.clear();

	assert(!edgeAnnotations.empty());
	int BettiNum = (edgeAnnotations.cbegin()->second).getBitSize();


	BitSet targetAnnotation(BettiNum); // the target annotation we should reach finally
	computeCycleAnnotation(inputCycle, cell2v_list, edgeAnnotations, targetAnnotation);

	int low = inputCycle.back();
	MatrixListType lowEdge;
	lowEdge = cell2v_list[low];
	int source = lowEdge[0];
	int target = lowEdge[1];
	if (source > target)
		SWAP(source, target);
	pair<int, int> key(source, target);
	map<pair<int, int>, BitSet>::const_iterator it_annotation;
	it_annotation = edgeAnnotations.find(key);
	if (it_annotation != edgeAnnotations.end()) // the pivot edge is a sentinel edge
		targetAnnotation ^= it_annotation->second; // exclude the pivot edge


	int targetCopy = convertAnnotation2Int(targetAnnotation);
	target = target + targetCopy * vertexNum;

	std::vector<vertex_t> previous; // used for backtracing the shortest path
	std::vector<weight_t> min_distance; // used for recording the shortest distance

	// Apply Dijkstra's algorithm
	DijkstraComputePaths(source, target, coveringGraph, min_distance, previous);

	// Backtrace the shortest path
	std::list<vertex_t> path = DijkstraGetShortestPathTo(target, previous);
	int ptr_1, ptr_2, resEdge;
	for (std::list<vertex_t>::iterator it = path.begin(); it != std::prev(path.end()); it++)
	{
		ptr_1 = (*it) % vertexNum;
		ptr_2 = (*(std::next(it))) % vertexNum;

		if (ptr_1 > ptr_2) // sort
			SWAP(ptr_1, ptr_2);

		std::pair<int, int> endPtrs(ptr_1, ptr_2);
		resEdge = edgeMap.at(endPtrs);
		resShortestCycle.push_back(resEdge);
	}

	// Add the pivot edge
	resShortestCycle.push_back(low);
	mysort(resShortestCycle);

	std::size_t physMemUsed; // for monitoring the memory footprint
	ofstream memoryFile(Globals::memoryFileName_ClassicalAlg, ios::out | ios::app);
	clock_t currTime;

	physMemUsed = getPeakRSS() >> 20;
	currTime = clock();
	// memoryFile << (float)(currTime - startTime) / CLOCKS_PER_SEC << " " << physMemUsed / 1024.0 / 1024.0 << endl;
	memoryFile << "Betti Number: " << BettiNum << "\t " << "Memory Footprint: " << physMemUsed << " (MB)" << endl;

	memoryFile.close();

	// free memory
	previous.clear();
	min_distance.clear();
	path.clear();
}


// interface for running classical annotation-based algorithm
template<int arrayDim, int vertexDim = arrayDim>
void reduceND_ExhaustiveSearch(blitz::Array<double, arrayDim> *phi, const vector<blitz::TinyVector<int, vertexDim>> & vList,
	const vector<int> & lowerCellList, const vector<CellNrType> & upperCellList, vector<MatrixListType> & boundaryMatrix,
	const std::map<std::pair<int, int>, int> & edgeMap, const vector<MatrixListType> &cell2v_list, int vertexNum,
	vector<int> &low_array)
{
	clock_t startClock, endClock;

	cout << "--- Using Classical Annotation Algorithm (Exhaustive Search) ---" << endl;
	ofstream memoryFile(Globals::memoryFileName_ClassicalAlg, ios::out | ios::trunc);
	memoryFile.close();

	std::map<int, MatrixListType> resCycles;
	for (int test = 0; test < boundaryMatrix.size(); test++)
	{
		if (boundaryMatrix[test].empty())
			continue;

		double birthTime, deathTime;
		double pers = computePersistence<arrayDim, vertexDim>(phi, vList, lowerCellList, upperCellList, boundaryMatrix, test, birthTime, deathTime);
		if (pers <= Globals::reduction_threshold)
			continue;

		// We first compute the annotations of all edges
		cout << "---------------------------------------" << endl;
		cout << "Compute edge annotations ..." << endl;
		map<pair<int, int>, BitSet> edgeAnnotations;
		computeAnnotations(boundaryMatrix, edgeMap, low_array, cell2v_list, test, vertexNum, edgeAnnotations);

		startClock = clock();
		// Construct convering graph
		cout << "Construct covering graph ..." << endl;
		adjacency_list_t coveringGraph;
		constructCoveringGraph(edgeAnnotations, boundaryMatrix[test].back(), vertexNum, cell2v_list, coveringGraph);

		cout << "Apply Exhaustive Search algorithm ..." << endl;
		MatrixListType resCycle;
		ExhaustiveSearch(coveringGraph, boundaryMatrix[test], cell2v_list, edgeAnnotations, vertexNum, edgeMap, resCycle);
		cout << "Size before: " << boundaryMatrix[test].size() << endl;
		cout << "Size after: " << resCycle.size() << endl;

		endClock = clock();
		cout << "Time consumed (s): " << (float)(endClock - startClock) / CLOCKS_PER_SEC << endl;
		cout << "---------------------------------------" << endl;

		resCycles.insert({ test, resCycle });
	}
	for (const auto & cycle : resCycles)
	{
		boundaryMatrix[cycle.first] = cycle.second;
	}
}

#endif // !PAPER2012_H

