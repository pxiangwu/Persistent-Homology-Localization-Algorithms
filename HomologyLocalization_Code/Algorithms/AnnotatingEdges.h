#ifndef ANNOTATING_EDGES_H
#define ANNOTATING_EDGES_H

#include <iostream>
#include <map>
#include <vector>
#include <stack>
#include <queue>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <mutex>

#include "../PersistenceIO.h"
#include "../Algorithms/DijkstraShortestPath.h"
#include "../STLUtils.h"
#include "../BitSet.h" // data structure for handling binary annotation
#include "../Globals.h"
 
std::mutex thread_mutex;

using namespace std;

/********************************************************************
* Description:	Given a graph, this function computes an arbitrary spanning tree
* Parameters:	
* - InputGraph:				the input graph
* - inputGraphEdges:		the edge set of input graph
* - resSpanningTree:		the result spanning tree (organized with adjacency list)
* - resSentinelEdges:		the result set containing sentinel edges (stored in a vector)
********************************************************************/
void computeSpanningTree(const adjacency_list_t & inputGraph, vector<pair<int, int>> & inputGraphEdges,
	adjacency_list_t & resSpanningTree, vector<pair<int, int>> & resSentinelEdges)
{
	std::stack<int> stk; // stack used for DFS
	size_t numV = inputGraph.size();
	std::vector<bool> isVisited(numV, false);
	int currentV, neighborV;
	int neighborSZ;
	int i;
	neighbor nb(0, 0);

	std::vector<std::pair<int, int>> treeEdges;
	resSentinelEdges.clear();
	resSpanningTree.clear();
	resSpanningTree.resize(numV);

	std::random_device rd;
	std::mt19937 g(rd());
	vector<int> order(4);
	vector<int> nodesOrder(numV);

	int temp = { -1 };
	std::generate(nodesOrder.begin(), nodesOrder.end(), [&temp] {return ++temp;});
	std::shuffle(nodesOrder.begin(), nodesOrder.end(), g);

	for (int h = 0; h < numV; h++) // DFS
	{
		i = nodesOrder[h];

		if (isVisited[i] == true || inputGraph[i].empty())
			continue;

		isVisited[i] = true;
		stk.push(i);

		while (!stk.empty())
		{
			currentV = stk.top();
			stk.pop();

			const std::vector<neighbor> & neighbors = inputGraph[currentV];
			neighborSZ = neighbors.size();

			int n = { -1 };
			order.resize(neighborSZ);
			std::generate(order.begin(), order.end(), [&n] {return ++n;});
			std::shuffle(order.begin(), order.end(), g);

			for (int h = 0; h < neighborSZ; ++h) // check all its neighbors
			{
				nb = neighbors[order[h]];

				neighborV = nb.target;
				if (isVisited[neighborV] == false)
				{
					isVisited[neighborV] = true;
					stk.push(neighborV);

					resSpanningTree[currentV].push_back(neighbor(neighborV)); // create a tree edge
					resSpanningTree[neighborV].push_back(neighbor(currentV));

					if (currentV > neighborV) // sort the endpoints into nondecreasing order
						treeEdges.push_back(std::make_pair(neighborV, currentV));
					else
						treeEdges.push_back(std::make_pair(currentV, neighborV));
				}
			}
		}
	} // end DFS

	// Next, compute the sentinel edges
	std::sort(treeEdges.begin(), treeEdges.end());

	std::sort(inputGraphEdges.begin(), inputGraphEdges.end());
	std::vector<std::pair<int, int>>::iterator it = std::unique(inputGraphEdges.begin(), inputGraphEdges.end());
	inputGraphEdges.resize(std::distance(inputGraphEdges.begin(), it));

	std::set_difference(inputGraphEdges.begin(), inputGraphEdges.end(), treeEdges.begin(), treeEdges.end(), std::back_inserter(resSentinelEdges));
}


/********************************************************************
* Description:	Given a sentinel edge, this function returns its corresponding unique
						sentinel cycle in the spanning tree. We use BFS to find such a cycle.
* Parameters:
* - inputSpanningTree:	the input spanning tree
* - edge:						the input sentinel edge
* - edgeMap:					a map, mapping two endpoints to an edge
* - resSentinelCycle:		the result sentinel cycle corresponding to the input edge
********************************************************************/
void computeSentinelCycle(const adjacency_list_t & inputSpanningTree, const pair<int, int> & edge, const map<pair<int, int>, int> & edgeMap,
	MatrixListType & resSentinelCycle)
{
	resSentinelCycle.clear();

	size_t numV = inputSpanningTree.size();
	std::queue<int> searchQ; // used for BFS
	std::map<int, int> previous; // used for backtracing the cycle path
	std::vector<bool> isVisited(numV, false);

	int source = edge.first;
	int target = edge.second;
	int neighborV, currentV;
	bool hasFound = false; // has found the target?

	// begin BFS
	searchQ.push(source);
	isVisited[source] = true;
	previous[source] = target;
	while (!hasFound)
	{
		currentV = searchQ.front();
		searchQ.pop();

		const std::vector<neighbor> & neighbors = inputSpanningTree[currentV];
		for (const auto & nb : neighbors)
		{
			neighborV = nb.target;
			if (isVisited[neighborV] == false)
			{
				isVisited[neighborV] = true;
				previous[neighborV] = currentV;
				searchQ.push(neighborV);

				if (neighborV == target)
				{
					hasFound = true;
					break;
				}
			}
		}
	} // end while

	// construct the sentinel cycle
	int ptr_1, ptr_2;
	int resEdge;
	currentV = target;
	do
	{
		neighborV = previous[currentV];
		if (currentV > neighborV) // sort the endpoints into nondecreasing order
		{
			ptr_1 = neighborV;
			ptr_2 = currentV;
		}
		else
		{
			ptr_1 = currentV;
			ptr_2 = neighborV;
		}

		resEdge = edgeMap.at(std::make_pair(ptr_1, ptr_2)); // query the associated edge
		resSentinelCycle.push_back(resEdge);

		currentV = neighborV;
	} while (currentV != target);

	mysort(resSentinelCycle);
}


/********************************************************************
* Description:	Given a persistent homology class, compute the corresponding betti number
* Parameters:
* - redBoundary:			reduced boundary matrix
* - birth:						the given birth time, which is a number from filtrationOrder
* - death:						the given death time
* - mapColorColumnIdx:	mapping color columns into consecutive indices. This is a by-product
********************************************************************/
int computeBettiNumber(const vector<MatrixListType> & redBoundary, int birth, int death, map<int, int> & mapColorColumnIdx)
{
	mapColorColumnIdx.clear();

	size_t cnt = 0;
	size_t low;
	size_t bettiNum = 0;
	size_t boundarySize = redBoundary.size();
	for (size_t i = 0; i < boundarySize; i++)
	{
		if (redBoundary[i].empty())
			continue;

		low = redBoundary[i].back();
		if (low <= birth && i >= death)
		{
			bettiNum++;
			mapColorColumnIdx[i] = cnt++; // record this column, convert it to consecutive index
		}
	}
	return bettiNum;
}


/********************************************************************
* Description:	Given a persistence homology, which is specified by parameter column,
						compute its persistence.
* Parameters:
* - phi:							the filter function
* - vList:						a list storing the coordinates of vertices
* - lowerCellList:			a list storing the maxValues of cells
* - upperCellList:			a list stroing the maxValues of cells, whose dimension is greater
									than that of lowerCellList by 1
* - boundaryMatrix:		the reduced boundary matrix
* - column:					specifies which homology class
********************************************************************/
template<int arrayDim, int vertexDim = arrayDim>
double computePersistence(blitz::Array<double, arrayDim> * phi, const vector<blitz::TinyVector<int, vertexDim>> &vList,
	const vector<int> & lowerCellList, const vector<CellNrType> & upperCellList, 
	const vector<MatrixListType> & boundaryMatrix, int column, double & birthTime, double & deathTime)
{
	int vBirth = lowerCellList[boundaryMatrix[column].back()];
	int vDeath = upperCellList[column];

	double tmp_death = (*phi)(vList[vDeath]);
	double tmp_birth = (*phi)(vList[vBirth]);
	double pers = tmp_death - tmp_birth;

	birthTime = tmp_birth;
	deathTime = tmp_death;

	return pers;
}


/********************************************************************
* Description:	thread for computing annotation of a given sentinel edge
* Parameters:
* - sentinelEdges:			the set of sentinel edges
* - spanningTree:			spanning tree
* - other parameters are self-explanatory
********************************************************************/
void threadComputeAnnotation(const vector<pair<int, int>> & sentinelEdges,
	const adjacency_list_t & spanningTree, const map<pair<int, int>, int> & edgeMap,
	const vector<int> & low_array, int death, int bettiNum, const vector<MatrixListType> & redBoundary,
	const map<int, int> & mapColorColumnIdx, map<pair<int, int>, BitSet> & resEdgeAnnotations)
{
	MatrixListType sentinelCycle;
	BitSet annotation(bettiNum);
	int low;

	for (const auto & edge : sentinelEdges)
	{
		annotation.reset(); // reinitialized as zeros
		sentinelCycle.clear();
		computeSentinelCycle(spanningTree, edge, edgeMap, sentinelCycle);

		// performe reduction on this sentinel cycle using all the colored columns, thereby obtaining the annotation.
		while (!sentinelCycle.empty()) // we continue the reduction until it is empty
		{
			low = sentinelCycle.back();
			assert(low_array[low] != Globals::BIG_INT);
			sentinelCycle = list_sym_diff(sentinelCycle, redBoundary[low_array[low]]);

			if (low_array[low] >= death)
			{
				int currentColumnIdx = mapColorColumnIdx.at(low_array[low]);
				int deathColumnIdx = mapColorColumnIdx.at(death);
				annotation.set(currentColumnIdx - deathColumnIdx); // set the corresponding bit to be 1 
			}
		}

		// finally, construct the map which associates the edge with its annotation
		thread_mutex.lock();
		resEdgeAnnotations.insert({ std::make_pair(edge.first, edge.second), annotation });
		thread_mutex.unlock();
	}
}


/********************************************************************
* Description:	For a certain birth time, compute all the annotations of the sentinel
						edges. The computed annotations are stored in a std::map structure.
						We use bit vector to operate the annotation.
* Parameters:
* - redBoundary:			reduced boundary matrix
* - edgeMap:					a map, mapping two endpoints to an edge
* - low_array:				an array storing the pivot information
* - cell2v_list:				a converter which projects cells to their corresponding constituent vertices
* - death:						the given death time, which is a number from filtrationOrder
* - vertexNum:				the number of vertices in the whole topological space
* - resEdgeAnnotations:	the result edge annotations
********************************************************************/
void computeAnnotations(const vector<MatrixListType> & redBoundary, const map<pair<int, int>, int> & edgeMap, const vector<int> & low_array,
	const vector<MatrixListType> & cell2v_list, int death, int vertexNum, map<pair<int, int>, BitSet> & resEdgeAnnotations)
{
	resEdgeAnnotations.clear(); // clear up old data

	int low = redBoundary[death].back();
	MatrixListType edge, pivot = cell2v_list[low];

	// build a graph, which consists of all edges above birth (inclusive)
	adjacency_list_t graph(vertexNum);
	vector<pair<int, int>> graphEdges;
	for (size_t i = 0; i <= low; i++)
	{
		edge = cell2v_list[i];
		graph[edge[0]].push_back(neighbor(edge[1]));
		graph[edge[1]].push_back(neighbor(edge[0]));

		if (edge[0] > edge[1]) // sort the endpoints into nondecreasing order
			graphEdges.push_back(std::make_pair(edge[1], edge[0]));
		else
			graphEdges.push_back(std::make_pair(edge[0], edge[1]));
	}

	// compute the spanning tree
	adjacency_list_t spanningTree;
	vector<pair<int, int>> sentinelEdges;
	computeSpanningTree(graph, graphEdges, spanningTree, sentinelEdges);

	// next, for each sentinel edge, find its unique sentinel cycle in the spanning tree.
	// Afterwards, compute the annotation of this sentinel edge
	MatrixListType sentinelCycle;
	map<int, int> mapColorColumnIdx;
	int bettiNum = computeBettiNumber(redBoundary, redBoundary[death].back(), death, mapColorColumnIdx);
	BitSet annotation(bettiNum); // annotation, organized with bit vector

	// dispatch tasks to workers
	int num_workers = Globals::num_threads;
	vector<vector<pair<int, int>>> batch_sentinelEdges;
	int batch_size = sentinelEdges.size() / num_workers;
	for (int i = 0; i < num_workers - 1; i++)
	{
		batch_sentinelEdges.push_back(vector<pair<int, int>>());
		std::copy(sentinelEdges.begin() + i * batch_size, sentinelEdges.begin() + (i + 1)*batch_size, std::back_inserter(batch_sentinelEdges[i]));
	}
	batch_sentinelEdges.push_back(vector<pair<int, int>>());
	std::copy(sentinelEdges.begin() + (num_workers - 1) * batch_size, sentinelEdges.end(), std::back_inserter(batch_sentinelEdges[num_workers - 1]));

	// creat workers
	std::vector<std::thread> threadList;
	for (int i = 0; i < num_workers; i++)
	{
		threadList.push_back(std::thread(threadComputeAnnotation, ref(batch_sentinelEdges[i]),
			ref(spanningTree), ref(edgeMap), ref(low_array), death, bettiNum, ref(redBoundary), ref(mapColorColumnIdx),
			ref(resEdgeAnnotations)));
	}

	// wait workers to finish
	std::for_each(threadList.begin(), threadList.end(), std::mem_fn(&std::thread::join));
}


/********************************************************************
* Description:	Given a cycle, this function computes its annotation
* Parameters:
* - inputCycle:				the input cycle
* - cell2v_list:				a converter which projects cells to their corresponding constituent vertices
* - edgeAnnotations:		the annotations for all sentinel edges
* - resAnnotation:			the result cycle annotation
********************************************************************/
void computeCycleAnnotation(const MatrixListType & inputCycle, const vector<MatrixListType> & cell2v_list,
	const map<pair<int, int>, BitSet> & edgeAnnotations, BitSet & resAnnotation)
{
    resAnnotation.reset(); // reinitialization

	int ptr_1, ptr_2;
	MatrixListType edge;
	pair<int, int> key; // for query of edgeAnnotations
	map<pair<int, int>, BitSet>::const_iterator it;
	for (const auto & edgeIdx : inputCycle)
	{
	    edge = cell2v_list[edgeIdx];
		ptr_1 = edge[0];
		ptr_2 = edge[1];
		if (ptr_1 > ptr_2)
			SWAP(ptr_1, ptr_2);

		key.first = ptr_1;
		key.second = ptr_2;

		it = edgeAnnotations.find(key);
		if (it != edgeAnnotations.end()) // it is a sentinel edge, and we find it
		{
			resAnnotation ^= it->second;
		}
	}
}


/********************************************************************
* Description:	Given edge annotations and entry index, construct the 1d-version convering graph.
						The input original graph is specified by cell2v_list and low index, which together
						define the edge set of original graph.
* Parameters:
* - edgeAnnotations:		the input edge annotations
* - cell2v_list:				a converter which projects cells to their corresponding constituent vertices
* - low:							the low index of a homology class
* - entryIdx:					the entry index
* - vertexNum:				the number of vertices
* - resCoveringGraph:	the result 1d covering graph
********************************************************************/
void constructCoveringGraph1D(const map<pair<int, int>, BitSet> & edgeAnnotations, const vector<MatrixListType> & cell2v_list,
	int low, int entryIdx, int vertexNum, adjacency_list_t & resCoveringGraph)
{
	resCoveringGraph.clear();
	resCoveringGraph.resize(2*vertexNum); // for 1d covering graph, we have two copies of original graph

	MatrixListType edge;
	pair<int, int> key;
	int ptr_1, ptr_2;
	map<pair<int, int>, BitSet>::const_iterator it_annotation;
	
	for (size_t i = 0; i < low; i++) // we do not need to consider the 'low' edge
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
		if (it_annotation == edgeAnnotations.end() || (it_annotation->second).checkBit(entryIdx) == false) // it is not a sentinel edge, or the bit at the given entry is 0
		{
			// first copy
			resCoveringGraph[edge[0]].push_back(neighbor(edge[1]));
			resCoveringGraph[edge[1]].push_back(neighbor(edge[0]));

			// second copy
			resCoveringGraph[edge[0] + vertexNum].push_back(neighbor(edge[1] + vertexNum));
			resCoveringGraph[edge[1] + vertexNum].push_back(neighbor(edge[0] + vertexNum));
		}
		else
		{
			// first copy
			resCoveringGraph[edge[0]].push_back(neighbor(edge[1] + vertexNum));
			resCoveringGraph[edge[1] + vertexNum].push_back(neighbor(edge[0]));

			// second copy
			resCoveringGraph[edge[0] + vertexNum].push_back(neighbor(edge[1]));
			resCoveringGraph[edge[1]].push_back(neighbor(edge[0] + vertexNum));
		}
	}
}


/********************************************************************
* Description:	Apply backward BFS on unweighted graph to find the length of shortest path
						between two vertices. During this process, we also record the distance from
						other vertices to the target, which will be used for building heuristic database.
* Parameters:
* - source:						the 'target' vertex in original graph
* - target:						the 'source' vertex in original graph
* - coveringGraph:			the input covering graph
* - heuristic_database:	the heuristic database
********************************************************************/
int BFS_Shortest_Path_Length(int source, int target, const adjacency_list_t & coveringGraph,
	unordered_map<int, int> & heuristic_database)
{
	// first, check if the distance from source to target has been computed before
	unordered_map<int, int>::const_iterator hash_it;
	hash_it = heuristic_database.find(target);
	if (hash_it != heuristic_database.end())
		return hash_it->second;

	// otherwise, we need to do BFS to find the path length
	int vertexNum = coveringGraph.size();

	unordered_set<int> hasVisited;
	queue<int> searchQ;

	searchQ.push(source); // backward BFS
	hasVisited.insert(source);
	heuristic_database[source] = 0;

	int currNode, neighborNode;
	int pathLen = -1;

	while (true)
	{
		currNode = searchQ.front();
		searchQ.pop();
		pathLen = heuristic_database.at(currNode);

		if (currNode == target)
			return pathLen;

		// visit its neighbors
		const std::vector<neighbor> & neighbors = coveringGraph[currNode];
		for (const auto & nb : neighbors)
		{
			if (hasVisited.find(nb.target) != hasVisited.end())
				continue;

			hasVisited.insert(nb.target);
			searchQ.push(nb.target);
			heuristic_database[nb.target] = pathLen + 1;
		}
	}
}


/********************************************************************
* Description:	Given source and target, compute the length of shortest path between them.
						This shortest path should have the same annotation as target annotation.
* Parameters:
* - coveringGraphs:		a vector containing different 1d covering graphs
* - targetAnnotation:		target annotation from source to target
* - source:						the given source
* - target:						the given target
* - BettiNum:					the Betti number
* - vertexNum:				the number of vertices
* - heuristic_database:	heuristic database
********************************************************************/
int computeHeuristic(const vector<adjacency_list_t> & coveringGraphs, const BitSet & targetAnnotation,
	int source, int target, int BettiNum, int vertexNum, vector<unordered_map<int, int>> & heuristic_database)
{
	assert(coveringGraphs.size() == BettiNum);

	int currH, maxH = -1;
	int u = target, v; // backward search from target to source

	for (int i = 0; i < BettiNum; i++)
	{
		if (targetAnnotation.checkBit(i) == false) // compute the path with even number of edges from Ei
			v = source;
		else
			v = source + vertexNum;

		currH = BFS_Shortest_Path_Length(u, v, coveringGraphs[i], heuristic_database[i]);
		if (maxH < currH)
			maxH = currH;
	}

	return maxH;
}

#endif
