#ifndef INCLUDED_PERSISTENCE_CALCULATOR_H
#define INCLUDED_PERSISTENCE_CALCULATOR_H

// Trick:
// Rather than operating on an existing matrix, we record the columns
// to be cleared, and clear them once we have calculated the matrix.
// We need to store only one matrix at a time.
#include <ctime>
#include "Algorithms/Reduction.h"
#include "Algorithms/AnnotatingEdges.h"
#include "Algorithms/AStar.h"
#include "BitSet.h"
#include "Globals.h"
#include "Algorithms/ExhaustiveSearch.h"
#include "Filtration/FullRipsFiltration.h"
#include "Filtration/CubicalFiltration.h"
#include "Filtration/SimComplexFiltration.h"

vector<vector<int>> birth_lists;
vector<vector<MatrixListType>> cell2v_lists;


// By switching FiltrationGeneratorType it should be possible to use for example simplicial complexes
template<int dim, int arrayDim = dim, int vertexDim = dim, typename FiltrationGeneratorType = CubicalFiltration<dim>, int type = 0 >
struct PersistenceCalculator 
{
	typedef blitz::TinyVector<int, vertexDim> Vertex;
	typedef vector<PersPair<Vertex> > PersResultContainer;

	// -- Save persistence infomation
	template<typename NDArray>
	void SavePersistence(NDArray * phi, const vector<Vertex> &vList, vector< int > & lowerCellList,	vector< int > & upperCellList, 
		vector< int > & low_array, const double pers_thd, PersResultContainer &veList,
		/* for reduction list*/
		vector< MatrixListType > & red_list, vector< MatrixListType > & red_cell2v_list, 	vector< MatrixListType > & final_red_list,
		vector< MatrixListType > & bd_list, 	vector< MatrixListType > & bd_cell2v_list, vector< MatrixListType > & final_boundary_list)
	{
		assert(final_red_list.empty());
		assert(final_boundary_list.empty());

		// output vertex-edge pairs whose persistence is bigger than pers_thd
		for (int i = 0; i < lowerCellList.size(); i++)
		{
			int tmp_int = low_array[i];
			if (tmp_int == Globals::BIG_INT) 
				continue;

			int vBirth = lowerCellList[i];
			int vDeath = upperCellList[tmp_int];

			double tmp_death = (*phi)(vList[vDeath]);
			double tmp_birth = (*phi)(vList[vBirth]);
			double tmp_pers = tmp_death - tmp_birth;

			assert(tmp_pers >= 0);

			if (tmp_pers > pers_thd) 
			{
				veList.push_back(PersPair<Vertex>(vList[vBirth],	vList[vDeath], tmp_pers, tmp_birth, tmp_death));		

				// save the reduction lists
				MatrixListType tmp_list;
				assert(!red_list[tmp_int].empty());
				for (MatrixListType::iterator tmpiter = red_list[tmp_int].begin(); tmpiter != red_list[tmp_int].end(); tmpiter++) 
				{
					tmp_list = list_union(tmp_list, red_cell2v_list[*tmpiter]);
				}
				final_red_list.push_back(tmp_list);

				// save the boundary lists
				MatrixListType tmp_boundary_list;
				assert(!bd_list[tmp_int].empty());
				for (MatrixListType::iterator tmpiter = bd_list[tmp_int].begin(); tmpiter != bd_list[tmp_int].end(); tmpiter++) 
				{
					tmp_boundary_list = list_union(tmp_boundary_list, bd_cell2v_list[*tmpiter]);
				}
				final_boundary_list.push_back(tmp_boundary_list);
			}
		}
	}


	// -- Perform persistence homology algorithm
	void calcPersistence(blitz::Array<double, arrayDim> *phi, const double pers_thd,
		vector<PersResultContainer> &result_lists, vector<Vertex> & _vList, const InputFileInfo &info)
	{
		time_t wholestart, wholeend, redstart, redend;
		double wholetime = 0, redtime = 0;

		time(&wholestart);

		vector<Vertex> *vList = &_vList;

		birth_lists.assign(dim + 1, vector<int>());
		cell2v_lists.assign(dim + 1, vector<MatrixListType>());

		int sizes[dim + 1] = { 0 };


		// initialize vertex lists, birth_lists, and cell2v_lists
		FiltrationGeneratorType filtration(phi, info);
		filtration.init(vList);

		assert(!vList->empty());
		for (int i = 0; i <= dim; i++)
		{
			sizes[i] = filtration.getSizeInDim(i);
			filtration.initList(&birth_lists[i], &cell2v_lists[i], i);
		}


		vector<vector<int>> low_arrays(dim + 1); // low_array[i] refers to the cell with "pivot cell" i in the boundary matrix
		for (int i = 1; i <= dim; i++)
			low_arrays[i].assign(sizes[i - 1], Globals::BIG_INT);

		vector<vector<MatrixListType>> boundaries(dim + 1); // boundary matrices
		vector<bool> willBeCleared(sizes[dim], false);
		

		// save for each negative simplex the simplices used to reduce it
		vector< MatrixListType > reduction_list;
		vector< MatrixListType > final_reduction_list;
		vector< MatrixListType > final_boundary_list;


		for (int d = dim; d >= 1; d--)
		{  
			filtration.calculateBoundaries(&boundaries[d], d, willBeCleared);

			willBeCleared.assign(sizes[d - 1], false);
			time(&redstart);

			// perform boundary matrix reduction
			reduceND(willBeCleared, birth_lists[d], boundaries[d], low_arrays[d], reduction_list);

			// for 1D homology, employ optimal shortest cycle algorihm to further reduce the boundary matrix
			if (d == 2 && Globals::use_optimal_alg == true)
			{
				std::map<std::pair<int, int>, int> edgeMap;
				constructMap_Edge2Ptr(edgeMap, cell2v_lists[d - 1]);

				switch (Globals::which_alg)
				{
				case Globals::Algorithm::HEURISTIC_BASED_ALG:
					reduceND_AStar<arrayDim, vertexDim>(phi, *vList, birth_lists[d - 1], birth_lists[d], boundaries[d], edgeMap, cell2v_lists[d - 1], sizes[0], low_arrays[d]);
					break;
				case Globals::Algorithm::CLASSICAL_ALG:
					reduceND_ExhaustiveSearch<arrayDim, vertexDim>(phi, *vList, birth_lists[d - 1], birth_lists[d], boundaries[d], edgeMap, cell2v_lists[d - 1], sizes[0], low_arrays[d]);
					break;
				default:
					break;
				}
			}// end if

			time(&redend);
			redtime += difftime(redend, redstart);
			cout << "Reduced dimension " << d << endl;

			// save persistence, boundaries, red_list for this dimension, such that the memory could be cleaned
			SavePersistence(phi, *vList, birth_lists[d - 1], birth_lists[d], low_arrays[d], pers_thd, result_lists[d - 1],
				reduction_list, cell2v_lists[d], final_reduction_list, boundaries[d], cell2v_lists[d - 1], final_boundary_list);
			
			// release memory
			cell2v_lists[d].clear();

			// save reduction and boundary files
			BinaryPersistentPairsSaver<dim, arrayDim, vertexDim> binSaver;
			stringstream output_red_file;
			output_red_file << info.input_path;
			output_red_file << ".red." << d;
			binSaver.saveOneDimReduction(final_reduction_list, (*vList), output_red_file.str().c_str(), type);

			stringstream output_boundary_file;
			output_boundary_file << info.input_path;
			output_boundary_file << ".bnd." << d;
			binSaver.saveOneDimReduction(final_boundary_list, (*vList), output_boundary_file.str().c_str(), type);

			cout << "Saved dimension " << d << endl << endl;

			reduction_list.clear();
			final_reduction_list.clear();
			boundaries[d].clear();
			final_boundary_list.clear();
		}// end for

		OUTPUT_MSG("Reduction done");

		time(&wholeend);
		wholetime = difftime(wholeend, wholestart);
		cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Filtration building time = " << setprecision(3) << (wholetime - redtime) / 60.0 << " Min" << endl;
		cout << "Reduction  time \t = " << setprecision(3) << redtime / 60.0 << " Min" << endl;
		cout << "+++++++++++++++++++++++++++++++++++++++++" << endl << endl;

		OUTPUT_MSG("Recording persistence pairs")
		OUTPUT_MSG("Finished");
	}

};

#endif
