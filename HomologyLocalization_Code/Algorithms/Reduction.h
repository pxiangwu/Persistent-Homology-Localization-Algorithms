#ifndef INCLUDED_REDUCTION_H
#define INCLUDED_REDUCTION_H

#include <stdexcept>
#include "Globals.h"

/*********************************************************************
* This function reduces a boundary matrix represented by its 'low_array'.
* 'continue_reduction' parameter chooses which mode to be used for reduction:
*		0 -- no continuous reduction
*		1 -- continuous reduction in a conservative manner, recommended
*		2 -- continuous reduction aggressively (not recommended)
*********************************************************************/
void reduceND(vector<bool> &willBeCleared, vector<CellNrType> &upperList, vector<MatrixListType> &boundary_upper, 
	vector<int> &low_array, vector< MatrixListType > & reduction_list, int continue_reduction = 0)
{
	OUTPUT_MSG("Reducing cells, total number = " << upperList.size());

	int total_reduction_effect = 0;

	for (size_t i = 0, sz = upperList.size(); i < sz; i++) 
	{
		reduction_list.push_back(MatrixListType());
		MY_ASSERT(reduction_list.size() == i + 1);

		if (boundary_upper[i].empty())
			continue;

		// initialize the reduction list first
		reduction_list[i].push_back(i);

		int low = boundary_upper[i].back();
		int column_used = 0;
		while (!boundary_upper[i].empty() && low_array[low] != Globals::BIG_INT)
		{
			assert(low_array[low] < i);
			assert(low == boundary_upper[low_array[low]].back());
			assert(!boundary_upper[i].empty());
			assert(!boundary_upper[low_array[low]].empty());

			boundary_upper[i] = list_sym_diff(boundary_upper[i], boundary_upper[low_array[low]]);

			// update the reduction list as well
			reduction_list[i] = list_sym_diff(reduction_list[i], reduction_list[low_array[low]]);
			if (!boundary_upper[i].empty()) 
			{
				int old_low = low;
				low = boundary_upper[i].back();
				assert(low<old_low);
			}

			column_used++;
		}

		if (!boundary_upper[i].empty()) 
		{
			assert(low >= 0);
			assert(low_array[low] == Globals::BIG_INT);
			low_array[low] = i;

			willBeCleared[low] = true;
		}


		// further reduction based on heuristics
		int before_reduction_size = boundary_upper[i].size();
		if (continue_reduction == 0) continue;

		MatrixListType tmp_bdry_v;
		int new_pivot = low;
		MatrixListType::iterator new_pivot_iter = lower_bound(boundary_upper[i].begin(),	boundary_upper[i].end(), new_pivot);
		assert(new_pivot_iter == boundary_upper[i].end() - 1);

		while (new_pivot_iter != boundary_upper[i].begin())
		{
			new_pivot_iter--;
			new_pivot = *new_pivot_iter;

			if (low_array[new_pivot] != Globals::BIG_INT)
			{
				// the new_pivot is paired with some column before ith column, use it for reduction
				assert(low_array[new_pivot] < i);
				if (continue_reduction == 1) 
				{
					tmp_bdry_v = list_sym_diff(boundary_upper[i], boundary_upper[low_array[new_pivot]]);
					if (tmp_bdry_v.size() < boundary_upper[i].size()) 
					{
						boundary_upper[i] = tmp_bdry_v;
						reduction_list[i] = list_sym_diff(reduction_list[i], reduction_list[low_array[new_pivot]]);
					}
				}
				else if (continue_reduction == 2) 
				{
					boundary_upper[i] = list_sym_diff(boundary_upper[i], boundary_upper[low_array[new_pivot]]);
					reduction_list[i] = list_sym_diff(reduction_list[i], reduction_list[low_array[new_pivot]]);
				}
			}

			new_pivot_iter = lower_bound(boundary_upper[i].begin(), boundary_upper[i].end(), new_pivot);
		}

		int after_reduction_size = boundary_upper[i].size();
		total_reduction_effect += after_reduction_size - before_reduction_size;
	}

	if (continue_reduction != 0) 
	{
		cout << "Total additional reduction = " << total_reduction_effect << endl << endl;
	}
}

#endif
