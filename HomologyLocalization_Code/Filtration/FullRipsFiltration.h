#ifndef FULL_RIPS_FILTRATION_H
#define FULL_RIPS_FILTRATION_H

#include "AbstractFiltration.h"
#include "InputFileInfo.h"

template<int maxDim>
class FullRipsFiltration : public AbstractFiltration<maxDim, 2, 2> // arrayDim = 2 for 2D distance matrix; vertexDim = 2 for 2D point index
{
	typedef blitz::TinyVector<int, 2> matrixIdx;

public:
	// constructor
	FullRipsFiltration(const blitz::Array<double, 2> *const p, const InputFileInfo &info) : mDistanceMatrix(p)
	{
		mPointsNum = p->extent(1);

		precomputeBinomials();
		precomputeCellNums();
	}

	// initialize the vertex list, which will be used for indexing birth time and death time
	void init(std::vector<matrixIdx> *vList)
	{
		if (!vList->empty())
			return;

		matrixIdx idx;
		for (int i = 0; i < mPointsNum; ++i)
		{
			for (int j = i ; j < mPointsNum;++j)
			{
				idx = i, j;
				vList->push_back(idx);
			}
		}
	}

	// return the number of cells in dimension d
	int getSizeInDim(int d)
	{
		return mBreakPoints[d + 1] - mBreakPoints[d];
	}

	// initialize the birth list, which records the birth time for each cell;
	// and the cellToVertex list, which records the mapping from each cell to the component vertices
	void initList(std::vector<int> *birth_list, vector<MatrixListType> *cell2v_list, int d)
	{
		typedef std::tuple<int, double, int, vector<int>> elemType;

		int cellNum = getSizeInDim(d);
		birth_list->resize(cellNum);
		birth_list->clear();

		vector<int> pointsIdxVector;
		vector<elemType> birthVector;
		double diameter;
		int matrixIntIdx;
		std::pair<int, int> idxPair;

		// compute the birth time index (integer) for each complex cell
		for (int i = mBreakPoints[d]; i < mBreakPoints[d+1]; ++i)
		{
			conversion(i, pointsIdxVector);
			diameter = getDiameter(pointsIdxVector.begin(), pointsIdxVector.end(), idxPair);
			matrixIntIdx = computeIntIndex(idxPair);
			std::sort(pointsIdxVector.begin(), pointsIdxVector.end()); // need to sort it in ascending order for future use

			birthVector.push_back(std::make_tuple(i, diameter, matrixIntIdx, pointsIdxVector));
		}

		// sort the cells in nondecreasing order
		std::sort(birthVector.begin(), birthVector.end(),
			[](elemType const & a, elemType const & b) {return std::get<1>(a) < std::get<1>(b);});

		// copy the birth time
		std::transform(birthVector.begin(), birthVector.end(), std::back_inserter(*birth_list),
			[](elemType const & a) {return std::get<2>(a);});

		// compute the cell_ to_vertices list
		std::transform(birthVector.begin(), birthVector.end(), std::back_inserter(*cell2v_list),
			[](elemType const & a) {return std::get<3>(a);});
	}

	// compute the d-dimensional boundary matrix 
	void calculateBoundaries(vector<MatrixListType> *boundary, int d, const vector<bool> & will_be_cleared)
	{
		typedef std::tuple<int, double, vector<int>> elemType;

		int cellNum = getSizeInDim(d);
		boundary->resize(cellNum);
		boundary->clear();

		vector<int> pointsIdxVector, subComplexIdxVector;
		vector<elemType> birthVector;
		vector<pair<int, double>> subcomplexVector;
		double diameter;
		int subComplexIdx;
		std::pair<int, int> idxPair;

		// first, list the subcomplex of dimension d-1
		for (int i = mBreakPoints[d-1]; i < mBreakPoints[d]; ++i)
		{
			conversion(i, pointsIdxVector);
			diameter = getDiameter(pointsIdxVector.begin(), pointsIdxVector.end());
			subcomplexVector.push_back(std::make_pair(i, diameter));
		}

		// sort
		std::sort(subcomplexVector.begin(), subcomplexVector.end(),
			[](pair<int, double> const & a, pair<int, double> const &b) {return a.second < b.second;});

		// create mapping
		std::map<int, int> map2order;
		for (int i = 0; i< subcomplexVector.size(); ++i)
		{
			map2order[subcomplexVector[i].first] = i;
		}
		

		for (int i = mBreakPoints[d]; i < mBreakPoints[d + 1]; ++i)
		{
			conversion(i, pointsIdxVector);
			diameter = getDiameter(pointsIdxVector.begin(), pointsIdxVector.end(), idxPair);
			subComplexIdxVector.clear();

			for (vector<int>::const_iterator skip = pointsIdxVector.begin(); skip != pointsIdxVector.end(); skip++)
			{
				subComplexIdx = conversion_with_skip<vector<int>::const_iterator>(pointsIdxVector.begin(), pointsIdxVector.end(), skip);
				subComplexIdxVector.push_back(map2order[subComplexIdx]);
			}

			birthVector.push_back(std::make_tuple(i, diameter, subComplexIdxVector));
		}

		std::sort(birthVector.begin(), birthVector.end(),
			[](elemType const & a, elemType const & b) {return std::get<1>(a) < std::get<1>(b);});

		std::transform(birthVector.begin(), birthVector.end(), std::back_inserter(*boundary),
			[](elemType const & a) {return std::get<2>(a);});

		for (std::vector<MatrixListType>::iterator it = boundary->begin(), end = boundary->end(); it != end; ++it)
			mysort(*it);
	}

private:
	// precompute the number of cells for each dimension
	void precomputeCellNums()
	{
		mBreakPoints.push_back(0);
		for (int i = 1; i <= maxDim + 1; ++i)
		{
			mBreakPoints.push_back(mBreakPoints[i - 1] + mBinomials[mPointsNum][i]);
		}
	}

	// precompute the binomials
	void precomputeBinomials()
	{
		for (int i = 0; i <= mPointsNum; ++i)
		{
			vector<int> binomRow;
			for (int j = 0; j <= maxDim + 1; ++j)
			{
				binomRow.push_back(binom(i, j)); // compute C(i, j)
			}
			mBinomials.push_back(binomRow);
		}
	}

	// compute binomial C(n, k)
	int binom(int n, int k) const
	{
		if (n < k)
			return 0;
		
		int result = 1;
		for (int i = n - k + 1; i <= n; ++i) 
		{
			result *= i;
		}

		for (int i = k; i >= 1; --i) 
		{
			result /= i;
		}

		return result;
	}

	// given an index, return its corresponding dimension
	int getLocalDim(int idx) const
	{
		int k = 0;
		while (k < maxDim + 1 && mBreakPoints[k + 1] <= idx) 
		{
			k++;
		}

		return k;
	}

	// index -> set of points (in decreasing index order)
	void conversion(int idx, vector<int> & out)
	{
		out.clear();

		int k = getLocalDim(idx);
		idx -= mBreakPoints[k];
		int bcoeff = mPointsNum - 1;

		for (; k >= 0; k--) 
		{
			while (mBinomials[bcoeff][k + 1] > idx)
			{
				bcoeff--;
			}
			out.push_back(bcoeff);
			idx -= mBinomials[bcoeff][k + 1];
			bcoeff--;
		}
	}

	// given a complex, compute its diameter. Besides, return the index into the distance matrix
	template<typename InputIterator>
	double getDiameter(InputIterator begin, InputIterator end, std::pair<int, int> & index) const
	{
		if (begin == end) 
		{
			index.first = 0;
			index.second = 0;
			return 0.;
		}

		double max = 0.;
		InputIterator curr = begin;

		do 
		{
			for (InputIterator run = curr + 1; run != end; run++) 
			{
				double cdist = (*mDistanceMatrix)(*curr, *run);
				if (cdist > max)
				{
					max = cdist;
					index.first = *curr;
					index.second = *run;
				}
			}
			curr++;
		} while (curr != end);

		return max;
	}

	// overload function
	template<typename InputIterator>
	double getDiameter(InputIterator begin, InputIterator end) const
	{
		if (begin == end)
		{
			return 0.;
		}

		double max = 0.;
		InputIterator curr = begin;

		do
		{
			for (InputIterator run = curr + 1; run != end; run++)
			{
				double cdist = (*mDistanceMatrix)(*curr, *run);
				if (cdist > max)
				{
					max = cdist;
				}
			}
			curr++;
		} while (curr != end);

		return max;
	}

	// convert the given complex into one of its subcomplexes, and return the corresponding index
	// Assume that the input sequence is sorted in decreasing order!
	template<typename InputIterator>
	int conversion_with_skip(InputIterator begin, InputIterator end, InputIterator skip) const
	{
		int dist = std::distance(begin, end);
		if (skip != end)
		{
			dist--;
		}
		assert(dist <= (int)maxDim);

		InputIterator it = begin;
		int ind = 0;

		for (int k = dist; k >= 1; k--) 
		{
			if (it == skip) 
			{
				it++;
			}
			ind += mBinomials[*it++][k];
		}
		ind += mBreakPoints[dist - 1];

		return ind;
	}

	// given an integer index, convert it to the pair index into the distance matrix
	std::pair<int, int> computePairIndex(int idx) const
	{
		int row = -1, col = -1;
		int rowNum = mPointsNum;

		while (idx >= 0)
		{
			row++;
			idx -= rowNum;
			rowNum--;
		}

		col = mPointsNum + idx;

		return std::make_pair(row, col);
	}

	// given a pair of indices, covert it to the integer index
	int computeIntIndex(const std::pair<int, int> & pairIdx) const
	{
		int resIdx = 0;
		int pointsPerRow = mPointsNum;

		int rowCnt = pairIdx.first;
		int colCnt = pairIdx.second;
		if (rowCnt > colCnt)
			SWAP(rowCnt, colCnt);

		while (rowCnt > 0)
		{
			resIdx += pointsPerRow;
			rowCnt--;
			pointsPerRow--;
		}

		resIdx += pointsPerRow - (mPointsNum - colCnt);

		return resIdx;
	}

private:
	vector<int> mBreakPoints; // cumulative cell numbers

	vector<vector<int>> mBinomials; // store the precomputed binomials

	const blitz::Array<double, 2> * const mDistanceMatrix; // pointer to the distance matrix data

	int mPointsNum; // number of cloud points
};

#endif // !FULL_RIPS_FILTRATION_H

