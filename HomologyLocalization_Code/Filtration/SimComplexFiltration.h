#ifndef SIM_COMPLEX_FILTRATION_H
#define SIM_COMPLEX_FILTRATION_H

#include "AbstractFiltration.h"
#include "InputFileInfo.h"

template<int dim>
class SimComplexFiltration : public AbstractFiltration<dim, 1, 1> // arrayDim = 1 for storing point values in 1D array; vertexDim = 1 for 1D index 
{
	typedef blitz::TinyVector<int, 1> valueIdx;

	template<typename ArrayType>
	struct PhiComparator // used in sort
	{
		const ArrayType &phi;
		PhiComparator(const ArrayType * const p) : phi(*p) {}

		template<typename VertexType>
		bool operator()(const VertexType &a, const VertexType &b) const
		{
			return phi(a) < phi(b);
		}
	};

public:
	// constructor
	SimComplexFiltration(const blitz::Array<double, 1> *const p, const InputFileInfo &info) : pointsVal(p), mFileInfo(info) 
	{
		readData();
	}

	// initialize the vertex list, which will be used for indexing birth time and death time
	void init(std::vector<valueIdx> *vList)
	{
		if (!vList->empty())
			return;

		vList->resize(pointsVal->extent(0));
		std::iota(vList->begin(), vList->end(), 0);

		// sort the list in increasing order
		//std::sort(vList->begin(), vList->end(), PhiComparator<blitz::Array<double, 1>>(this->pointsVal));
	}

	// return the number of cells in dimension d
	int getSizeInDim(int d)
	{
		return mCellNums[d];
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
		double filterValue;
		int maxIdx;

		for (int i = 0; i < cellNum; ++i)
		{
			filterValue = getFilterValue(mCells[d][i].begin(), mCells[d][i].end(), maxIdx);

			pointsIdxVector.clear();
			std::copy(mCells[d][i].begin(), mCells[d][i].end(), std::back_inserter(pointsIdxVector));
			std::sort(pointsIdxVector.begin(), pointsIdxVector.end()); // need to sort it in ascending order for future use
			
			birthVector.push_back(std::make_tuple(i, filterValue, maxIdx, pointsIdxVector));
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
		typedef std::pair<double, vector<int>> filterIdxType;
		typedef std::tuple<int, double, vector<int>> elemType;

		int cellNum = getSizeInDim(d);
		boundary->resize(cellNum);
		boundary->clear();

		vector<int> pointsIdxVector, subcomplexIdxVector, subcomplex;
		vector<filterIdxType> subcomplexVector;
		vector<elemType> birthVector;
		double filterValue;
		int maxIdx;

		// first, list the subcomplex of dimension d-1
		int subComplexNum = getSizeInDim(d - 1);
		for (int i = 0; i < subComplexNum; ++i)
		{
			pointsIdxVector.clear();
			std::copy(mCells[d-1][i].begin(), mCells[d-1][i].end(), std::back_inserter(pointsIdxVector));
			std::sort(pointsIdxVector.begin(), pointsIdxVector.end()); // sort it such that it could be used as a key in map

			filterValue = getFilterValue(mCells[d-1][i].begin(), mCells[d-1][i].end(), maxIdx);

			subcomplexVector.push_back(std::make_pair(filterValue, pointsIdxVector));
		}

		// sort
		std::sort(subcomplexVector.begin(), subcomplexVector.end(),
			[](filterIdxType const & a, filterIdxType const &b) {return a.first < b.first;});

		// create mapping
		std::map<vector<int>, int> map2order;
		for (int i = 0; i < subcomplexVector.size(); ++i)
		{
			map2order[subcomplexVector[i].second] = i;
		}


		// create the boundary matrix
		for (int i = 0; i < cellNum; ++i)
		{
			pointsIdxVector.clear();
			std::copy(mCells[d][i].begin(), mCells[d][i].end(), std::back_inserter(pointsIdxVector));
			std::sort(pointsIdxVector.begin(), pointsIdxVector.end());

			filterValue = getFilterValue(mCells[d][i].begin(), mCells[d][i].end(), maxIdx);
			
			subcomplexIdxVector.clear();
			for (auto skip : pointsIdxVector)
			{
				subcomplex.clear();
				std::copy_if(pointsIdxVector.begin(), pointsIdxVector.end(), std::back_inserter(subcomplex),
					[skip](const int & t) {return t != skip;});
				subcomplexIdxVector.push_back(map2order.at(subcomplex)); // get the integer index for this subcomlex, and store it
			}

			birthVector.push_back(std::make_tuple(i, filterValue, subcomplexIdxVector));
		}

		// sort the boundary matrix columns
		std::sort(birthVector.begin(), birthVector.end(),
			[](elemType const & a, elemType const & b) {return std::get<1>(a) < std::get<1>(b);});

		std::transform(birthVector.begin(), birthVector.end(), std::back_inserter(*boundary),
			[](elemType const & a) {return std::get<2>(a);});

		for (std::vector<MatrixListType>::iterator it = boundary->begin(), end = boundary->end(); it != end; ++it)
			mysort(*it);
	}

private:
	// read the input data file
	void readData()
	{
		std::ifstream f(mFileInfo.input_path.c_str(), std::ios::binary | std::ios::in);

		int fileType, maxDim, numPoints, dimPoints;
		f.read(reinterpret_cast<char*>(&fileType), sizeof(int));	 // read file type
		f.read(reinterpret_cast<char*>(&maxDim), sizeof(int));	 // read maximum dimension
		f.read(reinterpret_cast<char*>(&numPoints), sizeof(int)); // read number of points
		f.read(reinterpret_cast<char*>(&dimPoints), sizeof(int)); // read dimension of each point

		mCells.resize(maxDim + 1); // initialize for the cells of dimension 0
		for (int i = 0; i < numPoints; ++i)
		{
			mCells[0].push_back(vector<int>(1, i));
		}

		mCellNums.resize(maxDim + 1);
		mCellNums[0] = numPoints;

		vector<vector<double>> posPoints;
		posPoints.resize(numPoints, vector<double>(dimPoints)); // read point positions
		for (int i = 0; i < numPoints; ++i)
		{
			f.read(reinterpret_cast<char*>(posPoints[i].data()), sizeof(double)*dimPoints);
		}

		blitz::Array<double, 1> tempArr;
		tempArr.resize(numPoints);
		f.read(reinterpret_cast<char*>(tempArr.data()), sizeof(double)*tempArr.size());  // read point values

		int cellDim, cellNum;
		while (true)
		{
			f.read(reinterpret_cast<char*>(&cellDim), sizeof(int));	 // read cell dimension
			f.read(reinterpret_cast<char*>(&cellNum), sizeof(int)); // read number of cells in this dimension

			if (f.eof())
				break;

			mCellNums[cellDim] = cellNum;
			mCells[cellDim].resize(cellNum);
			for (int i = 0; i < cellNum; ++i)
			{
				mCells[cellDim][i].resize(cellDim + 1);
				for (int j = 0; j <= cellDim; ++j)
				{
					f.read(reinterpret_cast<char*>(&mCells[cellDim][i][j]), sizeof(int));	 // read cell dimension
				}
			}
		}

		f.close();
	}

	template<typename InputIterator>
	double getFilterValue(InputIterator begin, InputIterator end, int & maxIdx) const
	{
		InputIterator curr = begin;
		double currVal;

		double max = (*pointsVal)(*curr);
		maxIdx = *curr;
		curr++;

		while (curr != end)
		{
			currVal = (*pointsVal)(*curr);
			if (currVal > max)
			{
				max = currVal;
				maxIdx = *curr;
			}

			curr++;
		}

		return max;
	}

private:
	// The filter function
	const blitz::Array<double, 1> *const pointsVal;

	// The input file infomation
	const InputFileInfo mFileInfo;

	// Number of cells
	vector<int> mCellNums;

	// Cells (consisting of indices)
	vector<vector<vector<int>>> mCells;
};

#endif // !SIM_COMPLEX_FILTRATION_H

