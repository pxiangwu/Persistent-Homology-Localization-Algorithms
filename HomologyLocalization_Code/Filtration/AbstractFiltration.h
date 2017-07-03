#ifndef ABSTRACT_FILTRATION_H
#define ABSTRACT_FILTRATION_H


template<int dim, int arrayDim, int vertexDim>
class AbstractFiltration
{
public:
	// initialize the vertex list, which will be used for indexing birth time and death time
	virtual void init(std::vector<blitz::TinyVector<int, vertexDim>> *vList) = 0; 

	// return the number of cells in dimension d
	virtual int getSizeInDim(int d) = 0;

	// initialize the birth list, which records the birth time for each cell;
	// and the cellToVertex list, which records the mapping from each cell to the component vertices
	virtual void initList(std::vector<int> *birth_list, vector<MatrixListType> *cell2v_list, int d) = 0;

	// compute the d-dimensional boundary matrix 
	virtual void calculateBoundaries(vector<MatrixListType> *boundary, int d, const vector<bool> & will_be_cleared) = 0;
};

#endif // !ABSTRACT_FILTRATION_H

