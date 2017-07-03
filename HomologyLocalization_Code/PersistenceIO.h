#ifndef PERSISTENCEIO_H
#define PERSISTENCEIO_H

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <functional>
#include <vector>
#include <malloc.h>
#include <ctime>
#include "Debugging.h"

#define SWAP(a, b)  do { a ^= b; b ^= a; a ^= b; } while ( 0 )

typedef vector<int> MatrixListType;

/************************************************************************/
/* Save persistence values and reduction results and reduction results to binary file			 */
/************************************************************************/
void saveResults(const char *fileName,unsigned int *data, const vector<int> &header, int nrPairs);


template<int d>
struct TextPersistentPairsSaver
{
	template<typename ContT>
	void savePers(vector<ContT> &res, const char * output_fname)
	{
		OUTPUT_MSG( endl<< "---------- writing result ---------------" );

		cout << "writing text output to:  " << output_fname << endl;

		fstream output_filestr(output_fname, fstream::out | fstream::trunc);

		for (int i = 0; i < d; i++)
		{
			//sort(res[i].begin(), res[i].end());
			output_filestr << "[" << std::to_string(i) << "D-Cell" << ", " << std::to_string(i+1) << "D-Cell] Pairs = " << res[i].size() << endl;
			typename ContT::iterator veiter;
			for(veiter=res[i].begin(); veiter!=res[i].end(); veiter++)
				output_filestr << veiter->birth << "\t" << veiter->death << endl;
			output_filestr << endl;
		}	

		OUTPUT_MSG( endl<< "---------- writing result finished ---------------" );
	}
};


template<int d, int arrayDim = d, int vertexDim = d>
struct BinaryPersistentPairsSaver
{
	typedef blitz::TinyVector<int, vertexDim> Vertex;

	template<typename ContT>
	void saveOutput(vector<ContT> &res, vector< Vertex > & vList, const char * output_fname)
	{	
		assert(res.size() == d);	

		int count = 0;
		for (int i = 0; i < res.size(); i++){
			count += res[i].size();
		}

		for (int i = 0; i < d; i++)			
			cout << "[" << std::to_string(i) << "D-Cell" << ", " << std::to_string(i + 1) << "D-Cell] Pairs = " << res[i].size() << endl;

		unsigned int *pairArray = new unsigned int[count*2*vertexDim];

		int index = 0;

		vector<int> header(d);

		for (int i = 0; i < d; i++)
		{
			typename ContT::iterator veiter;
			for(veiter=res[i].begin(); veiter!=res[i].end(); veiter++)
			{
				Vertex birthcoord = veiter->birthV ;
				for( size_t ii = 0; ii < vertexDim; ii ++ )
					pairArray[index++] = birthcoord[vertexDim-1-ii] + 1;

				Vertex deathcoord = veiter->deathV;
				for( size_t ii = 0; ii < vertexDim; ii ++ )
					pairArray[index++] = deathcoord[vertexDim-1-ii] + 1;
			}
			header[i] = res[i].size();
		}	

		saveResults(output_fname, pairArray, header, index);

		delete[] pairArray;
	}

	void saveOneDimReduction(const vector< MatrixListType > & final_red_list, const vector< Vertex > & vList, const char * output_fname, int whichOne)
	{
		switch (whichOne)
		{
		case 0: // Cubical saver
			saveReduction(final_red_list, vList, output_fname);
			break;
		case 1: // Full Rips saver
			saveReduction(final_red_list, output_fname);
			break;
		case 2: // General Simplicial Complex saver
			saveReduction(final_red_list, output_fname);
			break;
		default:
			break;
		}
	}

	void saveReduction(const vector< MatrixListType > & final_red_list, const vector< Vertex > & vList, const char * output_fname)
	{	
		int count = 0;
		for( int j = 0; j < final_red_list.size(); j ++ )
		{
			assert( ! final_red_list[j].empty() );
			count += 1;
			count += final_red_list[j].size();
		}

		unsigned int *redArray = new unsigned int[count * d];
		int index = 0;
		vector<int> header(d);
		header[0] = final_red_list.size();

		for (vector< MatrixListType >::const_iterator red_list_iter = final_red_list.begin(); red_list_iter!=final_red_list.end(); red_list_iter++)
		{
			// size of the reduction list of this dot
			assert( red_list_iter->size() != 0 );
			redArray[ index++ ] = red_list_iter->size();
			for( int j = 1; j < d; j ++ )
				redArray[ index++ ] = 0;

			// for each vertex in the red list, write its coordinates
			for (MatrixListType::const_iterator cellid_iter = red_list_iter->begin(); cellid_iter != red_list_iter->end(); cellid_iter++)
			{
				Vertex coord = vList[ *cellid_iter ];	

				for (size_t ii = 0; ii < d; ii++)
				{
					redArray[ index++ ] = coord[ d-1-ii ] + 1;
				}
			}
		}

		saveResults(output_fname, redArray, header, index);

		delete[] redArray;
	}


	void saveReduction(const vector< MatrixListType > & final_red_list, const char * output_fname)
	{
		int count = 0;
		for (int j = 0; j < final_red_list.size(); j++)
		{
			assert(!final_red_list[j].empty());
			count += 1;
			count += final_red_list[j].size();
		}

		unsigned int *redArray = new unsigned int[count];
		int index = 0;
		vector<int> header(1);
		header[0] = final_red_list.size();

		for (vector< MatrixListType >::const_iterator red_list_iter = final_red_list.begin(); red_list_iter != final_red_list.end(); red_list_iter++)
		{
			assert(red_list_iter->size() != 0);
			redArray[index++] = red_list_iter->size(); // save the current reduction list size
			
			for (MatrixListType::const_iterator cellid_iter = red_list_iter->begin(); cellid_iter != red_list_iter->end(); cellid_iter++)
			{
				redArray[index++] = *cellid_iter + 1; // save the point index
			}
		}

		saveResults(output_fname, redArray, header, index);

		delete[] redArray;
	}
};


// Save persistence values and reduction results to binary file
void saveResults(const char *fileName, unsigned int *data, const vector<int> &header, int count)
{
	FILE *outFile;

	outFile = fopen(fileName, "wb");
	unsigned int dim = header.size();

	fwrite((void*)&dim, sizeof(unsigned int), 1, outFile);

	fwrite((void*)&header[0], sizeof(unsigned int), dim, outFile);

	int writeCount = fwrite((void *)data, sizeof(unsigned int), count, outFile);

	printf("Writing %d Data Reduction Results in File %s\n", writeCount, fileName);

	fclose(outFile);
}
#endif
