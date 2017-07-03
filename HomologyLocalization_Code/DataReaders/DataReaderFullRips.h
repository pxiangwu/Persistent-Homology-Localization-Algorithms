#ifndef DATA_READER_FULL_RIPS_H
#define DATA_READER_FULL_RIPS_H

#include <vector>
using namespace std;

struct RawDataReaderFullRips
{
	void read(const string &file_name, blitz::Array<double, 2> &arr)
	{
		std::ifstream f(file_name.c_str(), std::ios::binary | std::ios::in);

		int fileType, numPoints, dimPoints;
		f.read(reinterpret_cast<char*>(&fileType), sizeof(int));	 // read file type
		f.read(reinterpret_cast<char*>(&numPoints), sizeof(int)); // read number of points
		f.read(reinterpret_cast<char*>(&dimPoints), sizeof(int)); // read dimension of each point

		vector<vector<double>> posPoints;
		posPoints.resize(numPoints, vector<double>(dimPoints)); // read point positions
		for (int i = 0; i < numPoints; ++i)
		{
			f.read(reinterpret_cast<char*>(posPoints[i].data()), sizeof(double)*dimPoints);
		}

		arr.resize(numPoints, numPoints);
		f.read(reinterpret_cast<char*>(arr.data()), sizeof(double)*arr.size()); // read distance matrix
		
		f.close();
	}
};


struct TextDataReaderFullRips
{
	void read(const string &file_name, blitz::Array<double, 2> &arr)
	{
		std::ifstream f(file_name.c_str());

		int fileType, numPoints, dimPoints;
		f >> fileType >> numPoints >> dimPoints;

		vector<vector<double>> posPoints;
		posPoints.resize(numPoints, vector<double>(dimPoints)); // read point positions
		for (int i = 0; i < numPoints; ++i)
		{
			for (int j = 0; j < dimPoints; ++j)
			{
				f >> posPoints[i][j];
			}
		}

		arr.resize(numPoints, numPoints);
		for (blitz::Array<double, 2>::iterator it = arr.begin(), end = arr.end(); it != end; ++it)
		{
			f >> *it;
		}

		f.close();
	}
};



#endif // !DATA_READER_FULL_RIPS_H

