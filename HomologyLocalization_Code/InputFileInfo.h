#ifndef INPUT_FILE_INFO
#define INPUT_FILE_INFO


struct InputFileInfo
{
	// common info
	bool binary = true;				// if the input file is binary
	string input_path = "";		// input file path
	int file_type = 0;				// input file type (0: Image data; 1: Dense distance matrix; 2: sparse distance matrix)

	// info for cubical image data
	int dimension = 0;				// the data dimension, exclusively for cubical image data and general simplicial complex

	// info for dense distance matrix
	int numPoints = 0;				// number of points
	int dimPoints = 0;				// dimension of each feature point

	explicit InputFileInfo(const string &input_file)
	{
		input_path = input_file;

		if (input_path.find(".txt") != string::npos) // plain text file
		{
			binary = false;
			fstream f(input_path.c_str());
			
			if (!f.is_open())
				cout << "Cannot find the file " << input_path << endl;

			f >> file_type;
			if (file_type == 0) // Image data
			{
				f >> dimension;
				assert(dimension <= 8);
			}
			else if (file_type == 1) // Dense distance matrix
			{
				f >> numPoints;
				f >> dimPoints;
			}
			else if (file_type == 2) // General simplicial complex
			{
				f >> dimension;
				f >> numPoints;
				f >> dimPoints;
			}
			
			f.close();
		}
		else // binary data file
		{
			binary = true;
			fstream f(input_path.c_str(), ios::in | ios::binary);

			if (!f.is_open())
				cout << "Cannot find the file " << input_path << endl;

			f.read(reinterpret_cast<char*>(&file_type), sizeof(int));
			if (file_type == 0) // Image data
			{
				f.read(reinterpret_cast<char*>(&dimension), sizeof(int));
				assert(dimension <= 8);

				cout << "The dimension of input data: " << dimension << endl;
			}
			else if (file_type == 1) // Dense distance matrix
			{
				f.read(reinterpret_cast<char*>(&numPoints), sizeof(int));
				f.read(reinterpret_cast<char*>(&dimPoints), sizeof(int));
			}
			else if (file_type == 2) // General simplicial complex
			{
				f.read(reinterpret_cast<char*>(&dimension), sizeof(int));
				f.read(reinterpret_cast<char*>(&numPoints), sizeof(int));
				f.read(reinterpret_cast<char*>(&dimPoints), sizeof(int));
			}

			f.close();
		}
	}
};


#endif
