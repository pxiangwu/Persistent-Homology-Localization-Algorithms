#ifndef DATA_READER_CUBICAL_H
#define DATA_READER_CUBICAL_H


template<int dim, typename t>
struct TextDataReaderCubical
{
	void read(const string &file_name, blitz::Array<t, dim> &arr)
	{
		std::ifstream str(file_name.c_str());

		// eat the header infomation
		int fileType, d;
		str >> fileType >> d;

		cout << "reading " << dim << "-dimensional datafile" << endl;
		blitz::TinyVector<int, dim> dims;
		cout << "DATA Dims: ";
		for (int i = 0; i < dim; i++) {
			str >> dims[i];
			cout << dims[i] << " ";
		}
		cout << endl;

		cout << "found dimension: " << dims << endl;

		arr.resize(dims);

		for (typename blitz::Array<t, dim>::iterator it = arr.begin(), end = arr.end(); it != end; ++it)
		{
			str >> *it;
		}

		cout << "read the input" << endl;

		str.close();
	}
};

template<int dim, typename t>
struct RawDataReaderCubical
{
	// template<typename stream_t>
	void read(const string &file_name, blitz::Array<t, dim> &arr)
	{
		std::ifstream f(file_name.c_str(), std::ios::binary | std::ios::in);

		typedef double ElemT;

		typedef unsigned int HeaderElemT;

		// eat the header infomation
		int fileType, d;
		f.read(reinterpret_cast<char*>(&fileType), sizeof(int));
		f.read(reinterpret_cast<char*>(&d), sizeof(int));


		blitz::TinyVector<HeaderElemT, dim> cnt;

		f.read(reinterpret_cast<char*>(cnt.data()), sizeof(HeaderElemT) * dim); // read dimension infomation

		blitz::Array<ElemT, dim> a(cnt);
		f.read(reinterpret_cast<char*>(a.data()), sizeof(ElemT)*a.size());

		arr.resize(cnt);
		copy(a.begin(), a.end(), arr.begin());

		f.close();
	}
};


#endif // !DATA_READER_CUBICAL_H

