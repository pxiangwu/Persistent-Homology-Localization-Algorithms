#ifndef INCLUDED_INPUT_RUNNER_H
#define INCLUDED_INPUT_RUNNER_H

#include "Globals.h"
#include "InputFileInfo.h"
#include "DataReaders/DataReaderCubical.h"
#include "DataReaders/DataReaderFullRips.h"
#include "DataReaders/DataReaderSimComplex.h"
#include "PersistenceCalcRunner.h"

template<int dim>
struct InputRunnerCubical
{
	static void run(const InputFileInfo &info, double pers_thd) 
	{		
		blitz::Array<double, dim> phi;
		
		if (info.binary)
		{
			RawDataReaderCubical<dim, double> reader;
			reader.read(info.input_path, phi);
		} 
		else 
		{
			TextDataReaderCubical<dim, double> reader;
			reader.read(info.input_path, phi);
		}

		PersistenceCalcRunnerCubical<dim> calc; 
		calc.go(&phi, pers_thd, info);
	}
};

template<int maxDim>
struct InputRunnerFullRips
{
	static void run(const InputFileInfo &info, double pers_thd)
	{
		blitz::Array<double, 2> distMatrix;

		if (info.binary)
		{
			RawDataReaderFullRips reader;
			reader.read(info.input_path, distMatrix);
		}
		else
		{
			TextDataReaderFullRips reader;
			reader.read(info.input_path, distMatrix);
		}

		PersistenceCalcRunnerFullRips<maxDim> calc;
		calc.go(&distMatrix, pers_thd, info);
	}
};

template<int dim>
struct InputRunnerSimComplex
{
	static void run(const InputFileInfo &info, double pers_thd)
	{
		blitz::Array<double, 1> pointsVal;

		if(info.binary)
		{
			RawDataReaderSimComplex reader;
			reader.read(info.input_path, pointsVal);
		}
		else
		{
			TextDataReaderSimComplex reader;
			reader.read(info.input_path, pointsVal);
		}

		PersistenceCalcRunnerSimComplex<dim> calc;
		calc.go(&pointsVal, pers_thd, info);
	}
};



// This struct is used for solving the inconvenience of template compile-time code generation
template<int x, int to>
struct static_for_InputRunnerFullRips
{
	void operator()(int whichDim, const InputFileInfo & input_file_info)
	{
		if (whichDim == x)
			InputRunnerFullRips<2>::run(input_file_info, Globals::reduction_threshold);

		static_for_InputRunnerFullRips<x + 1, to>()(whichDim, input_file_info);
	}
};

template<int to>
struct static_for_InputRunnerFullRips<to, to>
{
	void operator()(int whichDim, const InputFileInfo & input_file_info)
	{}
};


// launch persistence homology calculation
void runPersistenceHomology(const InputFileInfo & input_file_info)
{
	if (input_file_info.dimension > 8 || Globals::max_dim > 8)
	{
		std::cerr << "\nThis code currently cannot deal with dimension higher than 8.\n";
		std::cerr << "Please modify the code of the function runPersistenceHomology() in file InputRunner.h\n";
		exit(EXIT_FAILURE);
	}

	if (input_file_info.file_type == Globals::FileType::IMAGE_DATA)
	{
		switch (input_file_info.dimension) // This is one approach to deal with template compile-time code generation, i.e., use switch
		{
		case 1:
			InputRunnerCubical<1>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 2:
			InputRunnerCubical<2>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 3:
			InputRunnerCubical<3>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 4:
			InputRunnerCubical<4>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 5:
			InputRunnerCubical<5>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 6:
			InputRunnerCubical<6>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 7:
			InputRunnerCubical<7>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 8:
			InputRunnerCubical<8>::run(input_file_info, Globals::reduction_threshold);
			break;
		}
	}
	else if (input_file_info.file_type == Globals::FileType::DENSE_DISTANCE_MATRIX)
	{
		// This is another approach to dealing with the template inconvenience. 
		// Note that if the range is too large (e.g., [1 100]), the compilation would take a lot of time.
		// The following code can deal with dimension from 1 to 8.
		static_for_InputRunnerFullRips<1, 9>()(Globals::max_dim, input_file_info);
	}
	else if (input_file_info.file_type == Globals::FileType::GENERAL_SIMPLICIAL_COMPLEX)
	{
		switch (input_file_info.dimension)
		{
		case 1:
			InputRunnerSimComplex<1>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 2:
			InputRunnerSimComplex<2>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 3:
			InputRunnerSimComplex<3>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 4:
			InputRunnerSimComplex<4>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 5:
			InputRunnerSimComplex<5>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 6:
			InputRunnerSimComplex<6>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 7:
			InputRunnerSimComplex<7>::run(input_file_info, Globals::reduction_threshold);
			break;
		case 8:
			InputRunnerSimComplex<8>::run(input_file_info, Globals::reduction_threshold);
			break;
		}
	}
}
#endif