# Persistent Homology Localization Algorithms #

Copyright 2017 Rutgers University and CUNY Queens College

## Description: ##

This page contains the code of persistent homology localization algorithm proposed in [1]. For a given homology class in a certain complex, this code computes the optimal (shortest) representative cycle.

Currently, there are three types of input data that are supported by this software package:

1. d-dimensional gray-scale image data. This data is internally discretized into cubical complex.
2. Distance matrix data. This data is internally interpreted as a Vietoris-Rips complex of as many points as there are columns in the given matrix. The distance between any two points is given by the input matrix data.
3. General simplical complex.

The output produced by this software consists of the optimal representative cycles along with the persistence diagram data.

The input and output file formats are specified below. This software package includes Matlab functions to create the input files and interpret output results.

## Setup: ##

1. **Windows**: 

     - This code has been tested on Visual Studio 2015. Since features of C++11 are utilized in this code, older C++ compilers are not supported (e.g., Visual Studio 2013 will not compile successfully).

     - To compile, first need to include the blitz library: `/Third_Party/`. To avoid the deprecation warnings, add `_CRT_SECURE_NO_WARNINGS` and `_SCL_SECURE_NO_WARNINGS` to the compile options, which can be set in Project Properties -> C/C++ -> Preprocessor -> Preprocessor Definitions.

2. **Linux/macOS**:
     - This code has been tested on GCC version 4.8.4, and any version greater than this should also work.
     - To compile, just use the `Makefile` provided.


## Usage: ##

In command line, run: `HomologyLocalization -f data_file_name [options]`. Here, `data_file_name` is the name of input data, which can be simplicial complex, Vietoris-Rips complex or cubical complex. The available options are:

- `-t` or `--threshold`: only homology classes with persistence greater than this threshold parameter will be considered for the computation of optimal cycles.
- `-a` or `--algorithm`: specify which algorithm to employ to find the optimal representative cycles. Currently there are two options: 0, A* search; 1, classical exhaustive search on the whole covering graph. For instance, the following command `HomologyLocalization -f data_file_name -a 0` will apply the A* search algorithm.
- `-d` or `--dimension`: specify the maximum dimension to be considered for computation. For example, by default `d=2`, and this means we will only compute the `1d` and `2d` boundary matrices.
- `-h` or `--help`: show help information.

If no options are offered, the program will not run the cycle optimization algorithm and will only naively reduce the boundary matrix. As a result, it just returns the possibly non-optimal cycles.

## File Formats: ##

For the input data, currently there are three file types that are supported.

1. d-dimensional gray-scale image data (interpreted as cubical complex). A simple Matlab file writer is provided in this software package `/Matlab/Save_Cubical_Image.m`
     - **File type**, which is `0` for cubical complex data.
     - **Data dimension**, which is followed by the exact size in each dimension. For example, for a `2D` image of size `400x400`, we have `2 400 400`.
     - **Pixel values**.
2. Vietoris-Rips complex. A simple Matlab file writer and demo can be found in `/Matlab/Save_Dense_Distance_Matrix.m` and `/Matlab/Demo_Full_Rips.m`.
     - **File type**, which is `1` for dense distance matrix data.
     - **Total number of points**.
     - **The dimension of each point**. For instance, for points in `3D` space, this number should be `3`.
     - **The positions of points**.
     - **Dense distance matrix**.
3. Simplicial complex. A simple Matlab file writer and demo are provided in `/Matlab/Save_General_Simplicial_Complex.m` and `/Matlab/Demo_General_Simplicial_Complex.m`.
     - **File type**, which is `2` for simplicial complex.
     - **Maximum dimension of simplicial complex**.
     - **Total number of vertices**.
     - **The dimension of each vertex**.
     - **Positions of vertices and their corresponding filtration values**.
     - **Write indices of edges, faces, etc ...**, just like adjacency matrix.

For the output, there are also three different kinds of file types: `.pers`, `.red` and `.bnd`. `.pers` stores the information of birth time and death time for each homology class, while `.red` and `.bnd` store the reduction process and the reduced boundary matrix, respectively. The number appended to the file extension `.bnd` indicates the dimension of the reduced boundary matrix, and this is similar for `.red`. For example, `file_name.bnd.1` is for `1D` reduced boundary matrix. The corresponding file readers are implemented in `/Matlab/Read_Pers_Results_Cubical.m`, `/Matlab/Read_Pers_Results_FullRips.m` and `/Matlab/Read_Pers_Results_General_SimComplex.m`.

## Examples: ##

- Command `HomologyLocalization -f filename.dat -a 1 -t 100` will compute the optimal cycles for homology class with persistence greater than 100, using classical exhaustive search.
- Command `HomologyLocalization -f filename.dat` just performs the column-wise Gaussian reduction for the boundary matrix, and returns the possibly non-optimal cycles.

## TODO: ##

- Parallelize the computation of edge annotations.

## License and Disclaimer: ##

Currently released under GPLv3 ([https://www.gnu.org/licenses/gpl.html](https://www.gnu.org/licenses/gpl.html "https://www.gnu.org/licenses/gpl.html")).

*The SOFTWARE PACKAGE provided in this page is provided "as is", without any guarantee made as to its suitability or fitness for any particular use. It may contain bugs, so use of this tool is at your own risk. We take no responsibility for any damage of any sort that may unintentionally be caused through its use.*

## Citation ##

If you find this code helpful, please cite our work [1] with the following bibtex:

    @inproceedings{ipmi/WuCWZYQMA17,
  		author = {Pengxiang Wu and 
               	Chao Chen and
               	Yusu Wang and
               	Shaoting Zhang and
               	Changhe Yuan and
               	Zhen Qian and
               	Dimitris N. Metaxas and
               	Leon Axel},
  		title = {Optimal Topological Cycles and Their Application in Cardiac Trabeculae Restoration},
  		booktitle = {Information Processing in Medical Imaging, {IPMI}},
  		pages = {80--92},
  		year = {2017}
	}

## Contacts: ##

If you have any questions regarding this code, please contact Pengxiang Wu (_pxiangwu@gmail.com_), or just leave a message below with Github (log-in is needed).

## References: ##

[1] P. Wu, C. Chen, Y. Wang, S. Zhang, C. Yuan, Z. Qian, D. Metaxas and L. Axel. "Optimal Topological Cycles and Their Application in Cardiac Trabeculae Restoration." In *International Conference on Information Processing in Medical Imaging (IPMI)*, 2017.

