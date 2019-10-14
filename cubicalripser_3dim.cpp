/*
CubicalRipser: C++ system for computation of Cubical persistence pairs
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
CubicalRipser is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

CubicalRipser is deeply depending on 'Ripser', software for Vietoris-Rips 
persitence pairs by Ulrich Bauer, 2015-2016.  We appreciate Ulrich very much.
We rearrange his codes of Ripser and add some new ideas for optimization on it 
and modify it for calculation of a Cubical filtration.

This part of CubicalRiper is a calculator of cubical persistence pairs for 
3 dimensional pixel data. The input data format conforms to that of DIPHA.
 See more descriptions in README.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#define FILE_OUTPUT

#include <fstream>
#include <iostream>
#include <algorithm>
#include <queue>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>

using namespace std;

#include "birthday_index.h"
#include "dense_cubical_grids.h"
#include "columns_to_reduce.h"
#include "simplex_coboundary_enumerator.h"
#include "union_find.h"
#include "write_pairs.h"
#include "joint_pairs.h"
#include "compute_pairs.h"


enum calculation_method { LINKFIND, COMPUTEPAIRS};

void print_usage_and_exit(int exit_code) {
	 cerr << "Usage: "
	      << "CR3 "
	      << "[options] [input_filename]" << endl
	      << endl
	      << "Options:" << endl
	      << endl
	      << "  --help           print this screen" << endl
	      << "  --format         use the specified file format for the input. Options are:" << endl
	      << "                     dipha          (3d array in DIPHA file format; default)" << endl
	      << "                     perseus        (3d array in Perseus file format)" << endl
	      << "                     numpy        (3d array in numpy file format)" << endl
	      << "  --threshold <t>  compute cubical complexes up to birth time <t>" << endl
	      << "  --maxdim <t>  compute persistent homology up to dimension <t>" << endl
	      << "  --method         method to compute the persistent homology of the cubical complexes. Options are" << endl
	      << "                     link_find      (calculating the 0-dim PP, use 'link_find' algorithm; default)" << endl
	      << "                     compute_pairs  (calculating the 0-dim PP, use 'compute_pairs' algorithm)" << endl
	      << "  --output         name of file that will contain the persistence diagram " << endl
	      << "  --print          print persistence pairs on your console" << endl
	      << "  --location          output birth location (currently, only for dim 0)" << endl
	      << endl;

	exit(exit_code);
}


int main(int argc, char** argv){

	const char* filename = nullptr;
	string output_filename = "answer_3dim.diagram"; //default output filename
	file_format format = DIPHA;
	calculation_method method = LINKFIND;
	double threshold = 99999;
	int maxdim = 2;  // compute PH up to this dimension
	bool print = false;
	bool location = false;

	for (int i = 1; i < argc; ++i) {
		const string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--threshold") {
			string parameter = string(argv[++i]);
			size_t next_pos;
			threshold = stod(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--maxdim") {
			maxdim = stoi(argv[++i]);
		} else if (arg == "--format") {
			string parameter = string(argv[++i]);
			if (parameter == "dipha") {
				format = DIPHA;
			} else if (parameter == "perseus") {
				format = PERSEUS;
			} else if (parameter == "numpy") {
				format = NUMPY;
			} else {
				print_usage_and_exit(-1);
			}
		} else if(arg == "--method") {
			string parameter = string(argv[++i]);
			if (parameter == "link_find") {
				method = LINKFIND;
			} else if (parameter == "compute_pairs") {
				method = COMPUTEPAIRS;
			} else {
				print_usage_and_exit(-1);
			}
		} else if (arg == "--output") {
			output_filename = string(argv[++i]);
		} else if (arg == "--print"){
			print = true;
		} else if (arg == "--location"){
			location = true;
		} else {
			if (filename) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}

	if (!filename) { print_usage_and_exit(-1); }
    ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		cerr << "couldn't open file " << filename << endl;
		exit(-1);
	}

	vector<WritePairs> writepairs; // dim birth death
	writepairs.clear();
	
	DenseCubicalGrids* dcg = new DenseCubicalGrids(filename, threshold, format);
	ColumnsToReduce* ctr = new ColumnsToReduce(dcg);
	
	switch(method){
		case LINKFIND:
		{
			JointPairs* jp = new JointPairs(dcg, ctr, writepairs, print);
			jp -> joint_pairs_main(); // dim0
			if(maxdim>0){
				ComputePairs* cp = new ComputePairs(dcg, ctr, writepairs, print);
				cp -> compute_pairs_main(); // dim1
				cp -> assemble_columns_to_reduce();
				if(maxdim>1){			
					cp -> compute_pairs_main(); // dim2
				}
			}
		break;
		}
		
		case COMPUTEPAIRS:
		{
			ComputePairs* cp = new ComputePairs(dcg, ctr, writepairs, print);
			cp -> compute_pairs_main(); // dim0
			if(maxdim>0){
				cp -> assemble_columns_to_reduce();
				cp -> compute_pairs_main(); // dim1
				if(maxdim>1){			
					cp -> assemble_columns_to_reduce();
					cp -> compute_pairs_main(); // dim2
				}
			}
		break;
		}
	}

#ifdef FILE_OUTPUT
	ofstream writing_file;
	writing_file.open(output_filename, ios::out);
	if(!writing_file.is_open()){
		cerr << " error: open file for output failed! " << endl;
	}
	int64_t p = writepairs.size();
	cout << "the number of pairs : " << p << endl;
	
	string extension = ".csv";
	if(equal(extension.rbegin(), extension.rend(), output_filename.rbegin()) == true){
		for(int64_t i = 0; i < p; ++i){
			int64_t d = writepairs[i].getDimension();
			writing_file << d << ",";
			writing_file << writepairs[i].getBirth() << ",";
			writing_file << writepairs[i].getDeath();
			if(location && d==0){
				writing_file << "," << writepairs[i].getBirthX() << "," << writepairs[i].getBirthY()<< "," << writepairs[i].getBirthZ();
			}
			writing_file << endl;
		}
		writing_file.close();
	} else {
		int64_t mn = 8067171840;
		writing_file.write((char *) &mn, sizeof( int64_t )); // magic number
		int64_t type = 2;
		writing_file.write((char *) &type, sizeof( int64_t )); // type number of PERSISTENCE_DIAGRAM
		writing_file.write((char *) &p, sizeof( int64_t )); // number of points in the diagram p
		for(int64_t i = 0; i < p; ++i){
			int64_t writedim = writepairs[i].getDimension();
			writing_file.write((char *) &writedim, sizeof( int64_t )); // dim

			double writebirth = writepairs[i].getBirth();
			writing_file.write((char *) &writebirth, sizeof( double )); // birth
			
			double writedeath = writepairs[i].getDeath();
			writing_file.write((char *) &writedeath, sizeof( double )); // death
		}
		writing_file.close();
	}
#endif

	return 0;
}