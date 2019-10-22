/*
This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


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
#include "simplex_coboundary_enumerator.h"
#include "union_find.h"
#include "write_pairs.h"
#include "joint_pairs.h"
#include "compute_pairs.h"
#include "npy.hpp"

enum calculation_method { LINKFIND, COMPUTEPAIRS};

void print_usage_and_exit(int exit_code) {
	 cerr << "Usage: "
	      << "cubicalripser "
	      << "[options] [input_filename]" << endl
	      << endl
	      << "Options:" << endl
	      << endl
	      << "  --help           print this screen" << endl
	      << "  --threshold <t>  compute cubical complexes up to birth time <t>" << endl
	      << "  --maxdim <t>     compute persistent homology up to dimension <t>" << endl
	      << "  --method         method to compute the persistent homology of the cubical complexes. Options are" << endl
	      << "                     link_find      (calculating the 0-dim PH, use 'link_find' algorithm; default)" << endl
	      << "                     compute_pairs  (calculating the 0-dim PH, use 'compute_pairs' algorithm)" << endl
	      << "  --output         name of file that will contain the persistence diagram " << endl
	      << "  --print          print persistence pairs on your console" << endl
	      << "  --location          output birth location (currently, only for dim 0)" << endl
	      << endl;

	exit(exit_code);
}


int main(int argc, char** argv){

	string filename = "";
	string output_filename = "output.csv"; //default output filename
	file_format format;
	calculation_method method = LINKFIND;
	double threshold = 99999;
	int maxdim = 2;  // compute PH up to this dimension
	bool print = false;
	bool location = false;

	// command-line argument parsing
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
			if (!filename.empty()) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}

	if (filename.empty()) { print_usage_and_exit(-1); }
    ifstream file_stream(filename);
	if (!filename.empty() && file_stream.fail()) {
		cerr << "couldn't open file " << filename << endl;
		exit(-1);
	}
	// infer input file type from its extention
	if(filename.find(".txt")!= std::string::npos){
		format = PERSEUS;
	}else if(filename.find(".npy")!= std::string::npos){
		format = NUMPY;
	}else{
		format = DIPHA;
	}
	vector<WritePairs> writepairs; // (dim birth death x y z)
	writepairs.clear();
	
	DenseCubicalGrids* dcg = new DenseCubicalGrids(filename, threshold, format);
	vector<BirthdayIndex> ctr;

	maxdim = std::min(maxdim, dcg->dim);

	// compute PH
	switch(method){
		case LINKFIND:
		{
			JointPairs* jp = new JointPairs(dcg, writepairs, print);
			jp -> joint_pairs_main(ctr); // dim0
			ComputePairs* cp = new ComputePairs(dcg, writepairs, print);
			if(maxdim>0){
				cp -> compute_pairs_main(ctr); // dim1
				cp -> assemble_columns_to_reduce(ctr,2);
				if(maxdim>1){			
					cp -> compute_pairs_main(ctr); // dim2
				}
			}
		break;
		}
		
		case COMPUTEPAIRS:
		{
			ComputePairs* cp = new ComputePairs(dcg, writepairs, print);
			cp -> assemble_columns_to_reduce(ctr, 0);
			cp -> compute_pairs_main(ctr); // dim0
			if(maxdim>0){
				cp -> assemble_columns_to_reduce(ctr,1);
				cp -> compute_pairs_main(ctr); // dim1
				if(maxdim>1){			
					cp -> assemble_columns_to_reduce(ctr,2);
					cp -> compute_pairs_main(ctr); // dim2
				}
			}
		break;
		}
	}

	// write to file
	ofstream writing_file;
	int64_t p = writepairs.size();
	cout << "the number of pairs : " << p << endl;
	if(output_filename.find(".csv")!= std::string::npos){
		writing_file.open(output_filename, ios::out);
		if(!writing_file.is_open()){
			cerr << " error: open file for output failed! " << endl;
		}
		for(int64_t i = 0; i < p; ++i){
			int64_t d = writepairs[i].dim;
			writing_file << d << ",";
			writing_file << writepairs[i].birth << ",";
			writing_file << writepairs[i].death;
			if(location){
				writing_file << "," << writepairs[i].birth_x << "," << writepairs[i].birth_y<< "," << writepairs[i].birth_z;
			}
			writing_file << endl;
		}
		writing_file.close();
	}else if(output_filename.find(".npy")!= std::string::npos){
		long unsigned leshape[] = {0,6};
		leshape[0] = p;
		vector<double> data(6*p);
		for(int64_t i = 0; i < p; ++i){
			data[6*i] = writepairs[i].dim;
			data[6*i+1] = writepairs[i].birth;
			data[6*i+2] = writepairs[i].death;
			data[6*i+3] = writepairs[i].birth_x;
			data[6*i+4] = writepairs[i].birth_y;
			data[6*i+5] = writepairs[i].birth_z;
		}
		npy::SaveArrayAsNumpy(output_filename, false, 2, leshape, data);
	} else { // DIPHA format
		writing_file.open(output_filename, ios::out | ios::binary);
		if(!writing_file.is_open()){
			cerr << " error: open file for output failed! " << endl;
		}
		int64_t mn = 8067171840;
		writing_file.write((char *) &mn, sizeof( int64_t )); // magic number
		int64_t type = 2;
		writing_file.write((char *) &type, sizeof( int64_t )); // type number of PERSISTENCE_DIAGRAM
		writing_file.write((char *) &p, sizeof( int64_t )); // number of points in the diagram p
		for(int64_t i = 0; i < p; ++i){
			int64_t writedim = writepairs[i].dim;
			writing_file.write((char *) &writedim, sizeof( int64_t ));
			writing_file.write((char *) &writepairs[i].birth, sizeof( double ));
			writing_file.write((char *) &writepairs[i].death, sizeof( double ));
		}
		writing_file.close();
	}

	return 0;
}