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
#include <cassert>
#include <chrono>

using namespace std;

#include "cube.h"
// switch T and V
#include "dense_cubical_grids.h"
#include "write_pairs.h"
#include "joint_pairs.h"
#include "compute_pairs.h"
#include "config.h"
#include "npy.hpp"


void print_usage_and_exit(int exit_code) {
	 cerr << "Usage: "
	      << "cubicalripser "
	      << "[options] [input_filename]" << endl
	      << endl
	      << "Options:" << endl
	      << endl
	      << "  --help, -h          print this screen" << endl
	      << "  --verbose, -v       " << endl
	      << "  --threshold <t>, -t compute cubical complexes up to birth time <t>" << endl
		  << "  --maxdim <t>, -m    compute persistent homology up to dimension <t>" << endl
	      << "  --algorithm, -a     algorithm to compute the 0-dim persistent homology:" << endl
	      << "                    		link_find      (default)" << endl
	      << "                    		compute_pairs  (slow in most cases)" << endl
	      << "  --min_recursion_to_cache, -mc  minimum number of recursion for a reduced column to be cached (the higher the slower but less memory)" << endl
	      << "  --cache_size, -c	maximum number of reduced columns to be cached (the lower the slower but less memory)" << endl
	      << "  --output, -o        name of the output file" << endl
	      << "  --print, -p         print persistence pairs on console" << endl
	      << "  --top_dim        	(not recommended) compute only for top dimension using Alexander duality (setting '--maxdim 0 --embedded' is generally faster for this purpose)" << endl
	      << "  --embedded, -e   	Take the Alexander dual (pad the image boundary with -infty and negate the pixel values)" << endl
	      << "  --location, -l   	whether creator/destroyer location is included in the output:" << endl
	      << "								yes      (default)" << endl
	      << "								none     " << endl
	      << endl;

	exit(exit_code);
}

/////////////////////////////////////////////
int main(int argc, char** argv){

	Config config;
	bool arg_embedded = false;
	// command-line argument parsing
	for (int i = 1; i < argc; ++i) {
		const string arg(argv[i]);
		if (arg == "--help" || arg == "-h") {
			print_usage_and_exit(0);
		} else if (arg == "--verbose" || arg == "-v") {
			config.verbose = true;
		} else if (arg == "--threshold" || arg == "-t") {
			string parameter = string(argv[++i]);
			size_t next_pos;
			config.threshold = stod(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--maxdim" || arg == "-m") {
			config.maxdim = stoi(argv[++i]);
		} else if(arg == "--algorithm" || arg == "-a") {
			string parameter = string(argv[++i]);
			if (parameter == "link_find") {
				config.method = LINKFIND;
			} else if (parameter == "compute_pairs") {
				config.method = COMPUTEPAIRS;
			} else {
				print_usage_and_exit(-1);
			}
		} else if (arg == "--output" || arg == "-o") {
			config.output_filename = string(argv[++i]);
		} else if (arg == "--min_recursion_to_cache" || arg == "-mc"){
            config.min_recursion_to_cache = stoi(argv[++i]);
		} else if (arg == "--cache_size" || arg == "-c"){
            config.cache_size = stoi(argv[++i]);
		} else if (arg == "--print" || arg == "-p"){
			config.print = true;
		} else if (arg == "--embedded" || arg == "-e"){
			arg_embedded = true;
		} else if (arg == "--top_dim") {
			config.method = ALEXANDER;
		} else if (arg == "--location" || arg == "-l"){
			string parameter = string(argv[++i]);
			if (parameter == "none") {
				config.location = LOC_NONE;
			}
		} else {
			if (!config.filename.empty()) { print_usage_and_exit(-1); }
			config.filename = argv[i];
		}
	}

	if (config.filename.empty()) { print_usage_and_exit(-1); }
	if(config.method == ALEXANDER){
		config.embedded = !arg_embedded;
	}else{
		config.embedded = arg_embedded;
	}
    ifstream file_stream(config.filename);
	if (!config.filename.empty() && file_stream.fail()) {
		cerr << "couldn't open file " << config.filename << endl;
		exit(-1);
	}
	// infer input file type from its extention
	if(config.filename.find(".txt")!= std::string::npos){
		config.format = PERSEUS;
	}else if(config.filename.find(".npy")!= std::string::npos){
		config.format = NUMPY;
	}else if(config.filename.find(".csv")!= std::string::npos){
		config.format = CSV;
	}else if(config.filename.find(".complex")!= std::string::npos){
		config.format = DIPHA;
	}else{
		cerr << "unknown input file format! (the filename extension should be one of npy, txt, complex): " << config.filename << endl;
		exit(-1);
	}

	vector<WritePairs> writepairs; // record (dim birth death x y z)
	writepairs.clear();
	
	DenseCubicalGrids* dcg = new DenseCubicalGrids(config);
	vector<Cube> ctr;

	// compute PH
    vector<uint64_t> betti(0);
	switch(config.method){
		case LINKFIND:
		{
			dcg->loadImage(config.embedded);
			config.maxdim = std::min<uint8_t>(config.maxdim, dcg->dim - 1);
            auto start0 = std::chrono::system_clock::now();
			JointPairs* jp = new JointPairs(dcg, writepairs, config);
			if(dcg->dim==1){
				jp -> enum_edges({0},ctr);
			}else if(dcg->dim==2){
				jp -> enum_edges({0,1},ctr);
			}else{
				jp -> enum_edges({0,1,2},ctr);
			}
            // dim0
			jp -> joint_pairs_main(ctr,0);
            auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start0).count();
            betti.push_back(writepairs.size());
            cout << "the number of pairs in dim 0: " << betti[0] << endl;
			if(config.verbose){
	            cout << "computation took " << msec << " [msec]" << endl;
			}
            // dim 1 and 2
			if(config.maxdim>0){
                auto start1 = std::chrono::system_clock::now();
				ComputePairs* cp = new ComputePairs(dcg, writepairs, config);
				cp -> compute_pairs_main(ctr); // dim1
                betti.push_back(writepairs.size() - betti[0]);
                auto msec1 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start1).count();
                cout << "the number of pairs in dim 1: " << betti[1] << endl;
				if(config.verbose){
    	            cout << "computation took " << msec1 << " [msec]" << endl;
				}
				if(config.maxdim>1){
                    auto start2 = std::chrono::system_clock::now();
					cp -> assemble_columns_to_reduce(ctr,2);
					cp -> compute_pairs_main(ctr); // dim2
                    auto msec2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start2).count();
                    betti.push_back(writepairs.size() - betti[0] - betti[1]);
                    cout << "the number of pairs in dim 2: " << betti[2] << endl;
					if(config.verbose){
        	            cout << "computation took " << msec2 << " [msec]" << endl;
					}
				}
			}
            auto mseca = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start0).count();
            cout << "the whole computation took " << mseca << " [msec]" << endl;
		break;
		}
		
		case COMPUTEPAIRS:
		{
			dcg->loadImage(config.embedded);
			config.maxdim = std::min<uint8_t>(config.maxdim, dcg->dim - 1);
			ComputePairs* cp = new ComputePairs(dcg, writepairs, config);
			cp -> assemble_columns_to_reduce(ctr,0);
			cp -> compute_pairs_main(ctr); // dim0
            betti.push_back(writepairs.size());
            cout << "the number of pairs in dim 0: " << betti[0] << endl;
			if(config.maxdim>0){
				cp -> assemble_columns_to_reduce(ctr,1);
				cp -> compute_pairs_main(ctr); // dim1
				betti.push_back(writepairs.size() - betti[0]);
				cout << "the number of pairs in dim 1: " << betti[1] << endl;
				if(config.maxdim>1){
					cp -> assemble_columns_to_reduce(ctr,2);
					cp -> compute_pairs_main(ctr); // dim2
					betti.push_back(writepairs.size() - betti[0] - betti[1]);
					cout << "the number of pairs in dim 2: " << betti[2] << endl;
				}
			}
		break;
		}

		case ALEXANDER: // only for top dim
		{
			if(config.tconstruction){
				cerr << "Alexander duality for T-construction is not implemented yet." << endl;
				exit(-9);
			}
			dcg->loadImage(config.embedded);
            auto start0 = std::chrono::system_clock::now();
			JointPairs* jp = new JointPairs(dcg, writepairs, config);
			if(dcg->dim==1){
				jp -> enum_edges({0},ctr);
				jp -> joint_pairs_main(ctr,0); // dim0
				cout << "the number of pairs in dim 0: " << writepairs.size() << endl;
			}else if(dcg->dim==2){
				jp -> enum_edges({0,1,3,4},ctr);
				jp -> joint_pairs_main(ctr,1); // dim1
				cout << "the number of pairs in dim 1: " << writepairs.size() << endl;
			}else if(dcg->dim==3){
				jp -> enum_edges({0,1,2,3,4,5,6,7,8,9,10,11,12},ctr);
				jp -> joint_pairs_main(ctr,2); // dim2
				cout << "the number of pairs in dim 2: " << writepairs.size() << endl;
			}
            auto mseca = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start0).count();
            cout << "computation took " << mseca << " [msec]" << endl;
		break;
		}		
	}

	// determine shift between dcg and the voxel coordinates
	auto pad_x = (dcg->ax - dcg->img_x)/2;
	auto pad_y = (dcg->ay - dcg->img_y)/2;
	auto pad_z = (dcg->az - dcg->img_z)/2;
	// cout << "dcg: " << dcg->ax <<", " << dcg->ay <<", " <<dcg->az <<", " << endl;
	// cout << "image: " << dcg->img_x <<", " << dcg->img_y <<", " <<dcg->img_z <<", " << endl;
	// cout << "padding: " << pad_x <<", " << pad_y <<", " <<pad_z <<", " << endl;

	// write to file
	ofstream writing_file;
	int64_t p = writepairs.size();
	cout << "the number of total pairs : " << p << endl;
	if(config.output_filename.find(".csv")!= std::string::npos){
		writing_file.open(config.output_filename, ios::out);
		if(!writing_file.is_open()){
			cerr << " error: open file for output failed! " << endl;
			exit(-3);
		}
		for(int64_t i = 0; i < p; ++i){
			int64_t d = writepairs[i].dim;
			writing_file << d << "," << writepairs[i].birth << "," << writepairs[i].death;
			if(config.location != LOC_NONE){
				writing_file << "," << writepairs[i].birth_x-pad_x << "," << writepairs[i].birth_y-pad_y << "," << writepairs[i].birth_z-pad_z;
				writing_file << "," << writepairs[i].death_x-pad_x << "," << writepairs[i].death_y-pad_y << "," << writepairs[i].death_z-pad_z;
			}
			writing_file << endl;
		}
		writing_file.close();
	}else if(config.output_filename.find(".npy")!= std::string::npos){ // output in npy
		long unsigned leshape[] = {0,9};
		leshape[0] = p;
		vector<double> data(9*p);
		for(int64_t i = 0; i < p; ++i){
			data[6*i] = writepairs[i].dim;
			data[6*i+1] = writepairs[i].birth;
			data[6*i+2] = writepairs[i].death;
			data[6*i+3] = writepairs[i].birth_x-pad_x;
			data[6*i+4] = writepairs[i].birth_y-pad_y;
			data[6*i+5] = writepairs[i].birth_z-pad_z;
			data[6*i+6] = writepairs[i].death_x-pad_x;
			data[6*i+7] = writepairs[i].death_y-pad_y;
			data[6*i+8] = writepairs[i].death_z-pad_z;
		}
		try{
			npy::SaveArrayAsNumpy(config.output_filename, false, 2, leshape, data);
		} catch (...) {
			cerr << " error: open file for output failed! " << endl;
		}
	}else if(config.output_filename.compare("none")==0){ // no output
		return 0;
	} else { // output in DIPHA format
		writing_file.open(config.output_filename, ios::out | ios::binary);
		if(!writing_file.is_open()){
			cerr << " error: open file for output failed! " << endl;
			exit(-3);
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
