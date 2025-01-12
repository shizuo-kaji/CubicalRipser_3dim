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

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

#include "cube.h"
#include "write_pairs.h"
#include "joint_pairs.h"
#include "compute_pairs.h"
#include "config.h"
#include "dense_cubical_grids.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

using namespace std;

namespace py = pybind11;

/////////////////////////////////////////////
py::array_t<double> computePH(py::array_t<double> img, int maxdim=0, bool top_dim=false, bool embedded=false, const std::string &location="yes"){
	// we ignore "location" argument
	Config config;
	config.format = NUMPY;

	vector<WritePairs> writepairs; // (dim birth death x y z)
	writepairs.reserve(1000);
	//writepairs.clear();
	
	auto dcg = std::make_unique<DenseCubicalGrids>(config);
	vector<Cube> ctr;

    const auto &buff_info = img.request();
    const auto &shape = buff_info.shape;
	dcg->dim = buff_info.ndim;
	config.maxdim = maxdim;
	config.maxdim = std::min<uint8_t>(config.maxdim, dcg->dim - 1);
	if(top_dim && dcg->dim > 1){
		config.method = ALEXANDER;
		config.embedded = !embedded;
	}else{
		config.embedded = embedded;
	}

	dcg->ax = shape[0];
	dcg->img_x = shape[0];
    dcg->ay = (dcg->dim > 1) ? shape[1] : 1;
    dcg->img_y = dcg->ay;
    dcg->az = (dcg->dim > 2) ? shape[2] : 1;
    dcg->img_z = dcg->az;
	bool fortran_order = img.flags() & py::array::f_style;
	dcg -> gridFromArray(&img.data()[0], embedded, fortran_order);
//	dense3[x][y][z] = -(*img.data(x-2, y-2, z-2));

	// padded boundary
    if(config.tconstruction){
        if(dcg->az>1) dcg->az++;
        dcg->ax++;
        dcg->ay++;
    }

	dcg -> axy = dcg->ax * dcg->ay;
	dcg -> ayz = dcg->ay * dcg->az;
	dcg -> axyz = dcg->ax * dcg->ay * dcg->az;
	

	// compute PH
	if(config.method==ALEXANDER){
		// compute PH
		auto jp = std::make_unique<JointPairs>(dcg.get(), writepairs, config);
		if(dcg->dim==1){
			jp -> enum_edges({0},ctr);
			jp -> joint_pairs_main(ctr,0); // dim0
		}else if(dcg->dim==2){
			jp -> enum_edges({0,1,3,4},ctr);
			jp -> joint_pairs_main(ctr,1); // dim1
		}else if(dcg->dim==3){
			jp -> enum_edges({0,1,2,3,4,5,6,7,8,9,10,11,12},ctr);
			jp -> joint_pairs_main(ctr,2); // dim2
		}
//		delete jp;
	}else{
        auto cp = std::make_unique<ComputePairs>(dcg.get(), writepairs, config);
        auto jp = std::make_unique<JointPairs>(dcg.get(), writepairs, config);
        std::vector<uint32_t> betti;
		if(dcg->dim==1){
			jp -> enum_edges({0},ctr);
		}else if(dcg->dim==2){
			jp -> enum_edges({0,1},ctr);
		}else{
			jp -> enum_edges({0,1,2},ctr);
		}
		jp -> joint_pairs_main(ctr,0); // dim0
		betti.push_back(writepairs.size());
		if(config.maxdim>0){
			cp -> compute_pairs_main(ctr); // dim1
			betti.push_back(writepairs.size() - betti[0]);
			if(config.maxdim>1){
				cp -> assemble_columns_to_reduce(ctr,2);
				cp -> compute_pairs_main(ctr); // dim2
				betti.push_back(writepairs.size() - betti[0] - betti[1]);
			}
		}
		// delete jp;
		// delete cp;
	}
//	delete dcg;

	// result
	// determine shift between dcg and the voxel coordinates
	auto pad_x = (dcg->ax - dcg->img_x)/2;
	auto pad_y = (dcg->ay - dcg->img_y)/2;
	auto pad_z = (dcg->az - dcg->img_z)/2;
	int64_t p = writepairs.size();
	vector<ssize_t> result_shape{p,9};
	py::array_t<double> data{result_shape};
	auto data_ptr = data.mutable_data();
	for(int64_t i = 0; i < p; ++i){
        data_ptr[i * 9 + 0] = writepairs[i].dim;
        data_ptr[i * 9 + 1] = writepairs[i].birth;
        data_ptr[i * 9 + 2] = writepairs[i].death;
        data_ptr[i * 9 + 3] = writepairs[i].birth_x - pad_x;
        data_ptr[i * 9 + 4] = writepairs[i].birth_y - pad_y;
        data_ptr[i * 9 + 5] = writepairs[i].birth_z - pad_z;
        data_ptr[i * 9 + 6] = writepairs[i].death_x - pad_x;
        data_ptr[i * 9 + 7] = writepairs[i].death_y - pad_y;
        data_ptr[i * 9 + 8] = writepairs[i].death_z - pad_z;
	}
	return data;
}