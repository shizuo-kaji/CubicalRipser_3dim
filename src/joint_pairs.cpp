/* joint_pairs.cpp

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdint>

#include "cube.h"
#include "dense_cubical_grids.h"
#include "coboundary_enumerator.h"
#include "union_find.h"
#include "write_pairs.h"
#include "joint_pairs.h"

using namespace std;


JointPairs::JointPairs(DenseCubicalGrids* _dcg, vector<WritePairs> &_wp, Config& _config){
	dcg = _dcg;
	config = &_config;
	wp = &_wp;
}


// enumerate all edges
void JointPairs::enum_edges(std::vector<uint8_t> types, vector<Cube>& ctr){
	ctr.clear();
	// the order of loop matters for performance!
	for (const auto& m : types) {
		for (uint32_t z = 0; z < dcg->az; ++z) {
			for (uint32_t y = 0; y < dcg->ay; ++y) {
				for(uint32_t x = 0; x < dcg->ax ; ++x){
					double birth = dcg -> getBirth(x,y,z,m, 1);
//                    cout << x << "," << y << "," << z << ", " << m << "," << birth << endl;
					if(birth < config->threshold){
						ctr.push_back(Cube(birth, x,y,z,m));
					}
				}
			}
		}
	}
	std::sort(ctr.begin(), ctr.end(), CubeComparator());
}

// compute H_0 by union find
void JointPairs::joint_pairs_main(vector<Cube>& ctr, int current_dim){
	UnionFind dset(dcg);
	uint64_t u,v=0;
	double min_birth = config->threshold;
	uint64_t min_idx=0;
	
    for (auto e = ctr.rbegin(), last = ctr.rend(); e != last; ++e) {
		// indexing scheme for union find is DIFFERENT from that of cubes
		uint64_t uind = e->x() + (dcg->ax)*e->y() + (dcg->axy)*e->z();
		uint64_t vind;
		// for each edge e, identify root indices u and v of the end points
		u = dset.find(uind);
		switch(e->m()){ //type
			case 0:
				vind = uind+1; // x+1
				break;
			case 1:
				vind = uind+(dcg->ax); // y+1
				break;
			case 2:
				vind = uind+(dcg->axy); // z+1
				break;
            // T-construction
			case 3:
				vind = uind+1+(dcg->ax); // x+1,y+1
				break;
			case 4:
				vind = uind+1-(dcg->ax); // x+1,y-1
				break;
			// 3D T-construction only
			case 5:
				vind = uind-(dcg->ax)+(dcg->axy); // y-1,z+1
				break;
			case 6:
				vind = uind+(dcg->ax)+(dcg->axy); // y+1,z+1
				break;
			case 7:
				vind = uind+1-(dcg->ax)+(dcg->axy); // x+1,y-1,z+1
				break;
			case 8:
				vind = uind+1+(dcg->axy); // x+1,z+1
				break;
			case 9:
				vind = uind+1+(dcg->ax)+(dcg->axy); // x+1,y+1,z+1
				break;
			case 10:
				vind = uind+1-(dcg->ax)-(dcg->axy); // x+1,y-1,z-1
				break;
			case 11:
				vind = uind+1-(dcg->axy); // x+1,z-1
				break;
			case 12:
				vind = uind+1+(dcg->ax)-(dcg->axy); // x+1,y+1,z-1
				break;
			default:
				exit(-1);
		}
//        cout << uind << "," << vind << endl;
		v = dset.find(vind); // (uind,under:u) meets (vind,under:v) 

		if(u != v){
			double birth;
			uint64_t birth_ind,death_ind;
			if(dset.birthtime[u] >= dset.birthtime[v]){
				birth = dset.birthtime[u]; // the younger component u is killed
				if(current_dim==0){
					if(dset.birthtime[uind]>dset.birthtime[vind]){
						death_ind = uind; // cell of which the location is recorded
					}else{
						death_ind = vind;
					}
					birth_ind = u; 
				}else{
					death_ind = u;
					if(dset.birthtime[uind]>dset.birthtime[vind]){
						birth_ind = uind; 
					}else{
						birth_ind = vind;
					}
				}
				if (dset.birthtime[v] < min_birth) {
					min_birth = dset.birthtime[v];
					min_idx = v;
				}
			}else{ // the younger component v is killed
				birth = dset.birthtime[v];
				if(current_dim==0){
					if(dset.birthtime[uind]>dset.birthtime[vind]){
						death_ind = uind; 
					}else{
						death_ind = vind;
					}
					birth_ind = v; 
				}else{
					death_ind = v;
					if(dset.birthtime[uind]>dset.birthtime[vind]){
						birth_ind = uind;
					}else{
						birth_ind = vind;
					}
				}
				if (dset.birthtime[u] < min_birth) {
					min_birth = dset.birthtime[u];
					min_idx = u;
				}
			}
			double death = e->birth;
			dset.link(u, v);
			if(birth != death){
				if(config->tconstruction){
					wp -> push_back(WritePairs(current_dim, Cube(birth, birth_ind%(dcg->ax), (birth_ind/(dcg->ax))%(dcg->ay), (birth_ind/(dcg->axy))%(dcg->az), 0), Cube(death, death_ind%(dcg->ax), (death_ind/(dcg->ax))%(dcg->ay), (death_ind/(dcg->axy))%(dcg->az), 0), dcg, config->print));
				}else{
					wp -> push_back(WritePairs(current_dim, birth, death, birth_ind%(dcg->ax), (birth_ind/(dcg->ax))%(dcg->ay), (birth_ind/(dcg->axy))%(dcg->az), death_ind%(dcg->ax), (death_ind/(dcg->ax))%(dcg->ay), (death_ind/(dcg->axy))%(dcg->az),config->print));
				}
			}
			// column clearing
	        e->index = NONE;
		}
	}
	// the base point component
	if(current_dim==0){
		wp -> push_back(WritePairs(current_dim, min_birth, dcg -> threshold, min_idx%(dcg->ax), (min_idx/(dcg->ax))%(dcg->ay), (min_idx/(dcg->axy))%(dcg->az), 0,0,0, config->print));
	}

	// remove unnecessary edges
	if(config->maxdim==0 || current_dim>0){
		return;
	}else{
		auto new_end = std::remove_if(ctr.begin(), ctr.end(),
								[](const Cube& e){ return e.index == NONE; });
		ctr.erase(new_end, ctr.end());
	//	cout << ctr.size() << endl;
	//	std::sort(ctr.begin(), ctr.end(), CubeComparator()); // we can skip sorting as it is already sorted
	}
}
