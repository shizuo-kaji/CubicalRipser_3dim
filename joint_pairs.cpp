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


// enumerate all edges (dim=1)
JointPairs::JointPairs(DenseCubicalGrids* _dcg, vector<Cube>& ctr, vector<WritePairs> &_wp, Config& _config){
	dcg = _dcg;
	config = &_config;
	wp = &_wp;
	ctr.clear();
	// the order of loop matters for performance!
	for (uint8_t m = 0; m < 3; ++m) {
		for (uint32_t z = 0; z < dcg->az; ++z) {
			for (uint32_t y = 0; y < dcg->ay; ++y) {
				for(uint32_t x = 0; x < dcg->ax ; ++x){
					double birth = dcg -> getBirth(x,y,z,m, 1);
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
void JointPairs::joint_pairs_main(vector<Cube>& ctr){
	UnionFind dset(dcg);
	uint64_t u,v=0;
	double min_birth = config->threshold;
	uint64_t min_idx=0;

	if(config->print == true){
		cout << "persistence intervals in dim " << 0 << ":" << endl;
	}
	
    for (auto e = ctr.rbegin(), last = ctr.rend(); e != last; ++e) {
		// indexing scheme for union find is DIFFERENT from that of cubes
		// for each edge e, identify root indices u and v of the end points
		uint64_t uind = e->x() + (dcg->ax)*e->y() + (dcg->axy)*e->z();
		uint64_t vind;
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
		}
		v = dset.find(vind); 

		if(u != v){
			double birth;
			int rcind;
			if(dset.birthtime[u] >= dset.birthtime[v]){
				birth = dset.birthtime[u]; // the component u is killed
				if(config->location==LOC_DEATH){
					rcind = uind; // cell of which the location is recorded
				}else{
					rcind = u; // cell of which the location is recorded
				}
				if (dset.birthtime[v] < min_birth) {
					min_birth = dset.birthtime[v];
					min_idx = v;
				}
			}else{
				birth = dset.birthtime[v]; 
				if(config->location==LOC_DEATH){
					rcind = vind; // cell of which the location is recorded
				}else{
					rcind = v; // cell of which the location is recorded
				}
				if (dset.birthtime[u] < min_birth) {
					min_birth = dset.birthtime[u];
					min_idx = u;
				}
			}
			double death = e->birth;
			dset.link(u, v);
			if(birth != death){
				wp -> push_back(WritePairs(0, birth, death, rcind%(dcg->ax), (rcind/(dcg->ax))%(dcg->ay), (rcind/(dcg->axy))%(dcg->az), config->print));
			}
			// column clearing
	        e->index = NONE;
		}
	}
	// the base point component
	wp -> push_back(WritePairs(0, min_birth, dcg -> threshold, min_idx%(dcg->ax), (min_idx/(dcg->ax))%(dcg->ay), (min_idx/(dcg->axy))%(dcg->az), config->print));
//	cout << ctr.size() << endl;

	// remove unnecessary edges
	auto new_end = std::remove_if(ctr.begin(), ctr.end(),
                              [](const Cube& e){ return e.index == NONE; });
	ctr.erase(new_end, ctr.end());
//	cout << ctr.size() << endl;
//	std::sort(ctr.begin(), ctr.end(), CubeComparator()); // we can skip sorting as it is already sorted
}
