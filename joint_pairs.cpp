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
JointPairs::JointPairs(DenseCubicalGrids* _dcg, vector<Cube>& ctr, vector<WritePairs> &_wp, const bool _print){
	dcg = _dcg;
	print = _print;
	wp = &_wp;
	ctr.clear();
	// the order of loop matters for performance!
	for (int m = 0; m < 3; ++m) {
		for (int z = 0; z < dcg->az; ++z) {
			for (int y = 0; y < dcg->ay; ++y) {
				for(int x = 0; x < dcg->ax ; ++x){
					double birthday = dcg -> getBirthday(x,y,z,m, 1);
					if(birthday < dcg -> threshold){
						long index = dcg->getIndex(x, y, z, m);
						ctr.push_back(Cube(birthday, index));
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
	long u,v=0;
	double min_birth = dcg -> threshold;
	long min_idx=0;

	if(print == true){
		cout << "persistence intervals in dim " << 0 << ":" << endl;
	}
	
    for (auto e = ctr.rbegin(), last = ctr.rend(); e != last; ++e) {
		// we have to modify here when indexing scheme is changed
		// identify end points
		long ind = e->index % dcg->axyz;
		u = dset.find(ind);
		switch(e->index / dcg->axyz){ //type
			case 0:
				v = dset.find(ind+1); // x+1
				break;
			case 1:
				v = dset.find(ind+(dcg->ax)); // y+1
				break;
			case 2:
				v = dset.find(ind+(dcg->axy)); // z+1
				break;
		}

		if(u != v){
			double birth;
			long idx;
			if(dset.birthtime[u] >= dset.birthtime[v]){
				birth = dset.birthtime[u];
				idx = u; // the one who dies to make a cycle
				if (dset.birthtime[v] < min_birth) {
					min_birth = dset.birthtime[v];
					min_idx = v;
				}
			}else{
				birth = dset.birthtime[v]; 
				idx = v; // the one who dies to make a cycle
				if (dset.birthtime[u] < min_birth) {
					min_birth = dset.birthtime[u];
					min_idx = u;
				}
			}
			double death = e->birthday;
			dset.link(u, v);
			if(birth != death){
				vector<int> loc(dcg->getXYZM(idx));
				wp -> push_back(WritePairs(0, birth, death, loc[0], loc[1], loc[2], print));
			}
			// If two values have same parent, these are potential edges which make a 2-simplex. Otherwise, remove.
	        e->index = -2;
		}
	}
	// the base point component
	vector<int> loc(dcg->getXYZM(min_idx));
	wp -> push_back(WritePairs(0, min_birth, dcg -> threshold, loc[0],loc[1],loc[2],print));
//	cout << ctr.size() << endl;

	// remove unnecessary edges
	auto new_end = std::remove_if(ctr.begin(), ctr.end(),
                              [](const Cube& e){ return e.index == -2; });
	ctr.erase(new_end, ctr.end());
//	cout << ctr.size() << endl;
//	std::sort(ctr.begin(), ctr.end(), CubeComparator());
}
