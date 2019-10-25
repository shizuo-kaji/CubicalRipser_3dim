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

#include "birthday_index.h"
#include "dense_cubical_grids.h"
#include "simplex_coboundary_enumerator.h"
#include "union_find.h"
#include "write_pairs.h"
#include "joint_pairs.h"

using namespace std;

JointPairs::JointPairs(DenseCubicalGrids* _dcg, vector<BirthdayIndex>& ctr, vector<WritePairs> &_wp, const bool _print){
	dcg = _dcg;
	print = _print;
	wp = &_wp;
	ctr.clear();

	for(int x = dcg->ax - 1; x >= 0; --x){
		for(int y = dcg->ay - 1; y >= 0; --y){
			for(int z = dcg->az - 1; z >= 0; --z){
				for(int m = 2; m >= 0 ; --m){
					double birthday = dcg -> getBirthday(x,y,z,m, 1);
					if(birthday < dcg -> threshold){
						long index = dcg->getIndex(x, y, z, m);
						ctr.push_back(BirthdayIndex(birthday, index, 1));
					}
				}
			}
		}
	}
	std::sort(ctr.rbegin(), ctr.rend(), BirthdayIndexComparator());
}

void JointPairs::joint_pairs_main(vector<BirthdayIndex>& ctr){
	UnionFind dset(dcg);
	long u,v=0;
	double min_birth = dcg -> threshold;
	long min_idx;

	if(print == true){
		cout << "persistence intervals in dim " << 0 << ":" << endl;
	}
	
	for(auto &e : ctr){ 
		// we have to modify here when indexing scheme is changed
		long ind = e.index % dcg->axyz;
		u = dset.find(ind);
		switch(e.index / dcg->axyz){
			case 0:
				v = dset.find(ind+1);
				break;
			case 1:
				v = dset.find(ind+(dcg->ax));
				break;
			case 2:
				v = dset.find(ind+(dcg->axy));
				break;
		}
			
		if(u != v){
			double birth;
			int idx;
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
			double death = e.birthday;
			dset.link(u, v);
			if(birth != death){
				vector<int> loc(dcg->getXYZM(idx));
				wp -> push_back(WritePairs(0, birth, death, loc[0], loc[1], loc[2], print));
			}
			// If two values have same parent, these are potential edges which make a 2-simplex. Otherwise, remove.
	        e.dim = -2;
		}
	}
	// the based point component
	vector<int> loc(dcg->getXYZM(min_idx));
	wp -> push_back(WritePairs(0, min_birth, dcg -> threshold, loc[0],loc[1],loc[2],print));
//	cout << ctr.size() << endl;
	auto new_end = std::remove_if(ctr.begin(), ctr.end(),
                              [](const BirthdayIndex& e){ return e.dim == -2; });
	ctr.erase(new_end, ctr.end());
//	cout << ctr.size() << endl;
	std::sort(ctr.begin(), ctr.end(), BirthdayIndexComparator());
}
