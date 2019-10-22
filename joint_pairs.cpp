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

JointPairs::JointPairs(DenseCubicalGrids* _dcg, vector<WritePairs> &_wp, const bool _print){
	dcg = _dcg;
	print = _print;
	wp = &_wp;

	for(int x = 0; x < dcg->ax; ++x){
		for(int y = 0; y < dcg->ay; ++y){
			for(int z = 0; z < dcg->az; ++z){
				for(int m = 0; m < 3; ++m){
					double birthday = dcg -> getBirthday(x,y,z,m, 1);
					if(birthday < dcg -> threshold){
						long index = dcg->getIndex(x, y, z, m);
						dim1_simplex_list.push_back(BirthdayIndex(birthday, index, 1));
					}
				}
			}
		}
	}
	std::sort(dim1_simplex_list.rbegin(), dim1_simplex_list.rend(), BirthdayIndexComparator());
}

void JointPairs::joint_pairs_main(vector<BirthdayIndex>& ctr){
	UnionFind dset(dcg);
	ctr.clear();
	long u,v=0;
	double min_birth = dcg -> threshold;
	long min_idx,midx;

	if(print == true){
		cout << "persistence intervals in dim " << 0 << ":" << endl;
	}
	
	for(auto e : dim1_simplex_list){
		vector<int> loc(dcg->getXYZM(e.index));
		switch(loc[3]){
			case 0:
				u = dset.find(dcg->getIndex(loc[0],loc[1],loc[2]));
				v = dset.find(dcg->getIndex(loc[0]+1, loc[1], loc[2]));
				break;
			case 1:
				u = dset.find(dcg->getIndex(loc[0], loc[1], loc[2]));
				v = dset.find(dcg->getIndex(loc[0], loc[1]+1, loc[2]));
			break;
			case 2:
				u = dset.find(dcg->getIndex(loc[0], loc[1], loc[2]));
				v = dset.find(dcg->getIndex(loc[0], loc[1], loc[2]+1));
				break;
		}
			
		if(u != v){
			double birth,mbirth;
			int idx;
			if(dset.birthtime[u] >= dset.birthtime[v]){
				birth = dset.birthtime[u];
				mbirth = dset.birthtime[v];
				midx = v;
				idx = u; // the one who dies to make a cycle
			}else{
				birth = dset.birthtime[v]; 
				mbirth = dset.birthtime[u];
				idx = v; // the one who dies to make a cycle
				midx = u;
			}
			if (mbirth < min_birth) {
				min_birth = mbirth;
				min_idx = midx;
			}
			double death = e.getBirthday();
			dset.link(u, v);
			if(birth != death){
				vector<int> loc(dcg->getXYZM(idx));
				wp -> push_back(WritePairs(0, birth, death, loc[0], loc[1], loc[2], print));
			}
		} else { // If two values have same "parent", these are potential edges which make a 2-simplex.
			ctr.push_back(e);
		}
	}

	// the based point component
	vector<int> loc(dcg->getXYZM(min_idx));
	wp -> push_back(WritePairs(0, min_birth, dcg -> threshold, loc[0],loc[1],loc[2],print));
	std::sort(ctr.begin(), ctr.end(), BirthdayIndexComparator());
}
