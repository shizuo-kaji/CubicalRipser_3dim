/* union_find.cpp

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <vector>
#include <algorithm>
#include "birthday_index.h"
#include "dense_cubical_grids.h"
#include "union_find.h"

using namespace std;

UnionFind::UnionFind(DenseCubicalGrids* _dcg) {
	dcg = _dcg;
	n = _dcg->ax * _dcg->ay * _dcg->az;
	parent.resize(n);
	birthtime.resize(n);
	time_max.resize(n);

	for(long i = 0; i < n; ++i){
		parent[i] = i;
		vector<int> loc(_dcg->getXYZM(i));
		birthtime[i] = dcg->getBirthday(loc[0], loc[1], loc[2], loc[3], 0);
		time_max[i] = dcg ->getBirthday(loc[0], loc[1], loc[2], loc[3], 0);
	}
}

// find the root of a node x (specified by the index)
long UnionFind::find(long x){
	long y = x, z = parent[y];
	while (z != y) {
		y = z;
		z = parent[y];
	}
	// reassign parents to the found root z
	y = parent[x];
	while (z != y) {
		parent[x] = z;
		x = y;
		y = parent[x];
	}
	return z;
}

// merge nodes x and y (they should be root nodes); older will be the new parent
void UnionFind::link(long x, long y){
	if (x == y) return;
	if (birthtime[x] >= birthtime[y]){
		parent[x] = y; 
		time_max[y] = std::max(time_max[x], time_max[y]);
	} else if(birthtime[x] < birthtime[y]) {
		parent[y] = x;
		time_max[x] = std::max(time_max[x], time_max[y]);
	}
}