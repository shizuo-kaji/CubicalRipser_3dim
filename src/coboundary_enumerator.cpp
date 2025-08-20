/* coboundary_enumerator.cpp
This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji
This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <algorithm>
#include <initializer_list>
#include <vector>

#include "cube.h"
#include "dense_cubical_grids.h"
#include "coboundary_enumerator.h"

using namespace std;

CoboundaryEnumerator::CoboundaryEnumerator(DenseCubicalGrids* _dcg, uint8_t _dim)
    : dcg(_dcg), dim(_dim), nextCoface(Cube()) {}

void CoboundaryEnumerator::setCoboundaryEnumerator(Cube& _s) {
	cube = _s;
	position = 0; // current position of coface search
}

bool CoboundaryEnumerator::hasNextCoface() {
	double birth=0;
	auto cx = cube.x();
	auto cy = cube.y();
	auto cz = cube.z();
	auto cw = cube.w();

	if (dcg->dim < 4) {
		switch (dim) {
			case 0: { // dim0
				const int8_t offsets[][4] = {{0,0,0,2}, {0,0,-1,2}, {0,0,0,1}, {0,-1,0,1}, {0,0,0,0}, {-1,0,0,0}};
				for (uint8_t i = position; i < 6; ++i) {
					int x = cx + offsets[i][0];
					int y = cy + offsets[i][1];
					int z = cz + offsets[i][2];
					int m = offsets[i][3];
					birth = dcg->getBirth(x,y,z,0,m,1);
					nextCoface = Cube(birth, x,y,z, 0, m);
					if (birth != dcg->threshold) { position = i + 1; return true; }
				}
				return false;
			}

			case 1: {// dim1
				const int8_t offsets[][4][4] = {
					{{0,0,0,1}, {0,0,-1,1}, {0,0,0,0}, {0,-1,0,0}},
					{{0,0,0,2}, {0,0,-1,2}, {0,0,0,0}, {-1,0,0,0}},
					{{0,0,0,2}, {0,-1,0,2}, {0,0,0,1}, {-1,0,0,1}},
				};
				for(uint8_t i = position; i < 4; ++i){
					int x = cx + offsets[cube.m()][i][0];
					int y = cy + offsets[cube.m()][i][1];
					int z = cz + offsets[cube.m()][i][2];
					int m = offsets[cube.m()][i][3];
					birth = dcg->getBirth(x,y,z,0, m,2);
					nextCoface = Cube(birth, x,y,z, 0, m);
					if (birth != dcg->threshold) { position = i + 1; return true; }
				}
				return false;
			}

			case 2: {// dim2
				const int8_t offsets[][2][4] = {
					{{0,0,0,0}, {0,0,-1,0}},
					{{0,0,0,0}, {0,-1,0,0}},
					{{0,0,0,0}, {-1,0,0,0}},
				};
				for(uint8_t i = position; i < 2; ++i){
					int x = cx + offsets[cube.m()][i][0];
					int y = cy + offsets[cube.m()][i][1];
					int z = cz + offsets[cube.m()][i][2];
					int m = offsets[cube.m()][i][3];
					birth = dcg->getBirth(x,y,z,0, m,3);
					nextCoface = Cube(birth, x,y,z, 0, m);
					if (birth != dcg->threshold) { position = i + 1; return true; }
				}
				return false;
			}
		}
	} else { // dim == 4
		// TODO: implement
	}
	return false;
}
