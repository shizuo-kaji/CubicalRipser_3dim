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
	// note the shift of indices to account for the boundary
	auto cx = cube.x();
	auto cy = cube.y();
	auto cz = cube.z();
	const int8_t offsets0[][4] = {{1,1,2,2}, {1,1,0,2}, {1,2,1,1}, {1,0,1,1}, {2,1,1,0}, {0,1,1,0}};
	switch (dim) {
		case 0: // dim0
			for (uint8_t i = position; i < 6; ++i) {
				int x = cx + offsets0[i][0];
				int y = cy + offsets0[i][1];
				int z = cz + offsets0[i][2];
				int m = offsets0[i][3];
				
				birth = max(cube.birth, dcg->dense3[x][y][z]);
				nextCoface = Cube(birth, cx + (m == 0 ? i/4 - 1 : 0), 
										cy + (m == 1 ? i/2 - 1 : 0), 
										cz + (m == 2 ? 1 - i%2 : 0), m);

				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;

		case 1: // dim1
			switch (cube.m()) {
				case 0: { // dim1 type0 (x-axis -> )
					const int8_t offsets[][4] = {{1,1,2,1}, {1,1,0,1}, {1,2,1,0}, {1,0,1,0}};
					for(uint8_t i = position; i < 4; ++i){
						birth = max(max(cube.birth, dcg->dense3[cx+offsets[i][0]][cy+offsets[i][1]][cz+offsets[i][2]]), dcg->dense3[cx+2][cy+offsets[i][1]][cz+offsets[i][2]]);
						nextCoface = Cube(birth, cx, cy - (i==3), cz - (i==1), offsets[i][3]);
						if (birth != dcg->threshold) { position = i + 1; return true; }
					}
					return false;
				}

				case 1: { // dim1 type1 (y-axis -> )
					const int8_t offsets[][4] = {{1,1,2,2}, {1,1,0,2}, {2,1,1,0}, {0,1,1,0}};
					for(uint8_t i = position; i < 4; ++i){
						birth = max(max(cube.birth, dcg->dense3[cx+offsets[i][0]][cy+offsets[i][1]][cz+offsets[i][2]]), dcg->dense3[cx + offsets[i][0]][cy+2][cz + offsets[i][2]]);
						nextCoface = Cube(birth, cx - (i==3), cy, cz - (i==1), offsets[i][3]);
						if (birth != dcg->threshold) { position = i + 1; return true; }
					}
					return false;
				}

				case 2: { // dim1 type2 (z-axis -> )
					const int8_t offsets[][4] = {{1,2,1,2}, {1,0,1,2}, {2,1,1,1}, {0,1,1,1}};
					for(uint8_t i = position; i < 4; ++i){
						birth = max(max(cube.birth, dcg->dense3[cx+offsets[i][0]][cy+offsets[i][1]][cz+offsets[i][2]]), dcg->dense3[cx + offsets[i][0]][cy+offsets[i][1]][cz + 2]);
						nextCoface = Cube(birth, cx - (i==3), cy - (i==1), cz, offsets[i][3]);
						if (birth != dcg->threshold) { position = i + 1; return true; }
					}
					return false;
				}
			}
			return false;

		case 2: // dim2
		switch (cube.m()) {
			case 0: // dim2 type0 (fix z)
			for(uint8_t i = position; i < 2; ++i){
				switch(i){
					case 0: // upper
						birth = max(max(max(max(cube.birth, dcg->dense3[cx+1][cy+1][cz + 2]), dcg->dense3[cx + 2][cy+1][cz + 2]),
							dcg->dense3[cx+1][cy + 2][cz + 2]),dcg->dense3[cx + 2][cy + 2][cz + 2]);
						nextCoface = Cube(birth, cx, cy, cz, 0);
						break;

					case 1: // lower
						birth = max(max(max(max(cube.birth, dcg->dense3[cx+1][cy+1][cz]), dcg->dense3[cx + 2][cy+1][cz]),
							dcg->dense3[cx+1][cy + 2][cz]),dcg->dense3[cx + 2][cy + 2][cz]);
						nextCoface = Cube(birth, cx, cy, cz - 1, 0);
						break;
				}

				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;

			case 1: // dim2 type1 (fix y)
			for(uint8_t i = position; i < 2; ++i){
				switch(i){
					case 0: // left
						birth = max(max(max(max(cube.birth, dcg->dense3[cx+1][cy + 2][cz+1]), dcg->dense3[cx + 2][cy + 2][cz+1]),
							dcg->dense3[cx+1][cy + 2][cz + 2]),dcg->dense3[cx + 2][cy + 2][cz + 2]);
						nextCoface = Cube(birth, cx, cy, cz, 0);
						break;

					case 1: //right
						birth = max(max(max(max(cube.birth, dcg->dense3[cx+1][cy][cz+1]), dcg->dense3[cx + 2][cy][cz+1]),
							dcg->dense3[cx+1][cy][cz + 2]),dcg->dense3[cx + 2][cy][cz + 2]);
						nextCoface = Cube(birth, cx, cy - 1, cz, 0);
						break;
				}

				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;

			case 2: // dim2 type2 (fix x)
			for(uint8_t i = position; i < 2; ++i){
				switch(i){
				case 0: // left
					birth = max(max(max(max(cube.birth, dcg->dense3[cx + 2][cy+1][cz+1]), dcg->dense3[cx + 2][cy + 2][cz+1]),
						dcg->dense3[cx + 2][cy+1][cz + 2]),dcg->dense3[cx + 2][cy +2][cz + 2]);
					nextCoface = Cube(birth, cx, cy, cz, 0);
					break;

				case 1: //right
					birth = max(max(max(max(cube.birth, dcg->dense3[cx][cy+1][cz+1]), dcg->dense3[cx][cy + 2][cz+1]),
						dcg->dense3[cx][cy+1][cz + 2]),dcg->dense3[cx][cy + 2][cz + 2]);
					nextCoface = Cube(birth, cx - 1, cy, cz, 0);
					break;
				}

				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;
		}
	}
	return false;
}