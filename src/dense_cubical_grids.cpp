/* dense_cubical_grids.cpp

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
#include <string>
#include <initializer_list>

// switch V and T constructions
#include "dense_cubical_grids.h"

using namespace std;


DenseCubicalGrids::DenseCubicalGrids(Config& _config)  {
	config = &_config;
	threshold = config->threshold;
	config->tconstruction = false;
}

// return filtlation value for a cube
// (cx,cy,cz) is the voxel coordinates in the original image
double DenseCubicalGrids::getBirth(uint32_t cx, uint32_t cy, uint32_t cz){
	return dense3[cx+1][cy+1][cz+1];
}

double DenseCubicalGrids::getBirth(uint32_t cx, uint32_t cy, uint32_t cz, uint32_t cw, uint8_t cm, uint8_t dim) {
	// beware of the shift due to the boundary
	if (this->dim < 4) {
		switch (dim) {
			case 0:
				return dense3[cx+1][cy+1][cz+1];
			case 1:
				switch (cm) {
				case 0:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 1][cz + 1]);
				case 1:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 1][cy + 2][cz + 1]);
				case 2:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 1][cy + 1][cz + 2]);
				case 3:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 2][cz + 1]);
				case 4:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 0][cz + 1]);
				// for 3d dual only
				case 5:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 1][cy + 0][cz + 2]);
				case 6:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 1][cy + 2][cz + 2]);
				case 7:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 0][cz + 2]);
				case 8:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 1][cz + 2]);
				case 9:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 2][cz + 2]);
				case 10:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 0][cz + 0]);
				case 11:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 1][cz + 0]);
				case 12:
					return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 2][cz + 0]);
				}
			case 2:
				switch (cm) {
				case 0: // x - y (fix z)
					return max({ dense3[cx+1][cy+1][cz+1], dense3[cx+2][cy+1][cz+1],
						dense3[cx+2][cy+2][cz+1], dense3[cx+1][cy+2][cz+1] });
				case 1: // z - x (fix y)
					return max({ dense3[cx+1][cy+1][cz+1], dense3[cx+1][cy+1][cz+2],
						dense3[cx+2][cy+1][cz+2], dense3[cx+2][cy+1][cz+1] });
				case 2: // y - z (fix x)
					return max({ dense3[cx+1][cy+1][cz+1], dense3[cx+1][cy+2][cz+1],
						dense3[cx+1][cy+2][cz+2], dense3[cx+1][cy+1][cz+2] });
				}
			case 3:
				return max({ dense3[cx+1][cy+1][cz+1], dense3[cx+2][cy+1][cz+1],
					dense3[cx+2][cy+2][cz+1], dense3[cx+1][cy+2][cz+1],
					dense3[cx+1][cy+1][cz+2], dense3[cx+2][cy+1][cz+2],
					dense3[cx+2][cy+2][cz+2], dense3[cx+1][cy+2][cz+2] });
			}
	} else { // dim == 4
		switch (dim) {
			case 0:
				return dense4[cx+1][cy+1][cz+1][cw+1];
			case 1:
				switch (cm) {
				case 0: // x
					return max(dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+2][cy+1][cz+1][cw+1]);
				case 1: // y
					return max(dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+1][cy+2][cz+1][cw+1]);
				case 2: // z
					return max(dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+1][cy+1][cz+2][cw+1]);
				case 3: // w
					return max(dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+1][cy+1][cz+1][cw+2]);
				}
			case 2:
				switch (cm) {
				case 0: // xy
					return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+2][cy+1][cz+1][cw+1], dense4[cx+1][cy+2][cz+1][cw+1], dense4[cx+2][cy+2][cz+1][cw+1] });
				case 1: // xz
					return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+2][cy+1][cz+1][cw+1], dense4[cx+1][cy+1][cz+2][cw+1], dense4[cx+2][cy+1][cz+2][cw+1] });
				case 2: // xw
					return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+2][cy+1][cz+1][cw+1], dense4[cx+1][cy+1][cz+1][cw+2], dense4[cx+2][cy+1][cz+1][cw+2] });
				case 3: // yz
					return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+1][cy+2][cz+1][cw+1], dense4[cx+1][cy+1][cz+2][cw+1], dense4[cx+1][cy+2][cz+2][cw+1] });
				case 4: // yw
					return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+1][cy+2][cz+1][cw+1], dense4[cx+1][cy+1][cz+1][cw+2], dense4[cx+1][cy+2][cz+1][cw+2] });
				case 5: // zw
					return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+1][cy+1][cz+2][cw+1], dense4[cx+1][cy+1][cz+1][cw+2], dense4[cx+1][cy+1][cz+2][cw+2] });
				}
			case 3:
				switch (cm) {
				case 0: // xyz
					return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+2][cy+1][cz+1][cw+1], dense4[cx+1][cy+2][cz+1][cw+1], dense4[cx+2][cy+2][cz+1][cw+1],
						dense4[cx+1][cy+1][cz+2][cw+1], dense4[cx+2][cy+1][cz+2][cw+1], dense4[cx+1][cy+2][cz+2][cw+1], dense4[cx+2][cy+2][cz+2][cw+1] });
				case 1: // xyw
					return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+2][cy+1][cz+1][cw+1], dense4[cx+1][cy+2][cz+1][cw+1], dense4[cx+2][cy+2][cz+1][cw+1],
						dense4[cx+1][cy+1][cz+1][cw+2], dense4[cx+2][cy+1][cz+1][cw+2], dense4[cx+1][cy+2][cz+1][cw+2], dense4[cx+2][cy+2][cz+1][cw+2] });
				case 2: // xzw
					return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+2][cy+1][cz+1][cw+1], dense4[cx+1][cy+1][cz+2][cw+1], dense4[cx+2][cy+1][cz+2][cw+1],
						dense4[cx+1][cy+1][cz+1][cw+2], dense4[cx+2][cy+1][cz+1][cw+2], dense4[cx+1][cy+1][cz+2][cw+2], dense4[cx+2][cy+1][cz+2][cw+2] });
				case 3: // yzw
					return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+1][cy+2][cz+1][cw+1], dense4[cx+1][cy+1][cz+2][cw+1], dense4[cx+1][cy+2][cz+2][cw+1],
						dense4[cx+1][cy+1][cz+1][cw+2], dense4[cx+1][cy+2][cz+1][cw+2], dense4[cx+1][cy+1][cz+2][cw+2], dense4[cx+1][cy+2][cz+2][cw+2] });
				}
			case 4:
				return max({ dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+2][cy+1][cz+1][cw+1], dense4[cx+1][cy+2][cz+1][cw+1], dense4[cx+2][cy+2][cz+1][cw+1],
					dense4[cx+1][cy+1][cz+2][cw+1], dense4[cx+2][cy+1][cz+2][cw+1], dense4[cx+1][cy+2][cz+2][cw+1], dense4[cx+2][cy+2][cz+2][cw+1],
					dense4[cx+1][cy+1][cz+1][cw+2], dense4[cx+2][cy+1][cz+1][cw+2], dense4[cx+1][cy+2][cz+1][cw+2], dense4[cx+2][cy+2][cz+1][cw+2],
					dense4[cx+1][cy+1][cz+2][cw+2], dense4[cx+2][cy+1][cz+2][cw+2], dense4[cx+1][cy+2][cz+2][cw+2], dense4[cx+2][cy+2][cz+2][cw+2] });
		}
	}
	return threshold;
}


// (x,y,z) of the voxel which defines the birthtime of the cube
vector<uint32_t> DenseCubicalGrids::ParentVoxel(uint8_t _dim, Cube &c){
	uint32_t cx = c.x();
	uint32_t cy = c.y();
	uint32_t cz = c.z();
	uint32_t cw = c.w();

	if (dim < 4) {
		if(c.birth == dense3[cx+1][cy+1][cz+1]){
			return {cx,cy,cz};
		}else if(c.birth == dense3[cx+2][cy+1][cz+1]){
			return {cx+1,cy,cz};
		}else if(c.birth == dense3[cx+2][cy+2][cz+1]){
			return {cx+1,cy+1,cz};
		}else if(c.birth == dense3[cx+1][cy+2][cz+1]){
			return {cx,cy+1,cz};
		}else if(c.birth == dense3[cx+1][cy+1][cz+2]){
			return {cx,cy,cz+1};
		}else if(c.birth == dense3[cx+2][cy+1][cz+2]){
			return {cx+1,cy,cz+1};
		}else if(c.birth == dense3[cx+1][cy+2][cz+2]){
			return {cx,cy+1,cz+1};
		}else if(c.birth == dense3[cx+2][cy+2][cz+2]){
			return {cx+1,cy+1,cz+1};
		}else if(c.birth == dense3[cx+2][cy+0][cz+1]){			// for 3d dual only
			return { cx + 1,cy - 1,cz };
		}else if(c.birth == dense3[cx+1][cy+0][cz+2]){
			return { cx,cy - 1,cz + 1 };
		}else if(c.birth == dense3[cx+2][cy+0][cz+2]){
			return { cx + 1,cy - 1,cz + 1 };
		}else if(c.birth == dense3[cx+2][cy+0][cz+0]){
			return { cx + 1,cy - 1,cz - 1 };
		}else if(c.birth == dense3[cx+2][cy+1][cz+0]){
			return { cx + 1,cy,cz - 1 };
		}else if(c.birth == dense3[cx+2][cy+2][cz+0]){
			return { cx + 1,cy + 1,cz - 1 };
		}else{
			cerr << "parent voxel not found!" << endl;
			return { 0,0,0 };
		}
	} else { // dim == 4
		if(c.birth == dense4[cx+1][cy+1][cz+1][cw+1]){
			return {cx,cy,cz,cw};
		} else if(c.birth == dense4[cx+2][cy+1][cz+1][cw+1]){
			return {cx+1,cy,cz,cw};
		} else if(c.birth == dense4[cx+1][cy+2][cz+1][cw+1]){
			return {cx,cy+1,cz,cw};
		} else if(c.birth == dense4[cx+2][cy+2][cz+1][cw+1]){
			return {cx+1,cy+1,cz,cw};
		} else if(c.birth == dense4[cx+1][cy+1][cz+2][cw+1]){
			return {cx,cy,cz+1,cw};
		} else if(c.birth == dense4[cx+2][cy+1][cz+2][cw+1]){
			return {cx+1,cy,cz+1,cw};
		} else if(c.birth == dense4[cx+1][cy+2][cz+2][cw+1]){
			return {cx,cy+1,cz+1,cw};
		} else if(c.birth == dense4[cx+2][cy+2][cz+2][cw+1]){
			return {cx+1,cy+1,cz+1,cw};
		} else if(c.birth == dense4[cx+1][cy+1][cz+1][cw+2]){
			return {cx,cy,cz,cw+1};
		} else if(c.birth == dense4[cx+2][cy+1][cz+1][cw+2]){
			return {cx+1,cy,cz,cw+1};
		} else if(c.birth == dense4[cx+1][cy+2][cz+1][cw+2]){
			return {cx,cy+1,cz,cw+1};
		} else if(c.birth == dense4[cx+2][cy+2][cz+1][cw+2]){
			return {cx+1,cy+1,cz,cw+1};
		} else if(c.birth == dense4[cx+1][cy+1][cz+2][cw+2]){
			return {cx,cy,cz+1,cw+1};
		} else if(c.birth == dense4[cx+2][cy+1][cz+2][cw+2]){
			return {cx+1,cy,cz+1,cw+1};
		} else if(c.birth == dense4[cx+1][cy+2][cz+2][cw+2]){
			return {cx,cy+1,cz+1,cw+1};
		} else if(c.birth == dense4[cx+2][cy+2][cz+2][cw+2]){
			return {cx+1,cy+1,cz+1,cw+1};
		} else {
			cerr << "parent voxel not found!" << endl;
			return {0,0,0,0};
		}
	}
}
