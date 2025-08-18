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
	config->tconstruction = true;
}


// return filtlation value for a cube
double DenseCubicalGrids::getBirth(uint32_t cx, uint32_t cy, uint32_t cz){
	return min({ dense3[cx][cy][cz], dense3[cx+1][cy][cz],
				dense3[cx+1][cy+1][cz], dense3[cx][cy+1][cz],
				dense3[cx][cy][cz+1], dense3[cx+1][cy][cz+1],
				dense3[cx+1][cy+1][cz+1], dense3[cx][cy+1][cz+1] });
}
double DenseCubicalGrids::getBirth(uint32_t cx, uint32_t cy, uint32_t cz, uint32_t cw, uint8_t cm, uint8_t dim) {
	// beware of the shift due to the boundary
	if (this->dim < 4) {
		switch (dim) {
			case 3:
				return dense3[cx+1][cy+1][cz+1];
			case 2:
				switch (cm) {
					case 0: // fix z
						return min(dense3[cx+1][cy+1][cz], dense3[cx + 1][cy + 1][cz + 1]);
					case 1: // fix y
						return min(dense3[cx+1][cy][cz+1], dense3[cx + 1][cy + 1][cz + 1]);
					case 2: // fix x
						return min(dense3[cx][cy+1][cz+1], dense3[cx + 1][cy + 1][cz + 1]);
					break;
				}
			case 1:
				switch (cm) {
					case 0: // x,x+1
						return min({ dense3[cx+1][cy+1][cz+1], dense3[cx+1][cy+1][cz],
							dense3[cx+1][cy][cz+1], dense3[cx+1][cy][cz] });
					case 1: // y,y+1
						return min({ dense3[cx+1][cy+1][cz+1], dense3[cx][cy+1][cz+1],
							dense3[cx+1][cy+1][cz], dense3[cx][cy+1][cz] });
					case 2: // z,z+1
						return min({ dense3[cx+1][cy+1][cz+1], dense3[cx][cy+1][cz+1],
							dense3[cx+1][cy][cz+1], dense3[cx][cy][cz+1] });
					break;
				}
			case 0:
				return min({ dense3[cx][cy][cz], dense3[cx+1][cy][cz],
					dense3[cx+1][cy+1][cz], dense3[cx][cy+1][cz],
					dense3[cx][cy][cz+1], dense3[cx+1][cy][cz+1],
					dense3[cx+1][cy+1][cz+1], dense3[cx][cy+1][cz+1] });
			}
	} else { // dim == 4
		switch (dim) {
			case 4:
				return dense4[cx+1][cy+1][cz+1][cw+1];
			case 3:
				switch (cm) {
					case 0: // fix w
						return min(dense4[cx+1][cy+1][cz+1][cw], dense4[cx+1][cy+1][cz+1][cw+1]);
					case 1: // fix z
						return min(dense4[cx+1][cy+1][cz][cw+1], dense4[cx+1][cy+1][cz+1][cw+1]);
					case 2: // fix y
						return min(dense4[cx+1][cy][cz+1][cw+1], dense4[cx+1][cy+1][cz+1][cw+1]);
					case 3: // fix x
						return min(dense4[cx][cy+1][cz+1][cw+1], dense4[cx+1][cy+1][cz+1][cw+1]);
				}
			case 2:
				switch (cm) {
					case 0: // yz
						return min({dense4[cx+1][cy+1][cz+1][cw], dense4[cx+1][cy+1][cz][cw], dense4[cx+1][cy][cz+1][cw], dense4[cx+1][cy][cz][cw]});
					case 1: // xz
						return min({dense4[cx+1][cy+1][cz+1][cw], dense4[cx+1][cy][cz+1][cw], dense4[cx][cy+1][cz+1][cw], dense4[cx][cy][cz+1][cw]});
					case 2: // xy
						return min({dense4[cx+1][cy+1][cz][cw], dense4[cx+1][cy][cz][cw], dense4[cx][cy+1][cz][cw], dense4[cx][cy][cz][cw]});
					case 3: // yw
						return min({dense4[cx+1][cy+1][cz+1][cw], dense4[cx+1][cy+1][cz][cw], dense4[cx+1][cy][cz+1][cw], dense4[cx+1][cy][cz][cw]});
					case 4: // xw
						return min({dense4[cx+1][cy+1][cz+1][cw], dense4[cx+1][cy][cz+1][cw], dense4[cx][cy+1][cz+1][cw], dense4[cx][cy][cz+1][cw]});
					case 5: // zw
						return min({dense4[cx+1][cy+1][cz][cw], dense4[cx+1][cy][cz][cw], dense4[cx][cy+1][cz][cw], dense4[cx][cy][cz][cw]});
				}
			case 1:
				return min({dense4[cx+1][cy+1][cz+1][cw], dense4[cx+1][cy+1][cz][cw], dense4[cx+1][cy][cz+1][cw], dense4[cx+1][cy][cz][cw],
				 dense4[cx][cy+1][cz+1][cw], dense4[cx][cy+1][cz][cw], dense4[cx][cy][cz+1][cw], dense4[cx][cy][cz][cw]});
			case 0:
				return min({dense4[cx+1][cy+1][cz+1][cw+1], dense4[cx+1][cy+1][cz][cw+1], dense4[cx+1][cy][cz+1][cw+1], dense4[cx+1][cy][cz][cw+1],
				 dense4[cx][cy+1][cz+1][cw+1], dense4[cx][cy+1][cz][cw+1], dense4[cx][cy][cz+1][cw+1], dense4[cx][cy][cz][cw+1],
				 dense4[cx+1][cy+1][cz+1][cw], dense4[cx+1][cy+1][cz][cw], dense4[cx+1][cy][cz+1][cw], dense4[cx+1][cy][cz][cw],
				 dense4[cx][cy+1][cz+1][cw], dense4[cx][cy+1][cz][cw], dense4[cx][cy][cz+1][cw], dense4[cx][cy][cz][cw]});
		}
	}
	return threshold; // dim > 4
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
		}else if(c.birth == dense3[cx+0][cy+1][cz+1]){			// T-construction
			return { cx - 1,cy,cz };
		}else if(c.birth == dense3[cx+0][cy+0][cz+1]){
			return { cx - 1,cy-1,cz };
		}else if(c.birth == dense3[cx+0][cy+0][cz+0]){
			return { cx - 1,cy-1,cz-1 };
		}else if(c.birth == dense3[cx+0][cy+1][cz+0]){
			return { cx - 1,cy,cz-1 };
		}else if(c.birth == dense3[cx+1][cy+0][cz+1]){
			return { cx,cy-1,cz};
		}else if(c.birth == dense3[cx+1][cy+0][cz+0]){
			return { cx,cy-1,cz-1};
		}else if(c.birth == dense3[cx+1][cy+1][cz+0]){
			return { cx,cy,cz-1};
		}else{
			cerr << "parent voxel not found!" << endl;
			return { 0,0,0 };
		}
	} else { // dim == 4
		if(c.birth == dense4[cx+1][cy+1][cz+1][cw+1]){
			return {cx,cy,cz,cw};
		} else if(c.birth == dense4[cx+0][cy+1][cz+1][cw+1]){
			return {cx-1,cy,cz,cw};
		} else if(c.birth == dense4[cx+0][cy+0][cz+1][cw+1]){
			return {cx-1,cy-1,cz,cw};
		} else if(c.birth == dense4[cx+0][cy+0][cz+0][cw+1]){
			return {cx-1,cy-1,cz-1,cw};
		} else if(c.birth == dense4[cx+0][cy+1][cz+0][cw+1]){
			return {cx-1,cy,cz-1,cw};
		} else if(c.birth == dense4[cx+1][cy+0][cz+1][cw+1]){
			return {cx,cy-1,cz,cw};
		} else if(c.birth == dense4[cx+1][cy+0][cz+0][cw+1]){
			return {cx,cy-1,cz-1,cw};
		} else if(c.birth == dense4[cx+1][cy+1][cz+0][cw+1]){
			return {cx,cy,cz-1,cw};
		} else if(c.birth == dense4[cx+1][cy+1][cz+1][cw+0]){
			return {cx,cy,cz,cw-1};
		} else if(c.birth == dense4[cx+0][cy+1][cz+1][cw+0]){
			return {cx-1,cy,cz,cw-1};
		} else if(c.birth == dense4[cx+0][cy+0][cz+1][cw+0]){
			return {cx-1,cy-1,cz,cw-1};
		} else if(c.birth == dense4[cx+0][cy+0][cz+0][cw+0]){
			return {cx-1,cy-1,cz-1,cw-1};
		} else if(c.birth == dense4[cx+0][cy+1][cz+0][cw+0]){
			return {cx-1,cy,cz-1,cw-1};
		} else if(c.birth == dense4[cx+1][cy+0][cz+1][cw+0]){
			return {cx,cy-1,cz,cw-1};
		} else if(c.birth == dense4[cx+1][cy+0][cz+0][cw+0]){
			return {cx,cy-1,cz-1,cw-1};
		} else if(c.birth == dense4[cx+1][cy+1][cz+0][cw+0]){
			return {cx,cy,cz-1,cw-1};
		} else {
			cerr << "parent voxel not found!" << endl;
			return {0,0,0,0};
		}
	}
}