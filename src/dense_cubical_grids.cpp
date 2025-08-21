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

// Explicit-shape constructor (V-construction)
DenseCubicalGrids::DenseCubicalGrids(Config& _config, uint8_t d, uint32_t x, uint32_t y, uint32_t z, uint32_t w)  {
    config = &_config;
    threshold = config->threshold;
    config->tconstruction = false;
    dim = d;
    ax = x; ay = y; az = z; aw = w;
    img_x = ax; img_y = ay; img_z = az; img_w = aw;
}

// return filtlation value for a cube
// (cx,cy,cz) is the voxel coordinates in the original image
double DenseCubicalGrids::getBirth(uint32_t cx, uint32_t cy, uint32_t cz){
	return (*dense3)(cx+1, cy+1, cz+1);
}

double DenseCubicalGrids::getBirth(uint32_t cx, uint32_t cy, uint32_t cz, uint32_t cw, uint8_t cm, uint8_t dim) {
	// beware of the shift due to the boundary
	if (this->dim < 4) {
		switch (dim) {
			case 0:
				return (*dense3)(cx+1, cy+1, cz+1);
			case 1: {
				static const int off[13][3] = {
					{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,-1,0},
					{0,-1,1},{0,1,1},{1,-1,1},{1,0,1},{1,1,1},
					{1,-1,-1},{1,0,-1},{1,1,-1}
				};
				if (cm < 13) {
					const double b = (*dense3)(cx+1, cy+1, cz+1);
					const int *o = off[cm];
					return max(b, (*dense3)(cx+1+o[0], cy+1+o[1], cz+1+o[2]));
				}
				// fallthrough on invalid cm
			}
			case 2:
				switch (cm) {
				case 0: // x - y (fix z)
					return max({ (*dense3)(cx+1, cy+1, cz+1), (*dense3)(cx+2, cy+1, cz+1),
						(*dense3)(cx+2, cy+2, cz+1), (*dense3)(cx+1, cy+2, cz+1) });
				case 1: // z - x (fix y)
					return max({ (*dense3)(cx+1, cy+1, cz+1), (*dense3)(cx+1, cy+1, cz+2),
						(*dense3)(cx+2, cy+1, cz+2), (*dense3)(cx+2, cy+1, cz+1) });
				case 2: // y - z (fix x)
					return max({ (*dense3)(cx+1, cy+1, cz+1), (*dense3)(cx+1, cy+2, cz+1),
						(*dense3)(cx+1, cy+2, cz+2), (*dense3)(cx+1, cy+1, cz+2) });
				}
			case 3:
				return max({ (*dense3)(cx+1, cy+1, cz+1), (*dense3)(cx+2, cy+1, cz+1),
					(*dense3)(cx+2, cy+2, cz+1), (*dense3)(cx+1, cy+2, cz+1),
					(*dense3)(cx+1, cy+1, cz+2), (*dense3)(cx+2, cy+1, cz+2),
					(*dense3)(cx+2, cy+2, cz+2), (*dense3)(cx+1, cy+2, cz+2) });
			}
	}
	return threshold;
}


// (x,y,z) or (x,y,z,w) of the voxel which defines the birthtime of the cube
vector<uint32_t> DenseCubicalGrids::ParentVoxel(uint8_t _dim, Cube &c){
	(void)_dim;
	uint32_t cx = c.x(), cy = c.y(), cz = c.z(), cw = c.w();

	if (dim < 4) {
		static const int rel[][3] = {
			{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},
			{1,-1,0},{0,-1,1},{1,-1,1},{1,-1,-1},{1,0,-1},{1,1,-1}
		};
		for (auto &r : rel) {
			int dx=r[0], dy=r[1], dz=r[2];
			if (c.birth == (*dense3)(cx+1+dx, cy+1+dy, cz+1+dz))
				return {cx+dx, cy+dy, cz+dz};
		}
		cerr << "parent voxel not found!" << endl;
		return {0,0,0};
	}
}
