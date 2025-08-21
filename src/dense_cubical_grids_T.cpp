/* dense_cubical_grids_T.cpp (for T-construction)

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

// Explicit-shape constructor (T-construction)
DenseCubicalGrids::DenseCubicalGrids(Config& _config, uint8_t d, uint32_t x, uint32_t y, uint32_t z, uint32_t w)  {
    config = &_config;
    threshold = config->threshold;
    config->tconstruction = true;
    dim = d;
    ax = x; ay = y; az = z; aw = w;
    img_x = ax; img_y = ay; img_z = az; img_w = aw;
}

// return filtlation value for a cube
double DenseCubicalGrids::getBirth(uint32_t cx, uint32_t cy, uint32_t cz){
	return min({ (*dense3)(cx, cy, cz), (*dense3)(cx+1, cy, cz),
				(*dense3)(cx+1, cy+1, cz), (*dense3)(cx, cy+1, cz),
				(*dense3)(cx, cy, cz+1), (*dense3)(cx+1, cy, cz+1),
				(*dense3)(cx+1, cy+1, cz+1), (*dense3)(cx, cy+1, cz+1) });
}
double DenseCubicalGrids::getBirth(uint32_t cx, uint32_t cy, uint32_t cz, uint32_t cw, uint8_t cm, uint8_t dim) {
	// beware of the shift due to the boundary
	if (this->dim < 4) {
		switch (dim) {
			case 3:
				return (*dense3)(cx+1, cy+1, cz+1);
			case 2:
				switch (cm) {
					case 0: // fix z
						return min((*dense3)(cx+1, cy+1, cz), (*dense3)(cx + 1, cy + 1, cz + 1));
					case 1: // fix y
						return min((*dense3)(cx+1, cy, cz+1), (*dense3)(cx + 1, cy + 1, cz + 1));
					case 2: // fix x
						return min((*dense3)(cx, cy+1, cz+1), (*dense3)(cx + 1, cy + 1, cz + 1));
					break;
				}
			case 1:
				switch (cm) {
					case 0: // x,x+1
						return min({ (*dense3)(cx+1, cy+1, cz+1), (*dense3)(cx+1, cy+1, cz),
							(*dense3)(cx+1, cy, cz+1), (*dense3)(cx+1, cy, cz) });
					case 1: // y,y+1
						return min({ (*dense3)(cx+1, cy+1, cz+1), (*dense3)(cx, cy+1, cz+1),
							(*dense3)(cx+1, cy+1, cz), (*dense3)(cx, cy+1, cz) });
					case 2: // z,z+1
						return min({ (*dense3)(cx+1, cy+1, cz+1), (*dense3)(cx, cy+1, cz+1),
							(*dense3)(cx+1, cy, cz+1), (*dense3)(cx, cy, cz+1) });
					break;
				}
			case 0:
				return min({ (*dense3)(cx, cy, cz), (*dense3)(cx+1, cy, cz),
					(*dense3)(cx+1, cy+1, cz), (*dense3)(cx, cy+1, cz),
					(*dense3)(cx, cy, cz+1), (*dense3)(cx+1, cy, cz+1),
					(*dense3)(cx+1, cy+1, cz+1), (*dense3)(cx, cy+1, cz+1) });
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
		static const int off[8][3] = {
			{1,1,1},{0,1,1},{0,0,1},{0,0,0},{0,1,0},{1,0,1},{1,0,0},{1,1,0}
		};
		for (const auto &o : off) {
			if (c.birth == (*dense3)(cx+o[0], cy+o[1], cz+o[2]))
				return { uint32_t(cx+o[0]-1), uint32_t(cy+o[1]-1), uint32_t(cz+o[2]-1) };
		}
		cerr << "parent voxel not found!" << endl;
		return {0,0,0};
	}
}
