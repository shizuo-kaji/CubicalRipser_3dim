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
	// current position of coface search
    if (dcg->az == 1 && dcg->config->tconstruction && dim < 2) {
        // For 2D images under T-construction, skip out-of-plane (z) cofaces
        position = 2;
    } else {
        position = 0;
    }
}

bool CoboundaryEnumerator::hasNextCoface() {
	double birth=0;
	auto cx = cube.x();
	auto cy = cube.y();
	auto cz = cube.z();
	auto cw = cube.w();

	if (dcg->dim < 4) {
		// offsets[dim][variant(m)][index][(dx,dy,dz,m')]
		static const int8_t offsets[3][3][6][4] = {
			// dim 0 (point) : 6 cofaces (duplicate rows for m-variants)
			{
				{{0,0,0,2},{0,0,-1,2},{0,0,0,1},{0,-1,0,1},{0,0,0,0},{-1,0,0,0}},
				{{0,0,0,0},{0,0, 0,0},{0,0,0,0},{0, 0,0,0},{0,0,0,0},{ 0,0,0,0}}, //pad
				{{0,0,0,0},{0,0, 0,0},{0,0,0,0},{0, 0,0,0},{0,0,0,0},{ 0,0,0,0}}, //pad
			},
			// dim 1 (edge) : 4 cofaces (pad remaining with zeros)
			{
				{{0,0,0,1},{0,0,-1,1},{0,0,0,0},{0,-1,0,0},{0,0,0,0},{0,0,0,0}},
				{{0,0,0,2},{0,0,-1,2},{0,0,0,0},{-1,0,0,0},{0,0,0,0},{0,0,0,0}},
				{{0,0,0,2},{0,-1,0,2},{0,0,0,1},{-1,0,0,1},{0,0,0,0},{0,0,0,0}},
			},
			// dim 2 (square) : 2 cofaces (pad)
			{
				{{0,0,0,0},{0,0,-1,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
				{{0,0,0,0},{0,-1,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
				{{0,0,0,0},{-1,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
			}
		};
		static const uint8_t counts[3] = {6,4,2};

		if (dim > 2) return false;
		uint8_t variant = cube.m();
		uint8_t cnt = counts[dim];

		for (uint8_t i = position; i < cnt; ++i) {
			int x = cx + offsets[dim][variant][i][0];
			int y = cy + offsets[dim][variant][i][1];
			int z = cz + offsets[dim][variant][i][2];
			int m = offsets[dim][variant][i][3];
			birth = dcg->getBirth(x,y,z,cw,m, dim+1);
			nextCoface = Cube(birth, x,y,z,cw, m);
			if (birth != dcg->threshold) { position = i + 1; return true; }
		}
		return false;
	} else { // dim == 4
		// TODO: implement
	}
	return false;
}
