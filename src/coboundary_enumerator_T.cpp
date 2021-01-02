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

CoboundaryEnumerator::CoboundaryEnumerator(DenseCubicalGrids* _dcg, uint8_t _dim){
	nextCoface = Cube();
	dcg = _dcg;
    dim = _dim;
}
	
void CoboundaryEnumerator::setCoboundaryEnumerator(Cube& _s) {
	cube = _s;
    if(dcg->az==1 && dim ==1){
        position = 2; // for 2d image, we can skip "out-of-the-plane" cofaces
    }else{
        position = 0; // current position of coface search
    }
}

bool CoboundaryEnumerator::hasNextCoface() {
	double birth=0;
	// note the shift of indices to account for the boundary
	auto cx = cube.x();
	auto cy = cube.y();
	auto cz = cube.z();
	switch (dim) {
		case 0: // dim0
		for (uint8_t i = position; i < 6; ++i) {
			switch (i){
			case 0:
				birth = dcg->getBirth(cx,cy,cz,2,1);
				nextCoface = Cube(birth, cx, cy, cz, 2);
				break;

			case 1:
				birth = dcg->getBirth(cx,cy,cz-1,2,1);
				nextCoface = Cube(birth, cx, cy, cz - 1, 2);
				break;

			case 2:
				birth = dcg->getBirth(cx,cy,cz,1,1);
				nextCoface = Cube(birth, cx, cy, cz, 1);
				break;

			case 3:
				birth = dcg->getBirth(cx,cy-1,cz,1,1);
				nextCoface = Cube(birth, cx, cy - 1, cz, 1);
				break;

			case 4:
				birth = dcg->getBirth(cx,cy,cz,0,1);
				nextCoface = Cube(birth, cx, cy, cz, 0);
				break;

			case 5:
				birth = dcg->getBirth(cx-1,cy,cz,0,1);
				nextCoface = Cube(birth, cx - 1, cy, cz, 0);
				break;
			}

			if (birth != dcg->threshold) {
				position = i + 1;
				return true;
			}
		}
		return false;

		case 1: // dim1
		switch (cube.m()) {
			case 0: // dim1 type0 (x-axis -> )
			for(uint8_t i = position; i < 4; ++i){
				switch(i){
				case 0:
					birth = dcg->getBirth(cx,cy,cz,1,2);
					nextCoface = Cube(birth, cx, cy, cz, 1);
					break;

				case 1:
					birth = dcg->getBirth(cx,cy,cz-1,1,2);
					nextCoface = Cube(birth, cx, cy, cz - 1, 1);
					break;

				case 2:
					birth = dcg->getBirth(cx,cy,cz,0,2);
					nextCoface = Cube(birth, cx, cy, cz, 0);
					break;

				case 3:
					birth = dcg->getBirth(cx,cy-1,cz,0,2);
					nextCoface = Cube(birth, cx, cy - 1, cz, 0);
					break;
				}
				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;

			case 1: // dim1 type1 (y-axis -> )
			for(uint8_t i = position; i < 4; ++i){
				switch(i){
				case 0:
					birth = dcg->getBirth(cx,cy,cz,2,2);
					nextCoface = Cube(birth, cx, cy, cz, 2);
					break;

				case 1:
					birth = dcg->getBirth(cx,cy,cz-1,2,2);
					nextCoface = Cube(birth, cx, cy, cz - 1, 2);
					break;

				case 2:
					birth = dcg->getBirth(cx,cy,cz,0,2);
					nextCoface = Cube(birth, cx, cy, cz, 0);
					break;

				case 3:
					birth = dcg->getBirth(cx-1,cy,cz,0,2);
					nextCoface = Cube(birth, cx - 1, cy, cz, 0);
					break;
				}

				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;

			case 2: // dim1 type2 (z-axis -> )
			for(uint8_t i = position; i < 4; ++i){
				switch(i){
					case 0:
						birth = dcg->getBirth(cx,cy,cz,2,2);
						nextCoface = Cube(birth, cx, cy, cz, 2);
						break;

					case 1:
						birth = dcg->getBirth(cx,cy-1,cz,2,2);
						nextCoface = Cube(birth, cx, cy - 1, cz, 2);
						break;

					case 2:
						birth = dcg->getBirth(cx,cy,cz,1,2);
						nextCoface = Cube(birth, cx, cy, cz, 1);
						break;

					case 3:
						birth = dcg->getBirth(cx-1,cy,cz,1,2);
						nextCoface = Cube(birth, cx - 1, cy, cz, 1);
						break;
				}

				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;
		}
		return false;

		case 2: // dim2
		switch (cube.m()) {
			case 0: // dim2 type0 (fix z)
			for(uint8_t i = position; i < 2; ++i){
				switch(i){
					case 0: // upper
						birth = dcg->getBirth(cx,cy,cz,0,3);
						nextCoface = Cube(birth, cx, cy, cz, 0);
						break;

					case 1: // lower
						birth = dcg->getBirth(cx,cy,cz-1,0,3);
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
						birth = dcg->getBirth(cx,cy,cz,0,3);
						nextCoface = Cube(birth, cx, cy, cz, 0);
						break;

					case 1: //right
						birth = dcg->getBirth(cx,cy-1,cz,0,3);
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
					birth = dcg->getBirth(cx,cy,cz,0,3);
					nextCoface = Cube(birth, cx, cy, cz, 0);
					break;

				case 1: //right
					birth = dcg->getBirth(cx-1,cy,cz,0,3);
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
