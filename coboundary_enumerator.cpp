/* cube_coboundary_enumerator.cpp

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
#include <vector>

#include "cube.h"
#include "dense_cubical_grids.h"
#include "coboundary_enumerator.h"

using namespace std;

CoboundaryEnumerator::CoboundaryEnumerator(DenseCubicalGrids* _dcg, int _dim){
	nextCoface = Cube();
	dcg = _dcg;
    dim = _dim;
}
	
void CoboundaryEnumerator::setCoboundaryEnumerator(Cube& _s) {
	cube = _s;
	position = 0; // current position of coface search
}

bool CoboundaryEnumerator::hasNextCoface() {
	double birth=0;
	// note the shift for the boundary
	switch (dim) {
		case 0: // dim0
		for (int i = position; i < 6; ++i) {
			switch (i){
			case 0:
				birth = max(cube.birth, dcg->dense3[cube.x()+1][cube.y()+1][cube.z()+2]);
				nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 2);
				break;

			case 1:
				birth = max(cube.birth, dcg->dense3[cube.x()+1][cube.y()+1][cube.z()]);
				nextCoface = Cube(birth, cube.x(), cube.y(), cube.z() - 1, 2);
				break;

			case 2:
				birth = max(cube.birth, dcg->dense3[cube.x()+1][cube.y() + 2][cube.z()+1]);
				nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 1);
				break;

			case 3:
				birth = max(cube.birth, dcg->dense3[cube.x()+1][cube.y()][cube.z()+1]);
				nextCoface = Cube(birth, cube.x(), cube.y() - 1, cube.z(), 1);
				break;

			case 4:
				birth = max(cube.birth, dcg->dense3[cube.x() + 2][cube.y()+1][cube.z()+1]);
				nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 0);
				break;

			case 5:
				birth = max(cube.birth, dcg->dense3[cube.x()][cube.y()+1][cube.z()+1]);
				nextCoface = Cube(birth, cube.x() - 1, cube.y(), cube.z(), 0);
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
			for(int i = position; i < 4; ++i){
				switch(i){
				case 0:
					birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y()+1][cube.z() + 2], dcg->dense3[cube.x() + 2][cube.y()+1][cube.z() + 2] });
					nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 1);
					break;

				case 1:
					birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y()+1][cube.z()], dcg->dense3[cube.x() + 2][cube.y()+1][cube.z()] });
					nextCoface = Cube(birth, cube.x(), cube.y(), cube.z() - 1, 1);
					break;

				case 2:
					birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y() + 2][cube.z()+1], dcg->dense3[cube.x() + 2][cube.y() + 2][cube.z()+1] });
					nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 0);
					break;

				case 3:
					birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y()][cube.z()+1], dcg->dense3[cube.x() + 2][cube.y()][cube.z()+1] });
					nextCoface = Cube(birth, cube.x(), cube.y() - 1, cube.z(), 0);
					break;
				}

				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;

			case 1: // dim1 type1 (y-axis -> )
			for(int i = position; i < 4; ++i){
				switch(i){
				case 0:
					birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y()+1][cube.z() + 2], dcg->dense3[cube.x()+1][cube.y() + 2][cube.z() + 2] });
					nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 2);
					break;

				case 1:
					birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y()+1][cube.z()], dcg->dense3[cube.x()+1][cube.y() + 2][cube.z()] });
					nextCoface = Cube(birth, cube.x(), cube.y(), cube.z() - 1, 2);
					break;

				case 2:
					birth = max({ cube.birth, dcg->dense3[cube.x() + 2][cube.y()+1][cube.z()+1], dcg->dense3[cube.x() + 2][cube.y() + 2][cube.z()+1] });
					nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 0);
					break;

				case 3:
					birth = max({ cube.birth, dcg->dense3[cube.x()][cube.y()+1][cube.z()+1], dcg->dense3[cube.x()][cube.y() + 2][cube.z()+1] });
					nextCoface = Cube(birth, cube.x() - 1, cube.y(), cube.z(), 0);
					break;
				}

				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;

			case 2: // dim1 type2 (z-axis -> )
			for(int i = position; i < 4; ++i){
				switch(i){
					case 0:
						birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y() + 2][cube.z()+1], dcg->dense3[cube.x()+1][cube.y() + 2][cube.z() + 2] });
						nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 2);
						break;

					case 1:
						birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y()][cube.z()+1], dcg->dense3[cube.x()+1][cube.y()][cube.z() + 2] });
						nextCoface = Cube(birth, cube.x(), cube.y() - 1, cube.z(), 2);
						break;

					case 2:
						birth = max({ cube.birth, dcg->dense3[cube.x() + 2][cube.y()+1][cube.z()+1], dcg->dense3[cube.x() + 2][cube.y()+1][cube.z() + 2] });
						nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 1);
						break;

					case 3:
						birth = max({ cube.birth, dcg->dense3[cube.x()][cube.y()+1][cube.z()+1], dcg->dense3[cube.x()][cube.y()+1][cube.z() + 2] });
						nextCoface = Cube(birth, cube.x() - 1, cube.y(), cube.z(), 1);
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
			for(int i = position; i < 2; ++i){
				switch(i){
					case 0: // upper
						birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y()+1][cube.z() + 2], dcg->dense3[cube.x() + 2][cube.y()+1][cube.z() + 2],
							dcg->dense3[cube.x()+1][cube.y() + 2][cube.z() + 2],dcg->dense3[cube.x() + 2][cube.y() + 2][cube.z() + 2] });
						nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 0);
						break;

					case 1: // lower
						birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y()+1][cube.z()], dcg->dense3[cube.x() + 2][cube.y()+1][cube.z()],
							dcg->dense3[cube.x()+1][cube.y() + 2][cube.z()],dcg->dense3[cube.x() + 2][cube.y() + 2][cube.z()] });
						nextCoface = Cube(birth, cube.x(), cube.y(), cube.z() - 1, 0);
						break;
				}

				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;

			case 1: // dim2 type1 (fix y)
			for(int i = position; i < 2; ++i){
				switch(i){
					case 0: // left
						birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y() + 2][cube.z()+1], dcg->dense3[cube.x() + 2][cube.y() + 2][cube.z()+1],
							dcg->dense3[cube.x()+1][cube.y() + 2][cube.z() + 2],dcg->dense3[cube.x() + 2][cube.y() + 2][cube.z() + 2] });
						nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 0);
						break;

					case 1: //right
						birth = max({ cube.birth, dcg->dense3[cube.x()+1][cube.y()][cube.z()+1], dcg->dense3[cube.x() + 2][cube.y()][cube.z()+1],
							dcg->dense3[cube.x()+1][cube.y()][cube.z() + 2],dcg->dense3[cube.x() + 2][cube.y()][cube.z() + 2] });
						nextCoface = Cube(birth, cube.x(), cube.y() - 1, cube.z(), 0);
						break;
				}

				if (birth != dcg->threshold) {
					position = i + 1;
					return true;
				}
			}
			return false;

			case 2: // dim2 type2 (fix x)
			for(int i = position; i < 2; ++i){
				switch(i){
				case 0: // left
					birth = max({ cube.birth, dcg->dense3[cube.x() + 2][cube.y()+1][cube.z()+1], dcg->dense3[cube.x() + 2][cube.y() + 2][cube.z()+1],
						dcg->dense3[cube.x() + 2][cube.y()+1][cube.z() + 2],dcg->dense3[cube.x() + 2][cube.y() +2][cube.z() + 2] });
					nextCoface = Cube(birth, cube.x(), cube.y(), cube.z(), 0);
					break;

				case 1: //right
					birth = max({ cube.birth, dcg->dense3[cube.x()][cube.y()+1][cube.z()+1], dcg->dense3[cube.x()][cube.y() + 2][cube.z()+1],
						dcg->dense3[cube.x()][cube.y()+1][cube.z() + 2],dcg->dense3[cube.x()][cube.y() + 2][cube.z() + 2] });
					nextCoface = Cube(birth, cube.x() - 1, cube.y(), cube.z(), 0);
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
