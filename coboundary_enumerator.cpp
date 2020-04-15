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
	nextCoface = Cube(0, -1);
	dcg = _dcg;
    dim = _dim;
}
	
void CoboundaryEnumerator::setCoboundaryEnumerator(Cube& _s) {
	cube = _s;
	position = 0; // current position of coface search
}

bool CoboundaryEnumerator::hasNextCoface() {
	double birthday=0;
    vector<int> loc = dcg -> getXYZM(cube.index);
    int cx = loc[0];
	int cy = loc[1];
	int cz = loc[2];
	int cm = loc[3];
	long index=0;
	// note the shift for the boundary
	switch (dim) {
		case 0: // dim0
		for (int i = position; i < 6; ++i) {
			switch (i){
			case 0:
				index= dcg->getIndex(cx, cy, cz, 2);
				birthday = max(cube.birthday, dcg->dense3[cx+1][cy+1][cz+2]);
				break;

			case 1:
				index= dcg->getIndex(cx, cy, cz - 1, 2);
				birthday = max(cube.birthday, dcg->dense3[cx+1][cy+1][cz]);
				break;

			case 2:
				index= dcg->getIndex(cx, cy, cz, 1);
				birthday = max(cube.birthday, dcg->dense3[cx+1][cy + 2][cz+1]);
				break;

			case 3:
				index= dcg->getIndex(cx, cy - 1, cz, 1);
				birthday = max(cube.birthday, dcg->dense3[cx+1][cy][cz+1]);
				break;

			case 4:
				index= dcg->getIndex(cx, cy, cz, 0);
				birthday = max(cube.birthday, dcg->dense3[cx + 2][cy+1][cz+1]);
				break;

			case 5:
				index= dcg->getIndex(cx - 1, cy, cz, 0);
				birthday = max(cube.birthday, dcg->dense3[cx][cy+1][cz+1]);
				break;
			}

			if (birthday != dcg->threshold) {
				position = i + 1;
				nextCoface = Cube(birthday, index);
				return true;
			}
		}
		return false;

		case 1: // dim1
		switch (cm) {
			case 0: // dim1 type0 (x-axis -> )
			for(int i = position; i < 4; ++i){
				switch(i){
				case 0:
					index= dcg->getIndex(cx, cy, cz, 1);
					birthday = max({ cube.birthday, dcg->dense3[cx+1][cy+1][cz + 2], dcg->dense3[cx + 2][cy+1][cz + 2] });
					break;

				case 1:
					index= dcg->getIndex(cx, cy, cz - 1, 1);
					birthday = max({ cube.birthday, dcg->dense3[cx+1][cy+1][cz], dcg->dense3[cx + 2][cy+1][cz] });
					break;

				case 2:
					index= dcg->getIndex(cx, cy, cz, 0);
					birthday = max({ cube.birthday, dcg->dense3[cx+1][cy + 2][cz+1], dcg->dense3[cx + 2][cy + 2][cz+1] });
					break;

				case 3:
					index= dcg->getIndex(cx, cy - 1, cz, 0);
					birthday = max({ cube.birthday, dcg->dense3[cx+1][cy][cz+1], dcg->dense3[cx + 2][cy][cz+1] });
					break;
				}

				if (birthday != dcg->threshold) {
					position = i + 1;
					nextCoface = Cube(birthday, index);
					return true;
				}
			}
			return false;

			case 1: // dim1 type1 (y-axis -> )
			for(int i = position; i < 4; ++i){
				switch(i){
				case 0:
					index= dcg->getIndex(cx, cy, cz, 2);
					birthday = max({ cube.birthday, dcg->dense3[cx+1][cy+1][cz + 2], dcg->dense3[cx+1][cy + 2][cz + 2] });
					break;

				case 1:
					index= dcg->getIndex(cx, cy, cz - 1, 2);
					birthday = max({ cube.birthday, dcg->dense3[cx+1][cy+1][cz], dcg->dense3[cx+1][cy + 2][cz] });
					break;

				case 2:
					index= dcg->getIndex(cx, cy, cz, 0);
					birthday = max({ cube.birthday, dcg->dense3[cx + 2][cy+1][cz+1], dcg->dense3[cx + 2][cy + 2][cz+1] });
					break;

				case 3:
					index= dcg->getIndex(cx - 1, cy, cz, 0);
					birthday = max({ cube.birthday, dcg->dense3[cx][cy+1][cz+1], dcg->dense3[cx][cy + 2][cz+1] });
					break;
				}

				if (birthday != dcg->threshold) {
					position = i + 1;
					nextCoface = Cube(birthday, index);
					return true;
				}
			}
			return false;

			case 2: // dim1 type2 (z-axis -> )
			for(int i = position; i < 4; ++i){
				switch(i){
					case 0:
						index= dcg->getIndex(cx, cy, cz, 2);
						birthday = max({ cube.birthday, dcg->dense3[cx+1][cy + 2][cz+1], dcg->dense3[cx+1][cy + 2][cz + 2] });
						break;

					case 1:
						index= dcg->getIndex(cx, cy - 1, cz, 2);
						birthday = max({ cube.birthday, dcg->dense3[cx+1][cy][cz+1], dcg->dense3[cx+1][cy][cz + 2] });
						break;

					case 2:
						index= dcg->getIndex(cx, cy, cz, 1);
						birthday = max({ cube.birthday, dcg->dense3[cx + 2][cy+1][cz+1], dcg->dense3[cx + 2][cy+1][cz + 2] });
						break;

					case 3:
						index= dcg->getIndex(cx - 1, cy, cz, 1);
						birthday = max({ cube.birthday, dcg->dense3[cx][cy+1][cz+1], dcg->dense3[cx][cy+1][cz + 2] });
						break;
				}

				if (birthday != dcg->threshold) {
					position = i + 1;
					nextCoface = Cube(birthday, index);
					return true;
				}
			}
			return false;
		}
		return false;

		case 2: // dim2
		switch (cm) {
			case 0: // dim2 type0 (fix z)
			for(int i = position; i < 2; ++i){
				switch(i){
					case 0: // upper
						index= dcg->getIndex(cx, cy, cz, 0);
						birthday = max({ cube.birthday, dcg->dense3[cx+1][cy+1][cz + 2], dcg->dense3[cx + 2][cy+1][cz + 2],
							dcg->dense3[cx+1][cy + 2][cz + 2],dcg->dense3[cx + 2][cy + 2][cz + 2] });
						break;

					case 1: // lower
						index= dcg->getIndex(cx, cy, cz - 1, 0);
						birthday = max({ cube.birthday, dcg->dense3[cx+1][cy+1][cz], dcg->dense3[cx + 2][cy+1][cz],
							dcg->dense3[cx+1][cy + 2][cz],dcg->dense3[cx + 2][cy + 2][cz] });
						break;
				}

				if (birthday != dcg->threshold) {
					position = i + 1;
					nextCoface = Cube(birthday, index);
					return true;
				}
			}
			return false;

			case 1: // dim2 type1 (fix y)
			for(int i = position; i < 2; ++i){
				switch(i){
					case 0: // left
						index= dcg->getIndex(cx, cy, cz, 0);
						birthday = max({ cube.birthday, dcg->dense3[cx+1][cy + 2][cz+1], dcg->dense3[cx + 2][cy + 2][cz+1],
							dcg->dense3[cx+1][cy + 2][cz + 2],dcg->dense3[cx + 2][cy + 2][cz + 2] });
						break;

					case 1: //right
						index= dcg->getIndex(cx, cy - 1, cz, 0);
						birthday = max({ cube.birthday, dcg->dense3[cx+1][cy][cz+1], dcg->dense3[cx + 2][cy][cz+1],
							dcg->dense3[cx+1][cy][cz + 2],dcg->dense3[cx + 2][cy][cz + 2] });
						break;
				}

				if (birthday != dcg->threshold) {
					position = i + 1;
					nextCoface = Cube(birthday, index);
					return true;
				}
			}
			return false;

			case 2: // dim2 type2 (fix x)
			for(int i = position; i < 2; ++i){
				switch(i){
				case 0: // left
					index= dcg->getIndex(cx, cy, cz, 0);
					birthday = max({ cube.birthday, dcg->dense3[cx + 2][cy+1][cz+1], dcg->dense3[cx + 2][cy + 2][cz+1],
						dcg->dense3[cx + 2][cy+1][cz + 2],dcg->dense3[cx + 2][cy +2][cz + 2] });
					break;

				case 1: //right
					index= dcg->getIndex(cx - 1, cy, cz, 0);
					birthday = max({ cube.birthday, dcg->dense3[cx][cy+1][cz+1], dcg->dense3[cx][cy + 2][cz+1],
						dcg->dense3[cx][cy+1][cz + 2],dcg->dense3[cx][cy + 2][cz + 2] });
					break;
				}

				if (birthday != dcg->threshold) {
					position = i + 1;
					nextCoface = Cube(birthday, index);
					return true;
				}
			}
			return false;
		}
	}
	return false;
}
