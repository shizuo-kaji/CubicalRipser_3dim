/* dense_cubical_grids.h

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <string>
#include <vector>

using namespace std;

enum file_format { DIPHA, PERSEUS, NUMPY };

class DenseCubicalGrids {

public:
	double*** dense3;
	double threshold;
	int dim;
	int ax, ay, az;
	long axy, axyz, ayz;

	DenseCubicalGrids(const std::string& filename, double _threshold, file_format format);

	double getBirthday(int x, int y, int z, int cm, int dim);
	vector<int> getXYZM(long index);

	// unique id for each simplex (unique only within a single dimension)
	long getIndex(int x, int y, int z, int cm=0){
		return(x + y * ax + z * axy + cm * axyz);
	}

	double get(int x, int y, int z) {
			return(dense3[x+1][y+1][z+1]);
	}

};