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
private:
	vector<double> dense3;

public:
	double threshold;
	int dim;
	int ax, ay, az;
	file_format format;

	DenseCubicalGrids(const std::string& filename, double _threshold, file_format _format);
	~DenseCubicalGrids();

	double getBirthday(int x, int y, int z, int cm, int dim);
	vector<int> getXYZM(long index);
	long getIndex(int x, int y, int z, int cm=0);
	void set(int x, int y, int z, double val) {
		dense3[x*ay*az + y * az + z] = val;
	}
	double get(int x, int y, int z) {
		//assert(-1 <= x && x <= ax && -1 <= y && y <= ay && -1 <= z && z <= az);
		if (x == -1 || y == -1 || z == -1 || x == ax || y == ay || z == az) {
			return(threshold);
		}
		else {
			return(dense3[x*ay*az + y * az + z]);
		}
	}

};