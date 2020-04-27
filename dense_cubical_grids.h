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
	double threshold;
	unsigned short dim;
	unsigned short ax, ay, az;
    int axy, axyz, ayz;
	double*** dense3;

	DenseCubicalGrids(const std::string& filename, double _threshold, file_format format);
	~DenseCubicalGrids();
	
	double getBirth(unsigned short x, unsigned short y, unsigned short z, unsigned short cm, unsigned short dim);
};
