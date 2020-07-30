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
#include "config.h"

using namespace std;

class DenseCubicalGrids{
public:
	Config *config;
	double threshold;
	uint8_t dim;
	uint32_t ax, ay, az;
    uint32_t axy, axyz, ayz;
	double*** dense3;

	DenseCubicalGrids(Config&);
	~DenseCubicalGrids();
	void loadImage(bool embedded);
	void gridFromArray(vector<double>& arr, bool embedded);	
	void gridFromNpyArray(const double *arr, bool embedded);	
	double ***alloc3d(uint32_t x, uint32_t y, uint32_t z);
	double getBirth(uint32_t x, uint32_t y, uint32_t z, uint8_t cm, uint8_t dim);
};
