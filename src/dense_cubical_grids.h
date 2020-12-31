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
#include <iostream>

#include "config.h"
#include "cube.h"

using namespace std;

class DenseCubicalGrids{
public:
	Config *config;
	double threshold;
	uint8_t dim;
	uint32_t img_x, img_y, img_z;
	uint32_t ax, ay, az;
    uint32_t axy, axyz, ayz;
	double*** dense3;

	DenseCubicalGrids(Config&);
	~DenseCubicalGrids();
	void loadImage(bool embedded);
	void gridFromArray(vector<double>& arr, bool embedded);	
	void gridFromNpyArray(const double *arr, bool embedded);	
	double getBirth(uint32_t x, uint32_t y, uint32_t z);
	double getBirth(uint32_t x, uint32_t y, uint32_t z, uint8_t cm, uint8_t dim);
	vector<uint32_t> ParentVoxel(uint8_t _dim, Cube &c);
	// allocate 3d array
	double ***alloc3d(uint32_t x, uint32_t y, uint32_t z) {
		double ***d = (double***)malloc(x * sizeof(double**));
		d[0] = (double**)malloc(x * y * sizeof(double*));
		d[0][0] = (double*)malloc(x*y*z * sizeof(double));
		for (uint32_t i = 0; i < x ; i++) {
			d[i] = d[0] + i * y;
			for (uint32_t j = 0; j < y; j++) d[i][j] = d[0][0] + i * y*z + j * z;
		}
		if (d == NULL) {
			cerr << "not enough memory!" << endl;
		}
		return d;
	}
};


