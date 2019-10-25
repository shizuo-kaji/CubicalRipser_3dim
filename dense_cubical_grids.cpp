/* dense_cubical_grids.cpp

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>

#include "dense_cubical_grids.h"
#include "npy.hpp"

using namespace std;

double ***alloc3d(int x, int y, int z) {
	double ***dense3 = (double***)malloc(x * sizeof(double**));
	dense3[0] = (double**)malloc(y * z * sizeof(double*));
	dense3[0][0] = (double*)malloc(x*y*z * sizeof(double));
	for (int i = 0; i < x ; i++) {
		dense3[i] = dense3[0] + i * x;
		for (int j = 0; j < y; j++) dense3[i][j] = dense3[0][0] + i * y*z + j * z;
	}
	if (dense3 == NULL) {
		cerr << "not enough memory!" << endl;
	}
	return dense3;
}
DenseCubicalGrids::DenseCubicalGrids(const string& filename, double _threshold, file_format format)  {

	threshold = _threshold;

	// read file
	cout << filename << endl;
	switch(format){
		case DIPHA:
		{
			ifstream fin( filename, ios::in | ios::binary );

			int64_t d;
			fin.read( ( char * ) &d, sizeof( int64_t ) ); // magic number
			assert(d == 8067171840);
			fin.read( ( char * ) &d, sizeof( int64_t ) ); // type number
			assert(d == 1);
			fin.read( ( char * ) &d, sizeof( int64_t ) ); //data num
			fin.read( ( char * ) &d, sizeof( int64_t ) ); // dim 
			dim = d;
			assert(dim < 4);
			fin.read( ( char * ) &d, sizeof( int64_t ) );
			ax = d;
			fin.read( ( char * ) &d, sizeof( int64_t ) );
			ay = d;
			if (dim == 3) {
				fin.read((char *)&d, sizeof(int64_t));
				az = d;
			}else {
				az = 1;
			}
			dense3 = alloc3d(ax, ay, az);
			cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;

			double dou;
			for (int z = 0; z < az + 2; ++z) {
				for (int y = 0; y < ay + 2; ++y) {
					for (int x = 0; x < ax + 2; ++x) {
						if (0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az) {
							if (!fin.eof()) {
								fin.read((char *)&dou, sizeof(double));
								dense3[x][y][z] = dou;
							}
							else {
								cerr << "file endof error " << endl;
							}
						}
						else {
							dense3[x][y][z] = threshold;
						}
					}
				}
			}
			fin.close();
			break;
		}

		case PERSEUS:
		{
			ifstream reading_file; 
			reading_file.open(filename.c_str(), ios::in); 

			string reading_line_buffer; 
			getline(reading_file, reading_line_buffer); 
			dim = atoi(reading_line_buffer.c_str());
			assert(dim < 4);
			getline(reading_file, reading_line_buffer);
			ax = atoi(reading_line_buffer.c_str()); 
			getline(reading_file, reading_line_buffer); 
			ay = atoi(reading_line_buffer.c_str()); 
			if (dim == 3) {
				getline(reading_file, reading_line_buffer);
				az = atoi(reading_line_buffer.c_str());
			}else {
				az = 1;
			}
			dense3 = alloc3d(ax, ay, az);
			cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;

			for (int z = 0; z < az + 2; ++z) {
				for (int y = 0; y < ay + 2; ++y) {
					for (int x = 0; x < ax + 2; ++x) {
						if (0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az) {
							if (!reading_file.eof()) {
								getline(reading_file, reading_line_buffer);
								double dou = atof(reading_line_buffer.c_str());
								if (dou == -1) {
									dense3[x][y][z] = threshold;
								}
								else {
									dense3[x][y][z] = dou;
								}
							}
						}
						else {
							dense3[x][y][z] = threshold;
						}
					} 
				}
			}
			reading_file.close();
			break;
		}
		case NUMPY:
		{
			vector<unsigned long> shape;
			vector<double> data;
			npy::LoadArrayFromNumpy(filename.c_str(), shape, data);
			if(shape.size() > 3){
				cerr << "Input array should be 2 or 3 dimensional " << endl;
				exit(-1);
			}
			dim = shape.size();
			ax = shape[0];
			ay = shape[1];
			if (dim == 3) {
				az = shape[2];
			}
			else {
				az = 1;
			}
			cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;
			int i = 0;
			for (int x = 0; x < ax + 2; ++x) {
				for (int y = 0; y <ay + 2; ++y) {
					for (int z = 0; z < az + 2; ++z) {
						if (0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az) {
							dense3[x][y][z] = data[i++];
						}
						else { // fill the boundary with the threashold value
							dense3[x][y][z] = threshold;
						}
					}
				}
			}
			break;
		}
	}
	
	axy = ax * ay;
	ayz = ay * az;
	axyz = ax * ay * az;
}

double DenseCubicalGrids::getBirthday(int cx, int cy, int cz, int cm, int dim) {
	switch (dim) {
		case 0:
			return get(cx,cy,cz);
		case 1:
			switch (cm) {
			case 0:
				return max(get(cx,cy,cz), get(cx+1, cy, cz));
			case 1:
				return max(get(cx,cy,cz), get(cx, cy+1, cz));
			case 2:
				return max(get(cx,cy,cz), get(cx, cy, cz+1));
			}
		case 2:
			switch (cm) {
			case 0: // x - y (fix z)
				return max({ get(cx,cy,cz), get(cx+1,cy,cz),
					get(cx+1,cy+1,cz), get(cx,cy+1,cz) });
			case 1: // z - x (fix y)
				return max({ get(cx,cy,cz), get(cx,cy,cz+1),
					get(cx+1,cy,cz+1), get(cx+1,cy,cz) });
			case 2: // y - z (fix x)
				return max({ get(cx,cy,cz), get(cx,cy+1,cz),
					get(cx,cy+1,cz+1), get(cx,cy,cz+1) });
			}
		case 3:
			return max({ get(cx,cy,cz), get(cx+1,cy,cz),
				get(cx+1,cy+1,cz), get(cx,cy+1,cz),
				get(cx,cy,cz+1), get(cx+1,cy,cz+1),
				get(cx+1,cy+1,cz+1), get(cx,cy+1,cz+1) });
		}
	return threshold; // dim > 3
}

// conversion from id to coordinates and type
vector<int> DenseCubicalGrids::getXYZM(long index) {
	vector<int> loc(4);   // (x,y,z,m)
	loc[0] = index % ax;
	loc[1] = (index / ax) % ay;
	loc[2] = (index / axy) % az;
	loc[3] = (index / axyz);
	return(loc);
}
