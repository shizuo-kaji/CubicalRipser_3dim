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
#include <initializer_list>

#include "dense_cubical_grids.h"
#include "npy.hpp"

using namespace std;


DenseCubicalGrids::DenseCubicalGrids(Config& _config)  {
	config = &_config;
	threshold = config->threshold;
}


// read from file
void DenseCubicalGrids::loadImage(bool embedded){
	// read file
	cout << "Reading " << config->filename << endl;
	switch(config->format){
		case DIPHA:
		{
			ifstream fin( config->filename, ios::in | ios::binary );

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
			dense3 = alloc3d(ax+2, ay+2, az+2);
			cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;

			double dou;
			if(embedded){ // dual complex
				if(az>1){
					dense3 = alloc3d(ax+4, ay+4, az+4);
					for (uint32_t z = 0; z < az + 4; ++z) {
						for (uint32_t y = 0; y < ay + 4; ++y) {
							for (uint32_t x = 0; x < ax + 4; ++x) {
								if (1 < x && x <= ax+1 && 1 < y && y <= ay+1 && 1 < z && z <= az+1) {
									if (!fin.eof()) {
										fin.read((char *)&dou, sizeof(double));
										dense3[x][y][z] = -dou;
									}
									else {
										cerr << "file endof error " << endl;
									}
								}else if (0 == x || x == ax+3 || 0 == y || y == ay+3 || 0 == z || z == az+3) {
									dense3[x][y][z] = config->threshold;
								}else{
									dense3[x][y][z] = -config->threshold;
								}
			//					cout << x << "," << y << "," << z << ": " << dense3[x][y][z] << endl;
							}
						}
					}
					az = az + 2;
				}else{
					dense3 = alloc3d(ax+4, ay+4, az+2);
					for (uint32_t z = 0; z < az + 2; ++z) { // z-axis remains the same
						for (uint32_t y = 0; y < ay + 4; ++y) {
							for (uint32_t x = 0; x < ax + 4; ++x) {
								if (1 < x && x <= ax+1 && 1 < y && y <= ay+1 && 1<=z && z<= az) {
									if (!fin.eof()) {
										fin.read((char *)&dou, sizeof(double));
										dense3[x][y][z] = -dou;
									}
									else {
										cerr << "file endof error " << endl;
									}
								}else if (0 == x || x == ax+3 || 0 == y || y == ay+3 || z==0 || z==az+1) {
									dense3[x][y][z] = config->threshold;
								}else{
									dense3[x][y][z] = -config->threshold;
								}
								cout << x << "," << y << "," << z << ": " << dense3[x][y][z] << endl;
							}
						}
					}
				}
				ax = ax + 2;
				ay = ay + 2;				
			}else{
				for (uint32_t z = 0; z < az + 2; ++z) {
					for (uint32_t y = 0; y < ay + 2; ++y) {
						for (uint32_t x = 0; x < ax + 2; ++x) {
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
								dense3[x][y][z] = config->threshold;
							}
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
			reading_file.open(config->filename.c_str(), ios::in); 

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
			dense3 = alloc3d(ax+2, ay+2, az+2);
			cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;
			if(embedded){ // dual complex
				if(az>1){
					dense3 = alloc3d(ax+4, ay+4, az+4);
					for (uint32_t z = 0; z < az + 4; ++z) {
						for (uint32_t y = 0; y < ay + 4; ++y) {
							for (uint32_t x = 0; x < ax + 4; ++x) {
								if (1 < x && x <= ax+1 && 1 < y && y <= ay+1 && 1 < z && z <= az+1) {
									if (!reading_file.eof()) {
										getline(reading_file, reading_line_buffer);
										double dou = atof(reading_line_buffer.c_str());
										if (dou == -1) {
											dense3[x][y][z] = config->threshold;
										}
										else {
											dense3[x][y][z] = -dou;
										}
									}
								}else if (0 == x || x == ax+3 || 0 == y || y == ay+3 || 0 == z || z == az+3) {
									dense3[x][y][z] = config->threshold;
								}else{
									dense3[x][y][z] = -config->threshold;
								}
			//					cout << x << "," << y << "," << z << ": " << dense3[x][y][z] << endl;
							}
						}
					}
					az = az + 2;
				}else{
					dense3 = alloc3d(ax+4, ay+4, az+2);
					for (uint32_t z = 0; z < az + 2; ++z) { // z-axis remains the same
						for (uint32_t y = 0; y < ay + 4; ++y) {
							for (uint32_t x = 0; x < ax + 4; ++x) {
								if (1 < x && x <= ax+1 && 1 < y && y <= ay+1 && 1<=z && z<= az) {
									if (!reading_file.eof()) {
										getline(reading_file, reading_line_buffer);
										double dou = atof(reading_line_buffer.c_str());
										if (dou == -1) {
											dense3[x][y][z] = config->threshold;
										}
										else {
											dense3[x][y][z] = -dou;
										}
									}
								}else if (0 == x || x == ax+3 || 0 == y || y == ay+3 || z==0 || z==az+1) {
									dense3[x][y][z] = config->threshold;
								}else{
									dense3[x][y][z] = -config->threshold;
								}
								cout << x << "," << y << "," << z << ": " << dense3[x][y][z] << endl;
							}
						}
					}
				}
				ax = ax + 2;
				ay = ay + 2;				
			}else{
				for (uint32_t z = 0; z < az + 2; ++z) {
					for (uint32_t y = 0; y < ay + 2; ++y) {
						for (uint32_t x = 0; x < ax + 2; ++x) {
							if (0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az) {
								if (!reading_file.eof()) {
									getline(reading_file, reading_line_buffer);
									double dou = atof(reading_line_buffer.c_str());
									if (dou == -1) {
										dense3[x][y][z] = config->threshold;
									}
									else {
										dense3[x][y][z] = dou;
									}
								}
							}
							else {
								dense3[x][y][z] = config->threshold;
							}
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
			npy::LoadArrayFromNumpy(config->filename.c_str(), shape, data);
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
			uint64_t i = 0;
			// note the order of axis
			if(embedded){ // dual complex
				if(az>1){
					dense3 = alloc3d(ax+4, ay+4, az+4);
					for (uint32_t x = 0; x < ax + 4; ++x) {
						for (uint32_t y = 0; y < ay + 4; ++y) {
							for (uint32_t z = 0; z < az + 4; ++z) {
								if (1 < x && x <= ax+1 && 1 < y && y <= ay+1 && 1 < z && z <= az+1) {
									dense3[x][y][z] = -data[i++];
								}else if (0 == x || x == ax+3 || 0 == y || y == ay+3 || 0 == z || z == az+3) {
									dense3[x][y][z] = config->threshold;
								}else{
									dense3[x][y][z] = -config->threshold;
								}
			//					cout << x << "," << y << "," << z << ": " << dense3[x][y][z] << endl;
							}
						}
					}
					az = az + 2;
				}else{
					dense3 = alloc3d(ax+4, ay+4, az+2);
					for (uint32_t x = 0; x < ax + 4; ++x) {
						for (uint32_t y = 0; y < ay + 4; ++y) {
							for (uint32_t z = 0; z < az + 2; ++z) { // z-axis remains the same
								if (1 < x && x <= ax+1 && 1 < y && y <= ay+1 && 1<=z && z<= az) {
									dense3[x][y][z] = -data[i++];
								}else if (0 == x || x == ax+3 || 0 == y || y == ay+3 || z==0 || z==az+1) {
									dense3[x][y][z] = config->threshold;
								}else{
									dense3[x][y][z] = -config->threshold;
								}
//								cout << x << "," << y << "," << z << ": " << dense3[x][y][z] << endl;
							}
						}
					}
				}
				ax = ax + 2;
				ay = ay + 2;				
			}else{
				dense3 = alloc3d(ax + 2, ay + 2, az + 2);
				for (uint32_t x = 0; x < ax + 2; ++x) {
					for (uint32_t y = 0; y <ay + 2; ++y) {
						for (uint32_t z = 0; z < az + 2; ++z) {
							if (0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az) {
								dense3[x][y][z] = data[i++];
							}
							else { // fill the boundary with the threashold value
								dense3[x][y][z] = config->threshold;
							}
						}
					}
				}
				break;
			}
		}
	}	
	axy = ax * ay;
	ayz = ay * az;
	axyz = ax * ay * az;
//	cout << ax << "," << ay << "," << az << endl;
}

// return filtlation value for a cube
double DenseCubicalGrids::getBirth(uint32_t cx, uint32_t cy, uint32_t cz, uint8_t cm, uint8_t dim) {
	// beware of the shift due to the boundary
	switch (dim) {
		case 0:
			return dense3[cx+1][cy+1][cz+1];
		case 1:
			switch (cm) {
			case 0:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 1][cz + 1]);
			case 1:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 1][cy + 2][cz + 1]);
			case 2:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 1][cy + 1][cz + 2]);
			case 3:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 2][cz + 1]);
			case 4:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 0][cz + 1]);
			// for 3d dual only
			case 5:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 1][cy + 0][cz + 2]);
			case 6:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 1][cy + 2][cz + 2]);
			case 7:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 0][cz + 2]);
			case 8:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 1][cz + 2]);
			case 9:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 2][cz + 2]);
			case 10:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 0][cz + 0]);
			case 11:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 1][cz + 0]);
			case 12:
				return max(dense3[cx+1][cy+1][cz+1], dense3[cx + 2][cy + 2][cz + 0]);
			}
		case 2:
			switch (cm) {
			case 0: // x - y (fix z)
				return max({ dense3[cx+1][cy+1][cz+1], dense3[cx+2][cy+1][cz+1],
					dense3[cx+2][cy+2][cz+1], dense3[cx+1][cy+2][cz+1] });
			case 1: // z - x (fix y)
				return max({ dense3[cx+1][cy+1][cz+1], dense3[cx+1][cy+1][cz+2],
					dense3[cx+2][cy+1][cz+2], dense3[cx+2][cy+1][cz+1] });
			case 2: // y - z (fix x)
				return max({ dense3[cx+1][cy+1][cz+1], dense3[cx+1][cy+2][cz+1],
					dense3[cx+1][cy+2][cz+2], dense3[cx+1][cy+1][cz+2] });
			}
		case 3:
			return max({ dense3[cx+1][cy+1][cz+1], dense3[cx+2][cy+1][cz+1],
				dense3[cx+2][cy+2][cz+1], dense3[cx+1][cy+2][cz+1],
				dense3[cx+1][cy+1][cz+2], dense3[cx+2][cy+1][cz+2],
				dense3[cx+2][cy+2][cz+2], dense3[cx+1][cy+2][cz+2] });
		}
	return threshold; // dim > 3
}

// allocate 3d array
double ***DenseCubicalGrids::alloc3d(uint32_t x, uint32_t y, uint32_t z) {
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

DenseCubicalGrids::~DenseCubicalGrids(){
	free(dense3[0][0]);
	free(dense3[0]);
	free(dense3);
}
