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
#include <fstream>
#include <cassert>

#include "config.h"
#include "cube.h"
#include "npy.hpp"

using namespace std;

class DenseCubicalGrids{
public:
	Config *config;
	double threshold;
	uint8_t dim;
	uint32_t img_x, img_y, img_z, img_w;
	uint32_t ax, ay, az, aw;
    uint32_t axy, axyz, ayz, azw, axw, ayw;
	double*** dense3;
	double**** dense4;

    DenseCubicalGrids(Config&);
    // Overloaded constructor allowing explicit shape initialization
    DenseCubicalGrids(Config&, uint8_t dim, uint32_t ax, uint32_t ay = 1, uint32_t az = 1, uint32_t aw = 1);
	~DenseCubicalGrids(){
		if (dim < 4){
			free(dense3[0][0]);
			free(dense3[0]);
			free(dense3);
		} else {
			free(dense4[0][0][0]);
			free(dense4[0][0]);
			free(dense4[0]);
			free(dense4);
		}
	}
	double getBirth(uint32_t x, uint32_t y, uint32_t z);
	double getBirth(uint32_t x, uint32_t y, uint32_t z, uint32_t w, uint8_t cm, uint8_t dim);
	vector<uint32_t> ParentVoxel(uint8_t _dim, Cube &c);
	// allocate 3d array
	double ***alloc3d(uint32_t x, uint32_t y, uint32_t z) {
		double ***d = (double***)malloc(x * sizeof(double**));
		d[0] = (double**)malloc(x * y * sizeof(double*));
		d[0][0] = (double*)malloc(x*y*z * sizeof(double));
		for (uint32_t i = 0; i < x; i++) {
			d[i] = d[0] + i * y;
			for (uint32_t j = 0; j < y; j++) d[i][j] = d[0][0] + i * y*z + j * z;
		}
		if (d == NULL) {
			cerr << "not enough memory!" << endl;
		}
		return d;
	}
	// allocate 4d array
	double ****alloc4d(uint32_t x, uint32_t y, uint32_t z, uint32_t w) {
		#include <stdexcept>
		throw std::runtime_error("not implemented yet");
		double ****d = (double****)malloc(x * sizeof(double***));
		d[0] = (double***)malloc(x * y * sizeof(double**));
		d[0][0] = (double**)malloc(x*y*z * sizeof(double*));
		d[0][0][0] = (double*)malloc(x*y*z*w * sizeof(double));
		for (uint32_t i = 0; i < x; i++) {
			d[i] = d[0] + i * y;
			for (uint32_t j = 0; j < y; j++) {
				d[i][j] = d[0][0] + (i * y + j) * z;
				for (uint32_t k = 0; k < z; k++) d[i][j][k] = d[0][0][0] + (i * y*z + j * z + k) * w;
			}
		}
		if (d == NULL) {
			cerr << "not enough memory!" << endl;
		}
		return d;
	}
	// load image array from file
	void loadImage(bool embedded){
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
				assert(dim < 5);
				fin.read( ( char * ) &d, sizeof( int64_t ) );
				ax = d;
				if (dim>1) {
					fin.read( ( char * ) &d, sizeof( int64_t ) );
					ay = d;
				}else{
					ay = 1;
				}
				if (dim>2) {
					fin.read((char *)&d, sizeof(int64_t));
					az = d;
				}else {
					az = 1;
				}
				if (dim>3) {
					fin.read((char *)&d, sizeof(int64_t));
					aw = d;
				}else {
					aw = 1;
				}
				double dou;
				vector<double> arr;
				arr.reserve(ax*ay*az*aw);
				while (!fin.eof()){
					fin.read((char *)&dou, sizeof(double));
					arr.push_back(dou);
				}
				fin.close();
				gridFromArray(&arr[0], embedded, true);
				break;
			}

			case PERSEUS:
			{
				ifstream reading_file;
				reading_file.open(config->filename.c_str(), ios::in);
				string reading_line_buffer;
				getline(reading_file, reading_line_buffer);
				dim = atoi(reading_line_buffer.c_str());
				assert(dim < 5);
				getline(reading_file, reading_line_buffer);
				ax = atoi(reading_line_buffer.c_str());
				if (dim>1) {
					getline(reading_file, reading_line_buffer);
					ay = atoi(reading_line_buffer.c_str());
				}else {
					ay = 1;
				}
				if (dim>2) {
					getline(reading_file, reading_line_buffer);
					az = atoi(reading_line_buffer.c_str());
				}else {
					az = 1;
				}
				if (dim>3) {
					getline(reading_file, reading_line_buffer);
					aw = atoi(reading_line_buffer.c_str());
				}else {
					aw = 1;
				}
				vector<double> arr;
				arr.reserve(ax*ay*az*aw);
				while(!reading_file.eof()){
					getline(reading_file, reading_line_buffer);
					double dou = atof(reading_line_buffer.c_str());
					if (dou != -1) {
						arr.push_back(dou);
					}else{
						arr.push_back(config->threshold);
					}
				}
				reading_file.close();
				gridFromArray(&arr[0], embedded, true);
				break;
			}

			case CSV:
			{
				dim = 2;
				vector<double> arr;
				ifstream reading_file;
				reading_file.open(config->filename.c_str(), ios::in);
				string line;
				ay = 0;
				while (getline(reading_file, line)) {
					istringstream stream(line);
					string field;
					ax = 0;
					while (getline(stream, field, ',')) {
						arr.push_back(stod(field));
						ax++;
					}
					ay++;
				}
				az = 1;
				aw = 1;
				gridFromArray(&arr[0], embedded, true);
				break;
			}

			case NUMPY:
			{
				vector<unsigned long> shape;
				vector<double> arr;
				bool fortran_order;
				try{
					npy::LoadArrayFromNumpy(config->filename.c_str(), shape, fortran_order, arr);
				} catch (...) {
					cerr << "The data type of an numpy array should be numpy.float64." << endl;
					exit(-2);
				}
				if(shape.size() > 4){
					cerr << "Input array should be 1,2,3 or 4 dimensional " << endl;
					exit(-1);
				}
				dim = shape.size();
				ax = shape[0];
				if (dim>1) {
					ay = shape[1];
				}else {
					ay = 1;
				}
				if (dim>2) {
					az = shape[2];
				}else {
					az = 1;
				}
				if (dim>3) {
					aw = shape[3];
				}else {
					aw = 1;
				}
				gridFromArray(&arr[0], embedded, fortran_order);
				break;
			}
		}
		if (dim < 4)
			cout << "x : y : z = " << img_x << " : " << img_y << " : " << img_z << endl;
		else
			cout << "x : y : z : w = " << img_x << " : " << img_y << " : " << img_z << " : " << img_w << endl;
		// T-construction (the number of vertices = that of the top cells in each dimension)
		if(config->tconstruction){
			if(dim>3) aw++;
			if(dim>2) az++;
			ax++;
			ay++;
		}
		//
		axy = ax * ay;
		ayz = ay * az;
		azw = az * aw;
		axw = ax * aw;
		ayw = ay * aw;
		axyz = ax * ay * az;
	}

	// construct volume with boundary
	void gridFromArray(const double *arr, bool embedded, bool fortran_order){
		img_x = ax;
		img_y = ay;
		img_z = az;
		img_w = aw;
		uint64_t i = 0;
		uint32_t x_shift = 2;
		uint32_t y_shift = 2;
		uint32_t z_shift = 2;
		uint32_t w_shift = 2;
		double sgn = 1;
		if(embedded){
			sgn = -1;
			x_shift = 4;
			y_shift = 4;
			if (az>1) z_shift = 4;
			if (aw>1) w_shift = 4;
		}
		if (dim < 4){
			dense3 = alloc3d(ax + x_shift, ay + y_shift, az + z_shift);
			if(fortran_order){
				for (uint32_t z = 0; z < az + z_shift; ++z) {
					for (uint32_t y = 0; y < ay + y_shift; ++y) {
						for (uint32_t x = 0; x < ax + x_shift; ++x) {
							if (x_shift/2-1 < x && x <= ax+x_shift/2-1 && y_shift/2-1 < y && y <= ay+y_shift/2-1 && z_shift/2-1 < z && z<= az + z_shift/2-1) {
								dense3[x][y][z] = sgn * arr[i++];
							}else if (0 == x || x == ax-1+y_shift || 0 == y || y == ay-1+y_shift || z==0 || z==az-1+z_shift) { // outer boundary
								dense3[x][y][z] = config->threshold;
							}else{  // only for embedded; inner boundary
								dense3[x][y][z] = -config->threshold;
							}
						}
					}
				}
			}else{
				for (uint32_t x = 0; x < ax + x_shift; ++x) {
					for (uint32_t y = 0; y < ay + y_shift; ++y) {
						for (uint32_t z = 0; z < az + z_shift; ++z) {
							if (x_shift/2-1 < x && x <= ax+x_shift/2-1 && y_shift/2-1 < y && y <= ay+y_shift/2-1 && z_shift/2-1 < z && z<= az + z_shift/2-1) {
								dense3[x][y][z] = sgn * arr[i++];
							}else if (0 == x || x == ax-1+y_shift || 0 == y || y == ay-1+y_shift || z==0 || z==az-1+z_shift) { // outer boundary
								dense3[x][y][z] = config->threshold;
							}else{  // only for embedded; inner boundary
								dense3[x][y][z] = -config->threshold;
							}
						}
					}
				}
			}
		} else { // dim==4
			dense4 = alloc4d(ax + x_shift, ay + y_shift, az + z_shift, aw + w_shift);
			if(fortran_order){
				for (uint32_t w = 0; w < aw + w_shift; ++w) {
					for (uint32_t z = 0; z < az + z_shift; ++z) {
						for (uint32_t y = 0; y < ay + y_shift; ++y) {
							for (uint32_t x = 0; x < ax + x_shift; ++x) {
								if (x_shift/2-1 < x && x <= ax+x_shift/2-1 && y_shift/2-1 < y && y <= ay+y_shift/2-1 && z_shift/2-1 < z && z<= az+z_shift/2-1 && w_shift/2-1 < w && w<= aw+w_shift/2-1) {
									dense4[x][y][z][w] = sgn * arr[i++];
								}else if (0==x || x==ax-1+y_shift || 0==y || y==ay-1+y_shift || 0==z || z==az-1+z_shift || 0==w || w==aw-1+w_shift) { // outer boundary
									dense4[x][y][z][w] = config->threshold;
								}else{  // only for embedded; inner boundary
									dense4[x][y][z][w] = -config->threshold;
								}
							}
						}
					}
				}
			}else{
				for (uint32_t x = 0; x < ax + x_shift; ++x) {
					for (uint32_t y = 0; y < ay + y_shift; ++y) {
						for (uint32_t z = 0; z < az + z_shift; ++z) {
							for (uint32_t w = 0; w < aw + w_shift; ++w) {
								if (x_shift/2-1 < x && x <= ax+x_shift/2-1 && y_shift/2-1 < y && y <= ay+y_shift/2-1 && z_shift/2-1 < z && z<= az+z_shift/2-1 && w_shift/2-1 < w && w<= aw+w_shift/2-1) {
									dense4[x][y][z][w] = sgn * arr[i++];
								}else if (0==x || x==ax-1+y_shift || 0==y || y==ay-1+y_shift || 0==z || z==az-1+z_shift || 0==w || w==aw-1+w_shift) { // outer boundary
									dense4[x][y][z][w] = config->threshold;
								}else{  // only for embedded; inner boundary
									dense4[x][y][z][w] = -config->threshold;
								}
							}
						}
					}
				}
			}
		}
		ax += x_shift-2;
		ay += y_shift-2;
		az += z_shift-2;
		aw += w_shift-2;
	}
};
