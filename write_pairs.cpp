/* write_pairs.cpp

Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.

This file is part of CubicalRipser_3dim.

CubicalRipser: C++ system for computation of Cubical persistence pairs
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
CubicalRipser is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

CubicalRipser is deeply depending on 'Ripser', software for Vietoris-Rips 
persitence pairs by Ulrich Bauer, 2015-2016.  We appreciate Ulrich very much.
We rearrange his codes of Ripser and add some new ideas for optimization on it 
and modify it for calculation of a Cubical filtration.

This part of CubicalRiper is a calculator of cubical persistence pairs for 
3 dimensional pixel data. The input data format conforms to that of DIPHA.
 See more descriptions in README.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "write_pairs.h"
#include <iostream>
#include "array_index.h"

WritePairs::WritePairs(int64_t _dim, double _birth, double _death, int64_t _birth_x,int64_t _birth_y,int64_t _birth_z, bool print){
	dim = _dim;
	birth = _birth;
	death = _death;
	birth_x = _birth_x;
	birth_y = _birth_y;
	birth_z = _birth_z;
	if (print == true) {
		std::cout << "[" << birth << "," << death << ")" << " birth loc (" << birth_x << "," << birth_y << "," << birth_z << ")" << std::endl;
	}
}

WritePairs::WritePairs(int64_t _dim, double _birth, double _death, long idx, bool print) {
	dim = _dim;
	birth = _birth;
	death = _death;
	ArrayIndex index(idx);
	birth_x = index.cx;
	birth_y = index.cy;
	birth_z = index.cz;
	if (print == true) {
		std::cout << "[" << birth << "," << death << ")" << " birth loc (" << birth_x << "," << birth_y << "," << birth_z << ")" << std::endl;
	}
}


int64_t WritePairs::getDimension(){
	return dim;
}

double WritePairs::getBirth(){
	return birth;
}

double WritePairs::getDeath(){
	return death;
}

int64_t WritePairs::getBirthX(){
	return birth_x;
}
int64_t WritePairs::getBirthY(){
	return birth_y;
}
int64_t WritePairs::getBirthZ(){
	return birth_z;
}
