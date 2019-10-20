/* vertices.cpp

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


#include "vertices.h"

Vertices::Vertices(){
	dim = 0;
	vertex = new ArrayIndex*[8];
	for(int d = 0; d < 8; ++d){
		vertex[d] = new ArrayIndex();
	}
}
void Vertices::setVertices(int _dim, long idx) {
	ArrayIndex index(idx);
	setVertices(_dim, index.cx, index.cy, index.cz, index.cm);
}

void Vertices::setVertices(int _dim, int _ox, int _oy, int _oz, int _om){
	dim = _dim;
	ox = _ox;
	oy = _oy;
	oz = _oz;
	type = _om;

	if(dim == 0){
		vertex[0] -> set(_ox, _oy, _oz);
	} else if(dim == 1){
		switch(_om){
			case 0:
			vertex[0] -> set(_ox, _oy, _oz);
			vertex[1] -> set(_ox + 1, _oy, _oz);
			break;

			case 1:
			vertex[0] -> set(_ox, _oy, _oz);
			vertex[1] -> set(_ox, _oy + 1, _oz);
			break;

			default:
			vertex[0] -> set(_ox, _oy, _oz);
			vertex[1] -> set(_ox, _oy, _oz + 1);
			break;
		}

	} else if(dim == 2){
		switch(_om){
			case 0: // x - y
			vertex[0] -> set(_ox, _oy, _oz);
			vertex[1] -> set(_ox + 1, _oy, _oz);
			vertex[2] -> set(_ox + 1, _oy + 1, _oz);
			vertex[3] -> set(_ox, _oy + 1, _oz);
			break;

			case 1: // z - x
			vertex[0] -> set(_ox, _oy, _oz);
			vertex[1] -> set(_ox, _oy, _oz + 1);
			vertex[2] -> set(_ox + 1, _oy, _oz + 1);
			vertex[3] -> set(_ox + 1, _oy, _oz);
			break;

			default: // y - z
			vertex[0] -> set(_ox, _oy, _oz);
			vertex[1] -> set(_ox, _oy + 1, _oz);
			vertex[2] -> set(_ox, _oy + 1, _oz + 1);
			vertex[3] -> set(_ox, _oy, _oz + 1);
			break;
		}

	} else if(dim == 3){ // cube
		vertex[0] -> set(_ox, _oy, _oz);
		vertex[1] -> set(_ox + 1, _oy, _oz);
		vertex[2] -> set(_ox + 1, _oy + 1, _oz);
		vertex[3] -> set(_ox, _oy + 1, _oz);
		vertex[4] -> set(_ox, _oy, _oz + 1);
		vertex[5] -> set(_ox + 1, _oy, _oz + 1);
		vertex[6] -> set(_ox + 1, _oy + 1, _oz + 1);
		vertex[7] -> set(_ox, _oy + 1, _oz + 1);
	}
}
