/* simplex_coboundary_enumerator.cpp

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


#include <algorithm>
#include <vector>

#include "birthday_index.h"
#include "dense_cubical_grids.h"
#include "columns_to_reduce.h"
#include "simplex_coboundary_enumerator.h"

using namespace std;

SimplexCoboundaryEnumerator::SimplexCoboundaryEnumerator(){
	vtx = new Vertices();
	nextCoface = BirthdayIndex(0, -1, 1);
}
	

void SimplexCoboundaryEnumerator::setSimplexCoboundaryEnumerator(BirthdayIndex _s, DenseCubicalGrids* _dcg) {
	simplex = _s;
	dcg = _dcg;
	_dcg -> GetSimplexVertices(simplex.index, simplex.dim, vtx);
	birthtime = simplex.birthday;
	ax = _dcg -> ax;
	ay = _dcg -> ay;
	az = _dcg -> az;
		
	threshold = _dcg -> threshold;
	count = 0;
}


bool SimplexCoboundaryEnumerator::hasNextCoface() {
	double birthday = 0;
	int m = 0;
	cx = vtx -> ox;
	cy = vtx -> oy;
	cz = vtx -> oz;
	ArrayIndex index;
	switch (vtx->dim) {
		case 0: // dim0
		for (int i = count; i < 6; ++i) {
			switch (i){
				case 0:
					index.set(cx, cy, cz, 2);
					birthday = max(birthtime, dcg -> dense3[cx][cy][cz + 1]);
				break;

				case 1:
					index.set(cx, cy, cz-1, 2);
					birthday = max(birthtime, dcg -> dense3[cx][cy][cz - 1]);
				break;

				case 2:
					index.set(cx, cy, cz, 1);
					birthday = max(birthtime, dcg -> dense3[cx][cy + 1][cz]);
				break;

				case 3:
					index.set(cx, cy-1, cz, 1);
					birthday = max(birthtime, dcg -> dense3[cx][cy - 1][cz]);
				break;

				case 4:
					index.set(cx, cy, cz, 0);
					birthday = max(birthtime, dcg -> dense3[cx + 1][cy][cz]);
				break;

				case 5:
					index.set(cx-1, cy, cz, 0);
					birthday = max(birthtime, dcg -> dense3[cx - 1][cy][cz]);
				break;
			}
			if (birthday != threshold) {
				count = i + 1;
				nextCoface = BirthdayIndex(birthday, index.getIndex(), 1);
				return true;
			}
		}
		return false;

		case 1: // dim1
		switch (vtx->type) {
			case 0: // dim1 type0 (x-axis -> )
			for(int i = count; i < 4; ++i){
				switch(i){
					case 0:
						index.set(cx, cy, cz, 1);
						birthday = max({birthtime, dcg -> dense3[cx][cy][cz + 1], dcg -> dense3[cx + 1][cy][cz + 1]});
					break;

					case 1:
						index.set(cx, cy, cz-1, 1);
						birthday = max({birthtime, dcg -> dense3[cx][cy][cz - 1], dcg -> dense3[cx + 1][cy][cz - 1]});
					break;

					case 2:
						index.set(cx, cy, cz, 0);
						birthday = max({birthtime, dcg -> dense3[cx][cy + 1][cz], dcg -> dense3[cx + 1][cy + 1][cz]});
					break;

					case 3:
						index.set(cx, cy-1, cz, 0);
						birthday = max({birthtime, dcg -> dense3[cx][cy - 1][cz], dcg -> dense3[cx + 1][cy - 1][cz]});
					break;
				}

				if (birthday != threshold) {
					count = i + 1;
					nextCoface = BirthdayIndex(birthday, index.getIndex(), 2);
					return true;
				}
			}
			return false;

			case 1: // dim1 type1 (y-axis -> )
			for(int i = count; i < 4; ++i){
				switch(i){
					case 0:
						index.set(cx, cy, cz, 2);
						birthday = max({birthtime, dcg -> dense3[cx][cy][cz + 1], dcg -> dense3[cx][cy + 1][cz + 1]});
					break;

					case 1:
						index.set(cx, cy, cz-1, 2);
						birthday = max({birthtime, dcg -> dense3[cx][cy][cz - 1], dcg -> dense3[cx][cy + 1][cz - 1]});
					break;

					case 2:
						index.set(cx, cy, cz, 0);
						birthday = max({birthtime, dcg -> dense3[cx + 1][cy][cz], dcg -> dense3[cx + 1][cy + 1][cz]});
					break;

					case 3:
						index.set(cx-1, cy, cz, 0);
						birthday = max({birthtime, dcg -> dense3[cx - 1][cy][cz], dcg -> dense3[cx - 1][cy + 1][cz]});
					break;
				}
				if (birthday != threshold) {
					count = i + 1;
					nextCoface = BirthdayIndex(birthday, index.getIndex(), 2);
					return true;
				}
			}
			return false;

			case 2: // dim1 type2 (z-axis -> )
			for(int i = count; i < 4; ++i){
				switch(i){
					case 0:
						index.set(cx, cy, cz, 2);
						birthday = max({birthtime, dcg -> dense3[cx][cy + 1][cz], dcg -> dense3[cx][cy + 1][cz + 1]});
					break;

					case 1:
						index.set(cx, cy-1, cz, 2);
						birthday = max({birthtime, dcg -> dense3[cx][cy - 1][cz], dcg -> dense3[cx][cy - 1][cz + 1]});
					break;

					case 2:
						index.set(cx, cy, cz, 1);
						birthday = max({birthtime, dcg -> dense3[cx + 1][cy][cz], dcg -> dense3[cx + 1][cy][cz + 1]});
					break;

					case 3:
						index.set(cx-1, cy, cz, 1);
						birthday = max({birthtime, dcg -> dense3[cx - 1][cy][cz], dcg -> dense3[cx - 1][cy][cz + 1]});
					break;
				}
				if (birthday != threshold) {
					count = i + 1;
					nextCoface = BirthdayIndex(birthday, index.getIndex(), 2);
					return true;
				}
			}
			return false;
		}
		return false;

		default: // dim2
		switch (vtx->type) {
			case 0: // dim2 type0 (fix z)
			for(int i = count; i < 2; ++i){
				switch(i){
				case 0: // upper
					index.set(cx, cy, cz, 0);
					birthday = max({birthtime, dcg -> dense3[cx][cy][cz + 1], dcg -> dense3[cx + 1][cy][cz + 1],
					dcg -> dense3[cx][cy + 1][cz + 1],dcg -> dense3[cx + 1][cy + 1][cz + 1]});
				break;

				case 1: // lower
					index.set(cx, cy, cz-1, 0);
					birthday = max({birthtime, dcg -> dense3[cx][cy][cz - 1], dcg -> dense3[cx + 1][cy][cz - 1],
					dcg -> dense3[cx][cy + 1][cz - 1],dcg -> dense3[cx + 1][cy + 1][cz - 1]});
				break;

				}
				if (birthday != threshold) {
					count = i + 1;
					nextCoface = BirthdayIndex(birthday, index.getIndex(), 3);
					return true;
				}
			}
			return false;

			case 1: // dim2 type1 (fix y)
			for(int i = count; i < 2; ++i){
				switch(i){
				case 0: // left
					index.set(cx, cy, cz, 0);
					birthday = max({birthtime, dcg -> dense3[cx][cy + 1][cz], dcg -> dense3[cx + 1][cy + 1][cz],
					dcg -> dense3[cx][cy + 1][cz + 1],dcg -> dense3[cx + 1][cy + 1][cz + 1]});
				break;

				case 1: //right
					index.set(cx, cy-1, cz, 0);
					birthday = max({birthtime, dcg -> dense3[cx][cy - 1][cz], dcg -> dense3[cx + 1][cy - 1][cz],
					dcg -> dense3[cx][cy - 1][cz + 1],dcg -> dense3[cx + 1][cy - 1][cz + 1]});
				break;

				}
				if (birthday != threshold) {
					count = i + 1;
					nextCoface = BirthdayIndex(birthday, index.getIndex(), 3);
					return true;
				}
			}
			return false;

			case 2: // dim2 type2 (fix x)
			for(int i = count; i < 2; ++i){
				switch(i){
				case 0: // left
					index.set(cx, cy, cz, 0);
					birthday = max({birthtime, dcg -> dense3[cx + 1][cy][cz], dcg -> dense3[cx + 1][cy + 1][cz],
					dcg -> dense3[cx + 1][cy][cz + 1],dcg -> dense3[cx + 1][cy + 1][cz + 1]});
				break;

				case 1: //right
					index.set(cx-1, cy, cz, 0);
					birthday = max({birthtime, dcg -> dense3[cx - 1][cy][cz], dcg -> dense3[cx - 1][cy + 1][cz],
					dcg -> dense3[cx - 1][cy][cz + 1],dcg -> dense3[cx - 1][cy + 1][cz + 1]});
				break;

				}
				if (birthday != threshold) {
					count = i + 1;
					nextCoface = BirthdayIndex(birthday, index.getIndex(), 3);
					return true;
				}
			}
			return false;
		}
		return false;
	}
}

BirthdayIndex SimplexCoboundaryEnumerator::getNextCoface() {
	return nextCoface;
}
