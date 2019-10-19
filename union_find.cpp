/* union_find.cpp

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


#include <vector>
#include <algorithm>
#include "birthday_index.h"
#include "dense_cubical_grids.h"
#include "union_find.h"

using namespace std;

UnionFind::UnionFind(int moi, DenseCubicalGrids* _dcg) : parent(moi), birthtime(moi), time_max(moi) {
	dcg = _dcg;
	max_of_index = moi;

	for(int i = 0; i < moi; ++i){
		parent[i] = i;
		birthtime[i] = dcg -> getBirthday(i, 0);
		time_max[i] = dcg -> getBirthday(i, 0);
	}
}

// find the root of a node x (specified by the index)
int UnionFind::find(int x){
	int y = x, z = parent[y];
	while (z != y) {
		y = z;
		z = parent[y];
	}
	// reassign parents to the found root z
	y = parent[x];
	while (z != y) {
		parent[x] = z;
		x = y;
		y = parent[x];
	}
	return z;
}

// merge nodes x and y
void UnionFind::link(int x, int y){
	x = find(x);
	y = find(y);
	if (x == y) return;
	if (birthtime[x] > birthtime[y]){
		parent[x] = y; 
		time_max[y] = std::max(time_max[x], time_max[y]);
	} else if(birthtime[x] < birthtime[y]) {
		parent[y] = x;
		time_max[x] = std::max(time_max[x], time_max[y]);
	} else { //birthtime[x] == birthtime[y]
		parent[x] = y;
		time_max[y] = std::max(time_max[x], time_max[y]);
	}
}