/* compute_pairs.cpp

This file is part of CubicalRipser_3dim.
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <algorithm>
#include <queue>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>
#include <time.h>

using namespace std;

#include "cube.h"
#include "dense_cubical_grids.h"
#include "coboundary_enumerator.h"
#include "union_find.h"
#include "write_pairs.h"
#include "joint_pairs.h"
#include "compute_pairs.h"
	
ComputePairs::ComputePairs(DenseCubicalGrids* _dcg, vector<WritePairs> &_wp, const bool _print){
	dcg = _dcg;
	dim = 1; //  default method is LINK_FIND, where we skip dim=0
	wp = &_wp;
	print = _print;
}

void ComputePairs::compute_pairs_main(vector<Cube>& ctr, bool no_cache){
	if(print == true){
		cout << "persistence intervals in dim " << dim << ":" << endl;
	}
	
	pivot_column_index = std::unordered_map<int, int>();
	vector<Cube> coface_entries;
	auto ctl_size = ctr.size();
	CoboundaryEnumerator cofaces(dcg,dim);
	unordered_map<int, priority_queue<Cube, vector<Cube>, CubeComparator>> recorded_wc;

	pivot_column_index.reserve(ctl_size);
	recorded_wc.reserve(ctl_size);
		
	for(int i = 0; i < ctl_size; ++i){ 
		priority_queue<Cube, vector<Cube>, CubeComparator> working_coboundary;
		double birth = ctr[i].birthday;
		long idx = ctr[i].index;

		int j = i;
		Cube pivot(0, -1);
		bool might_be_apparent_pair = true;
		bool found_persistence_pair = false;

		do {
			coface_entries.clear();
			cofaces.setCoboundaryEnumerator(ctr[j]);

			while (cofaces.hasNextCoface() && !found_persistence_pair) { // repeat while there remains a coface
				coface_entries.push_back(cofaces.nextCoface);
				if (might_be_apparent_pair && (ctr[j].birthday == cofaces.nextCoface.birthday)) { 
					if (pivot_column_index.find(cofaces.nextCoface.index) == pivot_column_index.end()) { // If coface is not in pivot list
						pivot.copyCube(cofaces.nextCoface); // I have a new pivot
						found_persistence_pair = true;
					} else { // If pivot list contains this coface,
						might_be_apparent_pair = false;
					}
				}
			}

			if (found_persistence_pair) { 
				double death = pivot.birthday;
				if (birth != death) {
					vector<int> loc(dcg->getXYZM(idx));
					wp->push_back(WritePairs(dim, birth, death, loc[0], loc[1], loc[2], print));
				}
				pivot_column_index.emplace(pivot.index, i);
				break;
			}else{
				auto findWc = recorded_wc.find(j); 

				if(findWc != recorded_wc.end()){ // If the pivot is old,
					auto wc = findWc -> second;
					while(!wc.empty()){ // push the old pivot's wc
						working_coboundary.push(wc.top());
						wc.pop();
					}
				} else { // If the pivot is new,
					for(auto e : coface_entries){
						working_coboundary.push(e);
					}
				}
				pivot = get_pivot(working_coboundary); // get a pivot from wc

				if (pivot.index != -1) { // When I have a pivot, ...
					auto pair = pivot_column_index.find(pivot.index);
					if (pair != pivot_column_index.end()) {	// If the pivot already exists, go on the loop 
						j = pair -> second;
						continue;
					} else { // If the pivot is new
						if(!no_cache){
							recorded_wc.emplace(i, working_coboundary);
						}
						double death = pivot.birthday;
						if (birth != death) {
							vector<int> loc(dcg->getXYZM(idx));
							wp->push_back(WritePairs(dim, birth, death, loc[0], loc[1], loc[2], print));
						}
						pivot_column_index.emplace(pivot.index, i);
						break;
					}
				} else { // If wc is empty
					if (birth != dcg->threshold) {
						vector<int> loc(dcg->getXYZM(idx));
						wp->push_back(WritePairs(dim, birth, dcg->threshold, loc[0], loc[1], loc[2], print));
					}
					break;
				}
			}			

		} while (true);
	}
}


Cube ComputePairs::pop_pivot(priority_queue<Cube, vector<Cube>, CubeComparator>&
	column){
	if (column.empty()) {
		return Cube(0, -1);
	} else {
		auto pivot = column.top();
		column.pop();

		while (!column.empty() && column.top().index == pivot.index) {
			column.pop();
			if (column.empty())
				return Cube(0, -1);
			else {
				pivot = column.top();
				column.pop();
			}
		}
		return pivot;
	}
}

Cube ComputePairs::get_pivot(priority_queue<Cube, vector<Cube>, CubeComparator>&
	column) {
	Cube result = pop_pivot(column);
	if (result.index != -1) {
		column.push(result);
	}
	return result;
}

void ComputePairs::assemble_columns_to_reduce(vector<Cube>& ctr, int _dim) {
	dim = _dim;
	ctr.clear();
	double birthday;
	long ind;
	if (dim == 0) {
		for (int z = 0; z < dcg->az ; ++z) {
			for (int y = 0; y < dcg->ay; ++y) {
				for (int x = 0; x < dcg->ax; ++x) {
					birthday = dcg->getBirthday(x,y,z,0,0);
					if (birthday < dcg->threshold) {
						ind = dcg->getIndex(x, y, z, 0);
						ctr.push_back(Cube(birthday, ind));
					}
				}
			}
		}
	}else{ 
		for (int m = 0; m < 3; ++m) {
			for(int z = 0; z < dcg->az; ++z){
				for (int y = 0; y < dcg->ay; ++y) {
					for (int x = 0; x < dcg->ax; ++x) {
						ind = dcg->getIndex(x,y,z,m);
						if (pivot_column_index.find(ind) == pivot_column_index.end()) {
							birthday = dcg -> getBirthday(x,y,z,m, dim);
							if (birthday < dcg -> threshold) {
								ctr.push_back(Cube(birthday, ind));
							}
						}
					}
				}
			}
		}
	}
    clock_t start = clock();
    sort(ctr.begin(), ctr.end(), CubeComparator());
    clock_t end = clock();
    const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000.0;
    cout << "Sorting Time: " <<  time << endl;
}
