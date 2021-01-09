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
#include "write_pairs.h"
#include "compute_pairs.h"


ComputePairs::ComputePairs(DenseCubicalGrids* _dcg, vector<WritePairs> &_wp, Config& _config){
	dcg = _dcg;
	dim = 1; //  default method is LINK_FIND, where we skip dim=0
	wp = &_wp;
	config = &_config;
#ifdef GOOGLE_HASH
	pivot_column_index.set_empty_key(0xffffffff); // for googlehash
#endif
}

void ComputePairs::compute_pairs_main(vector<Cube>& ctr){
	vector<Cube> coface_entries; // pivotIDs of cofaces
	auto ctl_size = ctr.size();
	if(config->verbose){
	    cout << "# columns to reduce: " << ctl_size << endl;
	}
	pivot_column_index.clear();
#ifdef GOOGLE_HASH
	pivot_column_index.resize(ctl_size); // googlehash
#else
    pivot_column_index.reserve(ctl_size);
#endif
	CoboundaryEnumerator cofaces(dcg,dim);
	unordered_map<uint32_t, CubeQue > recorded_wc;
	queue<uint32_t> cached_column_idx;
	recorded_wc.reserve(ctl_size);
    int num_apparent_pairs = 0;

	for(uint32_t i = 0; i < ctl_size; ++i){  // descending order of birth
        CubeQue working_coboundary;   // non-zero entries of the column
		double birth = ctr[i].birth;
//        cout << i << endl;  ctr[i].print();   // debug

		auto j = i;
		Cube pivot;
		bool might_be_apparent_pair = true;
		bool found_persistence_pair = false;
		int num_recurse = 0;

		while(true){
            bool cache_hit = false;
            if(i!=j){
                auto findWc = recorded_wc.find(j);
                if(findWc != recorded_wc.end()){ // If the reduced form of the pivot column is cached
                    cache_hit = true;
                    auto wc = findWc -> second;
                    while(!wc.empty()){ // add the cached pivot column
                        working_coboundary.push(wc.top());
                        wc.pop();
                    }
                }
//				assert(might_be_apparent_pair == false); // As there is always cell-coface pair with the same birthtime, the flag should be set by the next block.
			}
            if(!cache_hit){
                // make the column by enumerating cofaces
                coface_entries.clear();
                cofaces.setCoboundaryEnumerator(ctr[j]);
                while (cofaces.hasNextCoface()) {
                    coface_entries.push_back(cofaces.nextCoface);
    //                cout << "cf: " << j << endl;
//                    cofaces.nextCoface.print();
                    if (might_be_apparent_pair && (ctr[j].birth == cofaces.nextCoface.birth)) { // we cannot find this coface on the left (Short-Circuit Evaluation)
                        if (pivot_column_index.find(cofaces.nextCoface.index) == pivot_column_index.end()) { // If coface is not in pivot list
                            pivot.copyCube(cofaces.nextCoface);
                            found_persistence_pair = true;
							break;
                        } else {
                            might_be_apparent_pair = false;
                        }
                    }
                }
                if (found_persistence_pair) {
                    //pivot_column_index.emplace(pivot.index, i);
                    pivot_column_index[pivot.index] = i;
                    num_apparent_pairs++;
                    break;
                }
                for(auto e : coface_entries){
                    working_coboundary.push(e);
                }
            }
            pivot = get_pivot(working_coboundary);
            if (pivot.index != NONE){ // if the column is not reduced to zero
                auto pair = pivot_column_index.find(pivot.index);
                if (pair != pivot_column_index.end()) {	// found entry to reduce
                    j = pair -> second;
					num_recurse++;
//                        cout << i << " to " << j << " " << pivot.index << endl;
                    continue;
                } else { // If the pivot is new
                    if(num_recurse >= config->min_recursion_to_cache){
                        add_cache(i, working_coboundary, recorded_wc);
						cached_column_idx.push(i);
						if(cached_column_idx.size()>config->cache_size){
							recorded_wc.erase(cached_column_idx.front());
							cached_column_idx.pop();
						}
                    }
                    // pivot_column_index.emplace(pivot.index, i); // column i has the pivot
					pivot_column_index[pivot.index] = i;
                    double death = pivot.birth;
                    if (birth != death) {
						wp->push_back(WritePairs(dim, ctr[i], pivot, dcg, config->print));
                    }
//                        cout << pivot.index << ",f," << i << endl;
                    break;
                }
            } else { // the column is reduced to zero, which means it corresponds to a permanent cycle
                if (birth != dcg->threshold) {
					wp->push_back(WritePairs(dim, birth, dcg->threshold, ctr[i].x(), ctr[i].y(), ctr[i].z(),  0, 0, 0, config->print));
                }
                break;
            }
		}
	}
    if(config->verbose){
        cout << "# apparent pairs: " << num_apparent_pairs << endl;
    }
}

// cache a new reduced column after mod 2
void ComputePairs::add_cache(uint32_t i, CubeQue &wc, unordered_map<uint32_t, CubeQue>& recorded_wc){
	CubeQue clean_wc;
	while(!wc.empty()){
		auto c = wc.top();
		wc.pop();
		if(!wc.empty() && c.index==wc.top().index){
			wc.pop();
		}else{
			clean_wc.push(c);
		}
	}
	recorded_wc.emplace(i,clean_wc);
}

// get the pivot from a column after mod 2
Cube ComputePairs::pop_pivot(CubeQue& column){
    if (column.empty()) {
        return Cube();
    } else {
        auto pivot = column.top();
        column.pop();

        while (!column.empty() && column.top().index == pivot.index) {
            column.pop();
            if (column.empty())
                return Cube();
            else {
                pivot = column.top();
                column.pop();
            }
        }
        return pivot;
    }
}

Cube ComputePairs::get_pivot(CubeQue& column) {
	Cube result = pop_pivot(column);
	if (result.index != NONE) {
		column.push(result);
	}
	return result;
}

// enumerate and sort columns for a new dimension
void ComputePairs::assemble_columns_to_reduce(vector<Cube>& ctr, uint8_t _dim) {
	dim = _dim;
	ctr.clear();
	double birth;
    uint8_t max_m = 3;
	if (dim == 0) {
        max_m = 1;
        pivot_column_index.clear();
    }
    for (uint8_t m = 0; m < max_m; ++m) {
        for(uint32_t z = 0; z < dcg->az; ++z){
            for (uint32_t y = 0; y < dcg->ay; ++y) {
                for (uint32_t x = 0; x < dcg->ax; ++x) {
                    birth = dcg -> getBirth(x,y,z,m, dim);
//                        cout << x << "," << y << "," << z << ", " << m << "," << birth << endl;
                    Cube v = Cube(birth,x,y,z,m);
                    if (birth < dcg -> threshold && pivot_column_index.find(v.index) == pivot_column_index.end()) {
                        ctr.push_back(v);
                    }
                }
            }
        }
    }
    clock_t start = clock();
    sort(ctr.begin(), ctr.end(), CubeComparator());
	if(config->verbose){
		clock_t end = clock();
		const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000.0;
		cout << "Sorting took: " <<  time << endl;
	}
}
