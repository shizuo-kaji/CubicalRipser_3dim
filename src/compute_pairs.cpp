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
//	pivot_column_index.set_empty_key(0xffffffff); // for googlehash
}

void ComputePairs::compute_pairs_main(vector<Cube>& ctr){
	vector<Cube> coface_entries; // pivotIDs of cofaces
	auto ctl_size = ctr.size();
	pivot_column_index.clear();
	pivot_column_index.reserve(ctl_size);
//	pivot_column_index.resize(ctl_size); // googlehash
	CoboundaryEnumerator cofaces(dcg,dim);
	unordered_map<uint32_t, CubeQue > recorded_wc;
//	unordered_map<int, vector<Cube> > recorded_wc;   // for some unknown reasons, using vector directly with sorting when needed deteriorates performance.
	queue<uint32_t> cached_column_idx;

	recorded_wc.reserve(ctl_size);

	for(uint32_t i = 0; i < ctl_size; ++i){  // descending order of birth
		CubeQue working_coboundary;   // non-zero entries of the column
        //    vector<Cube> working_coboundary;
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
                    auto wc = findWc -> second;
                    while(!wc.empty()){ // add the cached pivot column
                        working_coboundary.push(wc.top());
                        wc.pop();
                    }
                    cache_hit = true;
                    /// If we use vector container
                    // for(auto e:wc){
                    //     auto idx = std::find(working_coboundary.begin(),working_coboundary.end(),e);
                    //     if(idx==working_coboundary.end()){
                    //         working_coboundary.push_back(e);
                    //     }else{
                    //         working_coboundary.erase(idx);
                    //     }
                    // }
                    /// working_coboundary.insert(working_coboundary.end(), wc.begin(), wc.end());
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
                    double death = pivot.birth;
                    if (birth != death) {
//						wp->push_back(WritePairs(dim, birth, death, ctr[i].x(), ctr[i].y(), ctr[i].z(), pivot.x(), pivot.y(), pivot.z(), config->print));
						wp->push_back(WritePairs(dim, ctr[i], pivot, dcg, config->print));
                    }
    //                cout << pivot.index << ",ap," << i << endl;
                    //pivot_column_index.emplace(pivot.index, i);
                    pivot_column_index[pivot.index] = i;
                    break;
                }
                for(auto e : coface_entries){
                    working_coboundary.push(e);
                }
                /// If we use vector container
                // for(auto e:coface_entries){
                //     auto idx = std::find(working_coboundary.begin(),working_coboundary.end(),e);
                //     if(idx==working_coboundary.end()){
                //         working_coboundary.push_back(e);
                //     }else{
                //         working_coboundary.erase(idx);
                //     }
                // }
                /// working_coboundary.insert(working_coboundary.end(), coface_entries.begin(), coface_entries.end());
            }
            
            /// If we use vector container
            // sort(working_coboundary.begin(),working_coboundary.end(),CubeComparator());
            pivot = get_pivot(working_coboundary);
            if (pivot.index != NONE){
//                if(!working_coboundary.empty()){
//                    pivot = working_coboundary.back();
                auto pair = pivot_column_index.find(pivot.index);
                if (pair != pivot_column_index.end()) {	// found entry to reduce
                    j = pair -> second;
					num_recurse++;
//                        cout << i << " to " << j << " " << pivot.index << endl;
                    continue;
                } else { // If the pivot is new
//                    if((int)working_coboundary.size() >= config->min_cache_size){
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
//						wp->push_back(WritePairs(dim, birth, death, ctr[i].x(), ctr[i].y(), ctr[i].z(), pivot.x(), pivot.y(), pivot.z(), config->print));
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

		};
	}
}

// cache newly found pivot after reducing by mod 2
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

// mod 2 operation
Cube ComputePairs::pop_pivot(vector<Cube>& column){
	if (column.empty()) {
		return Cube();
	} else {
		auto pivot = column.back();
		column.pop_back();
		while (!column.empty() && column.back().index == pivot.index) {
			column.pop_back();
			if (column.empty())
				return Cube();
			else {
				pivot = column.back();
			column.pop_back();
			}
		}
		return pivot;
	}
}

Cube ComputePairs::get_pivot(vector<Cube>& column) {
	Cube result = pop_pivot(column);
	if (result.index != NONE) {
		column.push_back(result);
	}
	return result;
}

// the same mod 2 operation for vector
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
	if (dim == 0) {
		for (uint32_t z = 0; z < dcg->az ; ++z) {
			for (uint32_t y = 0; y < dcg->ay; ++y) {
				for (uint32_t x = 0; x < dcg->ax; ++x) {
					birth = dcg->getBirth(x,y,z,0,0);
					if (birth < dcg->threshold) {
						ctr.push_back(Cube(birth, x,y,z,0));
					}
				}
			}
		}
	}else{ 
		for (uint8_t m = 0; m < 3; ++m) {
			for(uint32_t z = 0; z < dcg->az; ++z){
				for (uint32_t y = 0; y < dcg->ay; ++y) {
					for (uint32_t x = 0; x < dcg->ax; ++x) {
						birth = dcg -> getBirth(x,y,z,m, dim);
//                        cout << x << "," << y << "," << z << ", " << m << "," << birth << endl;
						Cube v = Cube(birth,x,y,z,m);
						if (pivot_column_index.find(v.index) == pivot_column_index.end()) {
							if (birth < dcg -> threshold) {
								ctr.push_back(v);
							}
						}
					}
				}
			}
		}
	}
//    clock_t start = clock();
    sort(ctr.begin(), ctr.end(), CubeComparator());
//    clock_t end = clock();
//    const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000.0;
//    cout << "Sorting Time: " <<  time << endl;
}
