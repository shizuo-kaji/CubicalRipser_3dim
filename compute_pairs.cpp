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

using namespace std;

#include "birthday_index.h"
#include "dense_cubical_grids.h"
#include "simplex_coboundary_enumerator.h"
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

void ComputePairs::compute_pairs_main(vector<BirthdayIndex>& ctr){
	if(print == true){
		cout << "persistence intervals in dim " << dim << ":" << endl;
	}
	
	pivot_column_index = std::unordered_map<int, int>();
	vector<BirthdayIndex> coface_entries;
	auto ctl_size = ctr.size();
	SimplexCoboundaryEnumerator cofaces(dcg);
	unordered_map<int, priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>> recorded_wc;

	pivot_column_index.reserve(ctl_size);
	recorded_wc.reserve(ctl_size);
		
	for(int i = 0; i < ctl_size; ++i){ 
		priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator> working_coboundary;
		double birth = ctr[i].getBirthday();
		long idx = ctr[i].getIndex();

		int j = i;
		BirthdayIndex pivot(0, -1, 0);
		bool might_be_apparent_pair = true;
		bool goto_found_persistence_pair = false;

		do {
			coface_entries.clear();
			cofaces.setSimplexCoboundaryEnumerator(ctr[j]);

			while (cofaces.hasNextCoface() && !goto_found_persistence_pair) { // repeat while there remains a coface
				BirthdayIndex coface = cofaces.getNextCoface();
				coface_entries.push_back(coface);
				if (might_be_apparent_pair && (ctr[j].getBirthday() == coface.getBirthday())) { 
					if (pivot_column_index.find(coface.getIndex()) == pivot_column_index.end()) { // If coface is not in pivot list
						pivot.copyBirthdayIndex(coface); // I have a new pivot
						goto_found_persistence_pair = true; // goto (B)
					} else { // If pivot list contains this coface,
						might_be_apparent_pair = false; // goto (A)
					}
				}
			}

			if (!goto_found_persistence_pair) { // (A) If pivot list contains this coface,
				auto findWc = recorded_wc.find(j); 

				if(findWc != recorded_wc.end()){ // If the pivot is old,
					auto wc = findWc -> second;
					while(!wc.empty()){ // push the data of the old pivot's wc
						auto e = wc.top(); // TODO: use list of pointers to avoid (de)construction
						working_coboundary.push(e);
						wc.pop();
					}
				} else { // If the pivot is new,
					for(auto e : coface_entries){
						working_coboundary.push(e);
					}
				}
				pivot = get_pivot(working_coboundary); // get a pivot from wc

				if (pivot.getIndex() != -1) { // When I have a pivot, ...
					auto pair = pivot_column_index.find(pivot.getIndex());
					if (pair != pivot_column_index.end()) {	// If the pivot already exists, go on the loop 
						j = pair -> second;
						continue;
					} else { // If the pivot is new, 
						// record this wc into recorded_wc, and 
						recorded_wc.insert(make_pair(i, working_coboundary));
						double death = pivot.getBirthday();
						if (birth != death) {
							vector<int> loc(dcg->getXYZM(idx));
							wp->push_back(WritePairs(dim, birth, death, loc[0], loc[1], loc[2], print));
						}
						pivot_column_index.insert(make_pair(pivot.getIndex(), i));
						break;
					}
				} else { // If wc is empty
					if (birth != dcg->threshold) {
						vector<int> loc(dcg->getXYZM(idx));
						wp->push_back(WritePairs(0, birth, dcg->threshold, loc[0], loc[1], loc[2], print));
					}
					break;
				}
			} else { // (B) I have a new pivot
				double death = pivot.getBirthday();
				if (birth != death) {
					vector<int> loc(dcg->getXYZM(idx));
					wp->push_back(WritePairs(dim, birth, death, loc[0], loc[1], loc[2], print));
				}
				pivot_column_index.insert(make_pair(pivot.getIndex(), i));
				break;
			}			

		} while (true);
	}
}


BirthdayIndex ComputePairs::pop_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
	column){
	if (column.empty()) {
		return BirthdayIndex(0, -1, 0);
	} else {
		auto pivot = column.top();
		column.pop();

		while (!column.empty() && column.top().index == pivot.getIndex()) {
			column.pop();
			if (column.empty())
				return BirthdayIndex(0, -1, 0);
			else {
				pivot = column.top();
				column.pop();
			}
		}
		return pivot;
	}
}

BirthdayIndex ComputePairs::get_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
	column) {
	BirthdayIndex result = pop_pivot(column);
	if (result.getIndex() != -1) {
		column.push(result);
	}
	return result;
}

void ComputePairs::assemble_columns_to_reduce(vector<BirthdayIndex>& ctr, int _dim) {
	dim = _dim;
	ctr.clear();
	double birthday;
	long ind;
	if (dim == 0) {
		for (int z = 0; z < dcg->az ; ++z) {
			for (int y = 0; y < dcg->ay; ++y) {
				for (int x = 0; x < dcg->ax; ++x) {
					birthday = dcg->get(x, y, z);
					if (birthday < dcg->threshold) {
						ind = dcg->getIndex(x, y, z, 0);
						ctr.push_back(BirthdayIndex(birthday, ind, 0));
					}
				}
			}
		}
	}else{ 
		for(int z = 0; z < dcg->az; ++z){
			for (int y = 0; y < dcg->ay; ++y) {
				for (int x = 0; x < dcg->ax; ++x) {
					for (int m = 0; m < 3; ++m) {
						ind = dcg->getIndex(x,y,z,m);
						if (pivot_column_index.find(ind) == pivot_column_index.end()) {
							birthday = dcg -> getBirthday(x,y,z,m, dim);
							if (birthday < dcg -> threshold) {
								ctr.push_back(BirthdayIndex(birthday, ind, dim));
							}
						}
					}
				}
			}
		}
	}
	sort(ctr.begin(), ctr.end(), BirthdayIndexComparator());
}