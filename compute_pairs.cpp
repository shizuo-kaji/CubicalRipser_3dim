/* compute_pairs.cpp

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
#include "columns_to_reduce.h"
#include "simplex_coboundary_enumerator.h"
#include "union_find.h"
#include "write_pairs.h"
#include "joint_pairs.h"
#include "compute_pairs.h"
	
ComputePairs::ComputePairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp, const bool _print){
	dcg = _dcg;
	ctr = _ctr;
	dim = _ctr -> dim;
	wp = &_wp;
	print = _print;

	ax = _dcg -> ax;
	ay = _dcg -> ay;
	az = _dcg -> az;
}

void ComputePairs::compute_pairs_main(){
	if(print == true){
		cout << "persistence intervals in dim " << dim << ":" << endl;
	}
	
	pivot_column_index = hash_map<int, int>();
	vector<BirthdayIndex> coface_entries;
	auto ctl_size = ctr -> columns_to_reduce.size();
	SimplexCoboundaryEnumerator cofaces;
	unordered_map<int, priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>> recorded_wc;

	pivot_column_index.reserve(ctl_size);
	recorded_wc.reserve(ctl_size);
		
	for(int i = 0; i < ctl_size; ++i){ 
		auto column_to_reduce = ctr -> columns_to_reduce[i]; 
		priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator> 
		working_coboundary;
		double birth = column_to_reduce.getBirthday();
		long idx = column_to_reduce.getIndex();

		int j = i;
		BirthdayIndex pivot(0, -1, 0);
		bool might_be_apparent_pair = true;
		bool goto_found_persistence_pair = false;

		do {
			auto simplex = ctr -> columns_to_reduce[j]; // get CTR[i] 
			coface_entries.clear();
			cofaces.setSimplexCoboundaryEnumerator(simplex, dcg);// make coface data

			while (cofaces.hasNextCoface() && !goto_found_persistence_pair) { // repeat there remains a coface
				BirthdayIndex coface = cofaces.getNextCoface();
				coface_entries.push_back(coface);
				if (might_be_apparent_pair && (simplex.getBirthday() == coface.getBirthday())) { // If bt is the same, go thru
					if (pivot_column_index.find(coface.getIndex()) == pivot_column_index.end()) { // If coface is not in pivot list
						pivot.copyBirthdayIndex(coface); // I have a new pivot
						goto_found_persistence_pair = true; // goto (B)
					} else { // If pivot list contains this coface,
						might_be_apparent_pair = false; // goto (A)
					}
				}
			}

			if (!goto_found_persistence_pair) { // (A) If pivot list contains this coface,
				auto findWc = recorded_wc.find(j); // we seek wc list by 'j'

				if(findWc != recorded_wc.end()){ // If the pivot is old,
					auto wc = findWc -> second;
					while(!wc.empty()){ // we push the data of the old pivot's wc
						auto e = wc.top();
						working_coboundary.push(e);
						wc.pop();
					}
				} else { // If the pivot is new,
					for(auto e : coface_entries){ // making wc here
						working_coboundary.push(e);
					}
				}
				pivot = get_pivot(working_coboundary); // getting a pivot from wc

				if (pivot.getIndex() != -1) { // When I have a pivot, ...
					auto pair = pivot_column_index.find(pivot.getIndex());
					if (pair != pivot_column_index.end()) {	// If the pivot already exists, go on the loop 
						j = pair -> second;
						continue;
					} else { // If the pivot is new, 
						// I record this wc into recorded_wc, and 
						recorded_wc.insert(make_pair(i, working_coboundary));
						// I output PP as Writepairs
						double death = pivot.getBirthday();
						outputPP(dim, birth, death, idx);
						pivot_column_index.insert(make_pair(pivot.getIndex(), i));
						break;
					}
				} else { // If wc is empty, I output a PP as [birth,) 
					outputPP(-1, birth, dcg -> threshold, idx);
					break;
				}
			} else { // (B) I have a new pivot and output PP as Writepairs 
				double death = pivot.getBirthday();
				outputPP(dim, birth, death, idx);
				pivot_column_index.insert(make_pair(pivot.getIndex(), i));
				break;
			}			

		} while (true);
	}
}

void ComputePairs::outputPP(int _dim, double _birth, double _death, int idx){
	if(_birth != _death){
		int x = idx & 0x01ff;
		int y = (idx >> 9) & 0x01ff;
		int z = (idx >> 18) & 0x01ff;
		if(_death != dcg -> threshold){
			if(print == true){
				cout << "[" << _birth << "," << _death << ")" << " birth loc (" << x << "," << y << "," << z << ")" << endl;
			}
			wp -> push_back(WritePairs(_dim, _birth, _death,x,y,z));
		} else {
			if(print == true){
				cout << "[" << _birth << ", )" << " birth loc (" << x << "," << y << "," << z << ")" << endl;
			}
			wp -> push_back(WritePairs(_dim, _birth, _death,x,y,z));
		}
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

void ComputePairs::assemble_columns_to_reduce() {
	++dim;
	ctr -> dim = dim;

	if (dim == 1) { 
		ctr -> columns_to_reduce.clear();
		for(int z = 1; z <= az; ++z){
			for (int y = 1; y <= ay; ++y) {
				for (int x = 1; x <= ax; ++x) {
					for (int m = 0; m < 3; ++m) { // the number of type
						double index = x | (y << 9) | (z << 18) | (m << 27);
						if (pivot_column_index.find(index) == pivot_column_index.end()) {
							double birthday = dcg -> getBirthday(index, 1);
							if (birthday != dcg -> threshold) {
								ctr -> columns_to_reduce.push_back(BirthdayIndex(birthday, index, 1));
							}
						}
					}
				}
			}
		}
	} else if(dim == 2){ 
		ctr -> columns_to_reduce.clear();
		for(int z = 1; z <= az; ++z){
			for (int y = 1; y <= ay; ++y) {
				for (int x = 1; x <= ax; ++x) {
					for (int m = 0; m < 3; ++m) { // the number of type
						double index = x | (y << 9) | (z << 18) | (m << 27);
						if (pivot_column_index.find(index) == pivot_column_index.end()) {
							double birthday = dcg -> getBirthday(index, 2);
							if (birthday != dcg -> threshold) {
								ctr -> columns_to_reduce.push_back(BirthdayIndex(birthday, index, 2));
							}
						}
					}
				}
			}
		}
	}
	sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndexComparator());
}