/* birthday_index.cpp

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include "birthday_index.h"

using namespace std;

BirthdayIndex::BirthdayIndex(){
	birthday = 0;
	index = -1;
	dim = 1;
}

BirthdayIndex::BirthdayIndex(double _b, long _index, int _d){
	birthday = _b;
	index = _index;
	dim = _d;
}

BirthdayIndex::BirthdayIndex(const BirthdayIndex& b){
	birthday = b.birthday;
	index = b.index;
	dim = b.dim;
}

void BirthdayIndex::copyBirthdayIndex(BirthdayIndex v){
	birthday = v.birthday;
	index = v.index;
	dim = v.dim;
}

double BirthdayIndex::getBirthday(){
	return birthday;
}

long BirthdayIndex::getIndex(){
	return index;
}

int BirthdayIndex::getDimension(){
	return dim;
}

void BirthdayIndex::print(){
	std::cout << "(dob:" << birthday << "," << index << ")" << std::endl;
}
	
bool BirthdayIndexComparator::operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const{
	if(o1.birthday == o2.birthday){
		if(o1.index < o2.index){
			return true;
		} else {
			return false;
		}
	} else {
		if(o1.birthday > o2.birthday){
			return true;
		} else {
			return false;
		}
	}
}

 	
bool BirthdayIndexInverseComparator::operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const{
	if(o1.birthday == o2.birthday){
		if(o1.index < o2.index){
			return false;
		} else {
			return true;
		}
	} else {
		if(o1.birthday > o2.birthday){
			return false;
		} else {
			return true;
		}
	}
}
	