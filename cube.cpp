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
#include "cube.h"

using namespace std;

Cube::Cube(){
	birthday = 0;
	index = -1;
}

Cube::Cube(double _b, long _index){
	birthday = _b;
	index = _index;
}

Cube::Cube(const Cube& b){
	birthday = b.birthday;
	index = b.index;
}

void Cube::copyCube(const Cube& v){
	birthday = v.birthday;
	index = v.index;
}

void Cube::print(){
	std::cout << "(dob:" << birthday << "," << index << ")" << std::endl;
}

// true when b1>b2 (tie break i1<i2)
bool CubeComparator::operator()(const Cube& o1, const Cube& o2) const{
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

// true when b1<b2 (tie break i1<i2, same as above)
bool CubeInverseComparator::operator()(const Cube& o1, const Cube& o2) const{
	if(o1.birthday == o2.birthday){
		if(o1.index > o2.index){
			return true;
		} else {
			return false;
		}
	} else {
		if(o1.birthday > o2.birthday){
			return false;
		} else {
			return true;
		}
	}
}
