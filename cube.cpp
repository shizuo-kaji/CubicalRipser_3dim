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
	birth = 0;
	index = NONE;
}

Cube::Cube(double _b, unsigned short _x, unsigned short _y, unsigned short _z, unsigned short _m){
	birth = _b;
	index = _x | (_y<<10) | (_z<<20) | (_m<<30);
}

Cube::Cube(const Cube& v){
	birth = v.birth;
	index = v.index;
}

void Cube::copyCube(const Cube& v){
	birth = v.birth;
	index = v.index;
}

unsigned short Cube::x(){
	return( (index) & 1023);
}

unsigned short Cube::y(){
	return( (index >> 10) & 1023);
}

unsigned short Cube::z(){
	return( (index >> 20) & 1023);
}

unsigned short Cube::m(){
	return( (index >> 30) & 3);
}

void Cube::print(){
	std::cout << birth << "," << index << "," << x() << "," << y() << "," << z() << "," << m() << std::endl;
}

bool Cube::operator==(const Cube& rhs) const{
    return(index == rhs.index);
}

// true when b1>b2 (tie break i1<i2)
bool CubeComparator::operator()(const Cube& o1, const Cube& o2) const{
	if(o1.birth == o2.birth){
		return(o1.index < o2.index);
	} else {
		return(o1.birth > o2.birth);
	}
}
