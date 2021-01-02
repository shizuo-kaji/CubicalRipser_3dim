/* cube.cpp

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

void Cube::copyCube(const Cube& v){
	birth = v.birth;
	index = v.index;
}

void Cube::print(){
	std::cout << birth << "," << index << "," << x() << "," << y() << "," << z() << "," << (uint32_t)m() << std::endl;
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
