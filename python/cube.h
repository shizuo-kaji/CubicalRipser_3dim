/* cube.h

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#define NONE 0xffffffffffffffff
#include <stdint.h>

class Cube
{
public:
	double birth;
	uint64_t index;

	Cube();
    
    // change the followings for bigger volume size
    Cube(double _b, uint32_t _x, uint32_t _y, uint32_t _z, uint8_t _m){
        birth = _b;
        index = (uint64_t)_x | ((uint64_t)_y<<20) | ((uint64_t)_z<<40) | ((uint64_t)_m<<60);
    };
    uint32_t x(){return( (index) & 0xfffff);}
    uint32_t y(){return( (index >> 20) & 0xfffff);}
    uint32_t z(){return( (index >> 40) & 0xfffff);}
    uint8_t m(){return( (index >> 60) & 0xf);}
    
    Cube(const Cube&);
	void copyCube(const Cube&);
	void print();

    
	bool operator==(const Cube& rhs) const;
};

struct CubeComparator{
	bool operator()(const Cube& o1, const Cube& o2) const;
};
