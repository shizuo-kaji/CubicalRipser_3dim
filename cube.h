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

#define NONE 0xffff

class Cube
{
public:
	double birth;
	uint32_t index;

	Cube();
    
    // change the followings for bigger volume size
    Cube(double _b, unsigned short _x, unsigned short _y, unsigned short _z, unsigned short _m){
        birth = _b;
        index = _x | (_y<<10) | (_z<<20) | (_m<<30);
    };
    unsigned short x(){return( (index) & 1023);}
    unsigned short y(){return( (index >> 10) & 1023);}
    unsigned short z(){return( (index >> 20) & 1023);}
    unsigned short m(){return( (index >> 30) & 3);}

    
    Cube(const Cube&);
	void copyCube(const Cube&);
	void print();

    
	bool operator==(const Cube& rhs) const;
};

struct CubeComparator{
	bool operator()(const Cube& o1, const Cube& o2) const;
};
