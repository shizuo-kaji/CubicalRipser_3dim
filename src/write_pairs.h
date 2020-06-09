/* write_pairs.h

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <cstdint>

class WritePairs
{
public:
	uint8_t dim;
	double birth;
	double death;
	uint32_t birth_x;
	uint32_t birth_y;
	uint32_t birth_z;

	WritePairs(uint8_t _dim, double _birth, double _death, uint32_t _birth_x,uint32_t _birth_y,uint32_t _birth_z, bool print = false){
        dim = _dim;
        birth = _birth;
        death = _death;
        birth_x = _birth_x;
        birth_y = _birth_y;
        birth_z = _birth_z;
        if (print == true) {
            std::cout << "[" << birth << "," << death << ")" << " location (" << birth_x << "," << birth_y << "," << birth_z << ")" << std::endl;
        }
    }

};
