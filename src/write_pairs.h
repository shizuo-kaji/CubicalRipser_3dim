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
#include "dense_cubical_grids.h"

class WritePairs
{
public:
	uint8_t dim;
	double birth;
	double death;
	uint32_t birth_x,birth_y,birth_z;
	uint32_t death_x,death_y,death_z;

	WritePairs(uint8_t _dim, double _birth, double _death, uint32_t _birth_x,uint32_t _birth_y,uint32_t _birth_z, uint32_t _death_x,uint32_t _death_y,uint32_t _death_z, bool print = false){
        dim = _dim;
        birth = _birth;
        death = _death;
        birth_x = _birth_x;
        birth_y = _birth_y;
        birth_z = _birth_z;
        death_x = _death_x;
        death_y = _death_y;
        death_z = _death_z;
        if (print == true) {
            std::cout << "[" << birth << "," << death << ")" << " birth loc. (" << birth_x << "," << birth_y << "," << birth_z << "), " << " death loc. (" << death_x << "," << death_y << "," << death_z << ")" << std::endl;
        }
    }
    WritePairs(uint8_t _dim, Cube _birthC, Cube _deathC, DenseCubicalGrids* _dcg, bool print = false){
        dim = _dim;
        birth = _birthC.birth;
        death = _deathC.birth;
        auto b =  _dcg->ParentVoxel(dim, _birthC);
        auto d =  _dcg->ParentVoxel(dim, _deathC);
        // vector<uint32_t> b = {0,0,0};
        // vector<uint32_t> d = {0,0,0};
        birth_x=b[0];
        birth_y=b[1];
        birth_z=b[2];
        death_x=d[0];
        death_y=d[1];
        death_z=d[2];
        if (print == true) {
            std::cout << "[" << birth << "," << death << ")" << " birth loc. (" << birth_x << "," << birth_y << "," << birth_z << "), " << " death loc. (" << death_x << "," << death_y << "," << death_z << ")" << std::endl;
        }
    }

};
