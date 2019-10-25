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
	int dim;
	double birth;
	double death;
	int birth_x;
	int birth_y;
	int birth_z;

	WritePairs(int _dim, double _birth, double _death, int _birth_x,int _birth_y,int _birth_z, bool print = false);
	
};