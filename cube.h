/* birthday_index.h

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


class Cube
{
	
public:
	double birthday;
	long index;

	Cube();
		
	Cube(double _b, long _index);

	Cube(const Cube& b);

	void copyCube(const Cube& v);

	void print();
};

struct CubeComparator
{
	bool operator()(const Cube& o1, const Cube& o2) const;
};

struct CubeInverseComparator
{
	bool operator()(const Cube& o1, const Cube& o2) const;
};
