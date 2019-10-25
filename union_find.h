/* union_find.h

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <vector>

using namespace std;

class UnionFind{
private:
	vector<int> parent;
	vector<double> time_max;
public:
	vector<double> birthtime;
	UnionFind(DenseCubicalGrids* _dcg);
	long find(long x); 
	void link(long x, long y);
};
