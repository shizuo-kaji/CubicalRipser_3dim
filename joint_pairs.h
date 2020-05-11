/* joint_pairs.h

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
#include <cstdint>
#include "config.h"

using namespace std;

class JointPairs{
private:
	DenseCubicalGrids* dcg;
	vector<WritePairs> *wp;
	Config* config;
public:
	JointPairs(DenseCubicalGrids* _dcg, vector<Cube>& ctr, vector<WritePairs> &_wp, Config&);
	void joint_pairs_main( vector<Cube>& ctr );
};
