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
	vector<WritePairs> *wp;
	Config* config;
public:
	DenseCubicalGrids* dcg;
	JointPairs(DenseCubicalGrids* _dcg, vector<WritePairs> &_wp, Config&);
	void enum_edges(std::vector<uint8_t>,vector<Cube>&);
	void enum_edges_alex(std::vector<uint8_t>,vector<Cube>&);
	void joint_pairs_main( vector<Cube>& ctr, int current_dim);
};
