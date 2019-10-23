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

using namespace std;

class JointPairs{

	DenseCubicalGrids* dcg;
	vector<WritePairs> *wp;
	bool print;

public:
	JointPairs(DenseCubicalGrids* _dcg, vector<BirthdayIndex>& ctr, vector<WritePairs> &_wp, const bool _print);
	void joint_pairs_main( vector<BirthdayIndex>& ctr );
};