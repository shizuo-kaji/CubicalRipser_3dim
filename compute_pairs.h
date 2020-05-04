/* compute_pairs.h

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
#include <unordered_map>

using namespace std;

typedef priority_queue<Cube, vector<Cube>, CubeComparator> CubeQue;

class ComputePairs{
private:
	DenseCubicalGrids * dcg;
	unordered_map<uint64_t, uint32_t> pivot_column_index;
	uint8_t dim;
	vector<WritePairs> *wp;
	bool print;

public:
	ComputePairs(DenseCubicalGrids* _dcg, vector<WritePairs> &_wp, const bool _print);
	void compute_pairs_main(vector<Cube>& ctr, int min_cache_size);
	void assemble_columns_to_reduce(vector<Cube>& ctr, uint8_t _dim);
	void add_cache(uint32_t i, CubeQue &wc, unordered_map<uint32_t, CubeQue>& recorded_wc);
	Cube pop_pivot(vector<Cube>& column);
	Cube get_pivot(vector<Cube>& column);
	Cube pop_pivot(CubeQue& column);
	Cube get_pivot(CubeQue& column);
};
