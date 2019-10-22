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

template <class Key, class T> class hash_map : public std::unordered_map<Key, T> {};

class ComputePairs
{
public:
	DenseCubicalGrids* dcg;
	hash_map<int, int> pivot_column_index;
	int ax, ay, az;
	int dim;
	vector<WritePairs> *wp;
	bool print;

	ComputePairs(DenseCubicalGrids* _dcg, vector<WritePairs> &_wp, const bool _print);

	void compute_pairs_main(vector<BirthdayIndex>& ctr);

	BirthdayIndex pop_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
		column);

	BirthdayIndex get_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
		column);

	void assemble_columns_to_reduce(vector<BirthdayIndex>& ctr, int _dim);
};