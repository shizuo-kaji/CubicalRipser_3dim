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


class BirthdayIndex
{
	
public:
	double birthday;
	long index;
	int dim;

	BirthdayIndex();
		
	BirthdayIndex(double _b, long _index, int _d);

	BirthdayIndex(const BirthdayIndex& b);

	void copyBirthdayIndex(const BirthdayIndex& v);

	void print();
};

struct BirthdayIndexComparator
{
	bool operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const; 
};

struct BirthdayIndexInverseComparator
{
	bool operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const;	
};