/* simplex_coboundary_enumerator.h

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


class SimplexCoboundaryEnumerator
{
private:
	int count;
	DenseCubicalGrids* dcg;
public:
	BirthdayIndex simplex;
	BirthdayIndex nextCoface;

	SimplexCoboundaryEnumerator(DenseCubicalGrids* _dcg);

	void setSimplexCoboundaryEnumerator(BirthdayIndex& _s); 

	bool hasNextCoface(); 
};