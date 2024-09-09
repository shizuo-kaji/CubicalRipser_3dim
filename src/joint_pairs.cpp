/* joint_pairs.cpp

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdint>

#include "cube.h"
#include "dense_cubical_grids.h"
#include "coboundary_enumerator.h"
#include "union_find.h"
#include "write_pairs.h"
#include "joint_pairs.h"

using namespace std;


// Constructor for JointPairs
JointPairs::JointPairs(DenseCubicalGrids* _dcg, vector<WritePairs>& _wp, Config& _config)
    : dcg(_dcg), wp(&_wp), config(&_config) {}

// Enumerate all edges based on given types
void JointPairs::enum_edges(const vector<uint8_t>& types, vector<Cube>& ctr) {
    ctr.clear();
    // Iterate over each type (order of loops matters for performance)
    for (const auto& m : types) {
        for (uint32_t z = 0; z < dcg->az; ++z) {
            for (uint32_t y = 0; y < dcg->ay; ++y) {
                for (uint32_t x = 0; x < dcg->ax; ++x) {
                    double birth = dcg->getBirth(x, y, z, m, 1);
                    // If birth value is below the threshold, add to the list
                    if (birth < config->threshold) {
                        ctr.emplace_back(birth, x, y, z, m);
                    }
                }
            }
        }
    }
    // Sort the cubes based on birth values
    std::sort(ctr.begin(), ctr.end(), CubeComparator());
}

// Compute H_0 by union-find
void JointPairs::joint_pairs_main(vector<Cube>& ctr, int current_dim) {
    UnionFind dset(dcg);
    uint64_t u, v = 0;
    double min_birth = config->threshold;
    uint64_t min_idx = 0;

    // Process cubes in reverse order (starting from the highest birth time)
    for (auto e = ctr.rbegin(); e != ctr.rend(); ++e) {
        // Calculate the linear index for the union-find structure
        uint64_t uind = e->x() + dcg->ax * e->y() + dcg->axy * e->z();
        uint64_t vind;

        // Determine the corresponding neighbor based on the cube's type
        switch (e->m()) {
            case 0: vind = uind + 1; break;                        // x+1
            case 1: vind = uind + dcg->ax; break;                  // y+1
            case 2: vind = uind + dcg->axy; break;                 // z+1
            case 3: vind = uind + 1 + dcg->ax; break;              // x+1, y+1
            case 4: vind = uind + 1 - dcg->ax; break;              // x+1, y-1
            case 5: vind = uind - dcg->ax + dcg->axy; break;       // y-1, z+1
            case 6: vind = uind + dcg->ax + dcg->axy; break;       // y+1, z+1
            case 7: vind = uind + 1 - dcg->ax + dcg->axy; break;   // x+1, y-1, z+1
            case 8: vind = uind + 1 + dcg->axy; break;             // x+1, z+1
            case 9: vind = uind + 1 + dcg->ax + dcg->axy; break;   // x+1, y+1, z+1
            case 10: vind = uind + 1 - dcg->ax - dcg->axy; break;  // x+1, y-1, z-1
            case 11: vind = uind + 1 - dcg->axy; break;            // x+1, z-1
            case 12: vind = uind + 1 + dcg->ax - dcg->axy; break;  // x+1, y+1, z-1
            default: exit(-1);  // Invalid type, exit with error
        }

        u = dset.find(uind);
        v = dset.find(vind);

        if (u != v) {  // If u and v are not already connected
            double birth;
            uint64_t birth_ind, death_ind;

            // Determine which component is younger and will be merged
            if (dset.birthtime[u] >= dset.birthtime[v]) {
                birth = dset.birthtime[u];
                birth_ind = current_dim == 0 ? u : (dset.birthtime[uind] > dset.birthtime[vind] ? uind : vind);
                death_ind = current_dim == 0 ? (dset.birthtime[uind] > dset.birthtime[vind] ? uind : vind) : u;
                if (dset.birthtime[v] < min_birth) {
                    min_birth = dset.birthtime[v];
                    min_idx = v;
                }
            } else {
                birth = dset.birthtime[v];
                birth_ind = current_dim == 0 ? v : (dset.birthtime[uind] > dset.birthtime[vind] ? uind : vind);
                death_ind = current_dim == 0 ? (dset.birthtime[uind] > dset.birthtime[vind] ? uind : vind) : v;
                if (dset.birthtime[u] < min_birth) {
                    min_birth = dset.birthtime[u];
                    min_idx = u;
                }
            }

            double death = e->birth;
            dset.link(u, v);  // Union the sets

            // Record the birth-death pair if they are not equal
            if (birth != death) {
                if (config->tconstruction) {
                    wp->emplace_back(current_dim, Cube(birth, birth_ind % dcg->ax, (birth_ind / dcg->ax) % dcg->ay, (birth_ind / dcg->axy) % dcg->az, 0),
                                     Cube(death, death_ind % dcg->ax, (death_ind / dcg->ax) % dcg->ay, (death_ind / dcg->axy) % dcg->az, 0), dcg, config->print);
                } else {
                    wp->emplace_back(current_dim, birth, death, birth_ind % dcg->ax, (birth_ind / dcg->ax) % dcg->ay, (birth_ind / dcg->axy) % dcg->az,
                                     death_ind % dcg->ax, (death_ind / dcg->ax) % dcg->ay, (death_ind / dcg->axy) % dcg->az, config->print);
                }
            }
            e->index = NONE;  // Mark edge as processed
        }
    }

    // Handle the base point component for H_0
    if (current_dim == 0) {
        wp->emplace_back(current_dim, min_birth, dcg->threshold, min_idx % dcg->ax, (min_idx / dcg->ax) % dcg->ay, (min_idx / dcg->axy) % dcg->az, 0, 0, 0, config->print);
    }

    // Remove unnecessary edges and optimize storage
    if (config->maxdim == 0 || current_dim > 0) {
        return;  // Skip further processing if we're not handling the highest dimension
    } else {
        auto new_end = std::remove_if(ctr.begin(), ctr.end(), [](const Cube& e) { return e.index == NONE; });
        ctr.erase(new_end, ctr.end());
        // No need to sort again since ctr was already sorted
		//	cout << ctr.size() << endl;
		//	std::sort(ctr.begin(), ctr.end(), CubeComparator()); // we can skip sorting as it is already sorted
	}
}