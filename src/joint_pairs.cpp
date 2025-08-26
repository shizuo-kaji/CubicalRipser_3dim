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
        for (uint32_t w = 0; w < dcg->aw; ++w) {
            for (uint32_t z = 0; z < dcg->az; ++z) {
                for (uint32_t y = 0; y < dcg->ay; ++y) {
                    for (uint32_t x = 0; x < dcg->ax; ++x) {
                        double birth = dcg->getBirth(x, y, z, w, m, 1);
                        // If birth value is below the threshold, add to the list
                        if (birth < config->threshold) {
                            ctr.emplace_back(birth, x, y, z, w, m);
                        }
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

    auto decode = [&](uint64_t idx, uint32_t& x, uint32_t& y, uint32_t& z, uint32_t& w) {
        uint64_t t = idx;
        x = t % dcg->ax; t /= dcg->ax;
        y = t % dcg->ay; t /= dcg->ay;
        z = t % dcg->az;
        w = t / dcg->az;
    };
    // Process cubes in reverse order (starting from the highest birth time)
    for (auto e = ctr.rbegin(); e != ctr.rend(); ++e) {
        // Calculate the linear index for the union-find structure
        uint64_t uind, vind;

        if (dcg->dim == 4) {
            // 4D indexing
            uind = e->x() + dcg->ax * e->y() + dcg->axy * e->z() + dcg->axyz * e->w();

            // 4D neighbor offsets for edge types
            static const int8_t dx4d[4] = {1, 0, 0, 0};  // x, y, z, w edges
            static const int8_t dy4d[4] = {0, 1, 0, 0};
            static const int8_t dz4d[4] = {0, 0, 1, 0};
            static const int8_t dw4d[4] = {0, 0, 0, 1};

            int m = e->m();
            if (m < 0 || m >= 4) std::exit(-1);

            vind = static_cast<uint64_t>((int64_t)uind +
                          dx4d[m] +
                          (int64_t)dcg->ax * dy4d[m] +
                          (int64_t)dcg->axy * dz4d[m] +
                          (int64_t)dcg->axyz * dw4d[m]);
        } else {
            // 3D indexing
            uind = e->x() + dcg->ax * e->y() + dcg->axy * e->z();

            static const int8_t dx[13]={1,0,0,1,1 ,0, 0,1, 1,1,1, 1, 1};
            static const int8_t dy[13]={0,1,0,1,-1,-1,1,-1,0,1,-1,0, 1};
            static const int8_t dz[13]={0,0,1,0,0, 1, 1,1, 1,1,-1,-1,-1};
            int m = e->m();
            if (m < 0 || m >= 13) std::exit(-1);
            vind = static_cast<uint64_t>( (int64_t)uind +
                          dx[m] +
                          (int64_t)dcg->ax * dy[m] +
                          (int64_t)dcg->axy * dz[m] );
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
                uint32_t bx, by, bz, bw, dx, dy, dz, dw;
                decode(birth_ind, bx, by, bz, bw);
                decode(death_ind, dx, dy, dz, dw);

                if (config->tconstruction) {
                    wp->emplace_back(current_dim,
                        Cube(birth, bx, by, bz, bw, 0),
                        Cube(death, dx, dy, dz, dw, 0),
                        dcg, config->print);
                } else {
                    wp->emplace_back(current_dim, birth, death,
                        bx, by, bz, bw, dx, dy, dz, dw, config->print);
                }
            }
            e->index = NONE;  // Mark edge as processed
        }
    }

    // Handle the base point component for H_0
    if (current_dim == 0) {
        uint32_t bx, by, bz, bw, dx, dy, dz, dw;
        decode(min_idx, bx, by, bz, bw);
        wp->emplace_back(current_dim, min_birth, dcg->threshold, bx, by, bz, bw, 0, 0, 0, 0, config->print);
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
