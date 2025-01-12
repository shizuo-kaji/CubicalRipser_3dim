/*
This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "cubicalripser_pybind.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

using namespace std;

namespace py = pybind11;

PYBIND11_MODULE(tcripser, m) {
    m.doc() = R"pbdoc(
        Cubical Ripser (T-construction) plugin
        -----------------------
        .. currentmodule:: tcripser
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

    m.def("computePH", &computePH, R"pbdoc(Compute Persistent Homology
    )pbdoc", py::arg("arr"),  py::arg("maxdim")=2, py::arg("top_dim")=false, py::arg("embedded")=false, py::arg("location")="yes");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}