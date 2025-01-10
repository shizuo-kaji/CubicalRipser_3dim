# CubicalRipser: Persistent Homology for 2D Image, 3D Voxel Data, and 1D Scalar Time Series

Written by 
- Takeki Sudo and Kazushi Ahara, Meiji University
- Shizuo Kaji, Kyushu University.

---

## Overview

CubicalRipser is an extension of [Ripser](http://ripser.org) by Ulrich Bauer, tailored for the efficient computation of persistent homology of cubical complexes.

### Key Features:
- **High Performance**: Among the fastest tools for computing persistent homology of 2D and 3D cubical complexes.
- **Flexible Filtrations**: Supports both **V-construction** and **T-construction** for cubical complexes ([details](#v-and-t-constructions)).
- **Binary Coefficients**: Computations are performed over the field with two elements.
- **Cross-Platform**: Python module and standalone command-line executable available.

For description, refer to the paper:  
*[Cubical Ripser: Software for Computing Persistent Homology of Image and Volume Data](https://arxiv.org/abs/2005.12692)*  
by Shizuo Kaji, Takeki Sudo, and Kazushi Ahara.

---

## License

CubicalRipser is open-source software licensed under the GNU Lesser General Public License v3.0 or later.  
Refer to the [LICENSE](LICENSE) file for more details.

---

## Getting Started

### Try Online
- **Google Colab Demo**: [CubicalRipser in Action](https://colab.research.google.com/github/shizuo-kaji/CubicalRipser_3dim/blob/master/demo/cubicalripser.ipynb)  
- **Topological Data Analysis (TDA) Tutorial**: [Hands-On Guide](https://colab.research.google.com/github/shizuo-kaji/TutorialTopologicalDataAnalysis/blob/master/TopologicalDataAnalysisWithPython.ipynb)  
- **Applications in Deep Learning**:  
  - [Example 1: Homology-enhanced CNNs](https://github.com/shizuo-kaji/HomologyCNN)  
  - [Example 2: Pretraining CNNs without Data](https://github.com/shizuo-kaji/PretrainCNNwithNoData)  

### Installation

#### Using `pip` (Recommended)
Install the Python module directly:  
```bash
pip install -U cripser
```

If you encounter architecture compatibility issues, try:  
```bash
pip uninstall cripser
pip install --no-binary cripser cripser
```

#### Building from Source
Requires a C++11-compatible compiler (e.g., GCC, Clang, MSVC).

1. Build the command-line executable:  
   ```bash
   cd build
   cmake ..
   make
   ```
   The executable `cubicalripser` will be created.

2. Alternatively, without `cmake`:  
   ```bash
   cd src
   make all
   ```
   Modify the `Makefile` if needed.

3. Install the Python module:  
   ```bash
   pip install .
   ```

#### Windows Notes
- Use a 64-bit compiler to match Python's architecture, e.g.:  
  ```bash
  cmake .. -G"Visual Studio 15 2017 Win64"
  cmake --build . --target ALL_BUILD --config Release
  ```
- Fix potential `ssize_t` issues in `pybind11` by adding:  
  ```cpp
  typedef SSIZE_T ssize_t;
  ```
  in `pybind11/include/pybind11/numpy.h` after `#if defined(_MSC_VER)`.

---

## Usage

### Python Module

Cubical Ripser accepts 1D/2D/3D Numpy arrays.

```python
import cripser
import numpy as np

arr = np.load("input.npy").astype(np.float64)
pd = cripser.computePH(arr, maxdim=2)
```
**Result**: A NumPy array of shape `(n, 9)` where each row contains:  
`dim, birth, death, x1, y1, z1, x2, y2, z2`.

They indicate the dimension of the cycle, birth-time, death-time, location (x1,y1,z1) of the cell giving birth to the cycle, and location (x2,y2,z2) of the cell destroying the cycle.

- To use the **T-construction**:  
  ```python
  import tcripser
  pd = tcripser.computePH(arr, maxdim=2)
  ```

### Command-Line Executable
```bash
./cubicalripser --print --maxdim 2 --output out.csv demo/3dimsample.txt
```
**Result**: `out.csv` with rows formatted as:  
`dim, birth, death, x1, y1, z1, x2, y2, z2`.

Each line consists of nine numbers indicating
the dimension of the cycle, birth-time, death-time, the creator location (x,y,z), and the destroyer location (x,y,z). 

For **Numpy arrays**:  
```bash
./cubicalripser --output result.csv input.npy
```

---

## Input Formats

### Supported Formats (command-line version)
- **NUMPY (.npy)**: Native format for both Python and CLI.  
- **Perseus Text (.txt)**: [Specification](http://people.maths.ox.ac.uk/nanda/perseus/).  
- **CSV (.csv)**: Simplified input for 2D images.
- **DIPHA (.complex)**: [Specification](https://github.com/DIPHA/dipha#file-formats).

### Image to Array Conversion
A small utility is included that converts images in various formats into NUMPY arrays.
- Convert images to `.npy`:  
  ```bash
  python demo/img2npy.py input.jpg output.npy
  ```
  A series of image files such as JPEG and PNG files (as long as the Pillow library can handle them)
can also be made into a volume in a similar way:
  ```bash
    python demo/img2npy.py input*.jpg volume.npy 
  ```
    Note that here we rely on the shell's path expansion. If your shell does not support it, you can manually specify file names as in the following:
  ```bash
    python demo/img2npy.py input00.dcm input01.dcm input02.dcm volume.npy 
  ```

- Handle DICOM volumes:
Given a series of DICOM files named **input00.dcm**, **input01.dcm**, **input02.dcm**... under the directory **dicom**,
we can convert the DICOM files to a single 3D Numpy array **volume.npy**
that is compatible with Cubical Ripser by
  ```bash
  python demo/img2npy.py dicom/*.dcm output.npy
  ```
  Or, we can compute persistent homology directly by
  ```bash
    python demo/cr.py dicom  --sort -it dcm -o output.csv
  ```
  by reading .dcm files from the directry **dicom** in a sorted order.

### DIPHA file
The filename should end with ".complex".
Look at [DIPHA binary format](https://github.com/DIPHA/dipha#file-formats) for specification.

We can convert input and output files between Cubical Ripser and DIPHA.
- to convert an Numpy array **img.npy** into DIPHA's format **img.complex**
  ```bash
    python dipha2npy.py img.npy img.complex 
  ```
- the other way around
  ```bash
    python dipha2npy.py img.complex img.npy
  ```
- convert DIPHA's output **result.output** into an Numpy array **result.npy**
  ```bash
    python dipha2npy.py result.output result.npy 
  ```
### 1D time series
A scalar time-series can be considered as a 1D image,
so Cubical Ripser can compute its persistent homology.
Note that other software would be more efficient for this purpose.

An example of regressing the frequency of noisy sine curves
is demonstrated [here](https://github.com/shizuo-kaji/TutorialTopologicalDataAnalysis).

---

## V and T Constructions

- **V-Construction**: Pixels represent 0-cells (4-neighbor connectivity in 2D).  
- **T-Construction**: Pixels represent top-cells (8-neighbor connectivity in 2D).

Use the appropriate executable for your needs:  
- **V-construction**: `cubicalripser` (Python module: `cripser`).  
- **T-construction**: `tcubicalripser` (Python module: `tcripser`).

By the Alexander duality, the following two give essentially the same results:

    ./cubicalripser input.npy
    ./tcubicalripser --embedded input.npy

The difference is in the sign of the filtration and the permanent cycle.
Here, (--embedded) converts the input I to -I^\infty described in the paper below.

For more details, see 
*[Duality in Persistent Homology of Images](https://arxiv.org/abs/2005.04597)* by Adélie Garin et al.

---

## Creator and Destroyer cells
The creator of a cycle is the cell which gives birth to the cycle. 
For example, the voxel in a connected component with the lowest filtration value creates a 0-dimensional cycle,
and the voxel which connects two separate connected components destroys the component with a higher birth time.
The creator and the destroyer cells are not uniquely determined, but they provide useful information to localise the cycle.
Cubical Ripser adopts the following convention on the location of these cells:
when the lifetime of a cycle is finte,

    arr[x2,y2,z2] - arr[x1,y1,z1] = death - birth = lifetime

where arr is the image, (x1,y1,z1) is the location of the creator cell, and (x2,y2,z2) is the location of the destroyer cell.
Note that when computed with the (--embedded) option, the roles of creator and destroyer are switched:

    arr[x1,y1,z1] - arr[x2,y2,z2] = death - birth = lifetime


The authors thank Nicholas Byrne for suggesting the convention and providing a test code.


---

## Deep Learning Integration

- **Lifetime Enhanced Image**: Adds topological features as additional channels for CNNs.  
  ```bash
  ./cubicalripser --output result.npy input.npy
  python demo/stackPH.py result.npy -o lifetime_image.npy -i input.npy
  ```
  In **lifetime_image.npy**, persistent homology is encoded as the extra channels so that it can be used as input for CNNs.

  Please look at the example section of [our paper](https://arxiv.org/abs/2005.12692).

- **Persistent Histogram Image**:  
  Similarly, the *persistent histogram image* can be obtained by
  ```bash
  python demo/stackPH.py result.npy -o hist_image.npy -t hist -i input.npy
  ```

  For practical examples, see [HomologyCNN](https://github.com/shizuo-kaji/HomologyCNN).

---





## Other software for persistent homology of cubical complexes
We give a referece to various software for persistent homology of images.
The comments are based on our limited understanding and tests, and hence, could be wrong.

- [Cubicle](https://bitbucket.org/hubwag/cubicle/src/master/) by Hubert Wagner

It computes for the V-construction of the image.
Its parallelised algorithm offers faster computation on multi-core machines.
Also, it reads the input image in small chunks so that it requires much less memory footprint.

- [HomcCube](https://i-obayashi.info/software.html) By Ippei Obayashi.

It computes for the V-construction of the image.
It is integrated into Homcloud developed by the same author.

- [DIPHA](https://github.com/DIPHA/dipha) by Ulrich Bauer and Michael Kerber

It computes for the V-construction of the image.
It is parallelised with MPI so it works on a cluster. 
The software has been used in various projects.
The memory footprint is relatively large.

- [GUDHI](http://gudhi.gforge.inria.fr/) developed at INRIA

It computes for the V- and T-construction of an array of any dimension.
It is well-documented and offers a well-organised and easy to use interface.
It focuses more on usability than performance.

- [diamorse](https://github.com/AppliedMathematicsANU/diamorse) developed at The Australian National University.

It computes for the V-construction of the image.

- [Perseus](http://people.maths.ox.ac.uk/nanda/perseus/) by Vidit Nanda

It computes for the V-construction of the image.

## Release Notes
- (v0.0.8) fixed memory leak in Python bindings (pointed out by Nicholas Byrne)
- (v0.0.7) slight speed up
- (v0.0.6) changes in the [definition of birth/death location](#Creator-and-Destroyer-cells) (suggested by Nicholas Byrne)
- (up to v0.0.5, difference from the [original version](https://github.com/CubicalRipser/CubicalRipser_3dim)
    - optimised codes (much less memory footprint, much faster for certain data; sometimes more than 100 times.)
    - Python friendly: see the Jupyter Notebook example found under the demo directory.
    - virtually infinite input size (compared to 510x510x510)
    - cache control
    - option to use the Alexander duality for the highest degree persistent homology
    - V and T construction for building cubical complexes from an image
    - output birth/death location
