# CubicalRipser : Persistent homology for 2D image and 3D voxel data (and 1D scalar timeseries)

copyright by Takeki Sudo and Kazushi Ahara, Meiji University, 2018

modified by Shizuo Kaji, Kyushu University, 2019

## Description
CubicalRipser is an adaptation of [Ripser](http://ripser.org) by Ulrich Bauer for a filtered cubical complex.
For 2 and 3 dimensional cubical complexes, we believe CubicalRipser is currently the fastest program for computing persistent homology.

For details, please look at our preprint
[Cubical Ripser: Software for computing persistent homology of image and volume data](https://arxiv.org/abs/2005.12692)
by Shizuo Kaji, Takeki Sudo, Kazushi Ahara.

## License
CubicalRipser is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

## Get Started
- You can try Cubical Ripser on [Google Colaboratory](https://colab.research.google.com/github/shizuo-kaji/CubicalRipser_3dim/blob/master/demo/cubicalripser.ipynb)
- You may also want to look at [A guide through TDA tools](https://colab.research.google.com/github/shizuo-kaji/TutorialTopologicalDataAnalysis/blob/master/TopologicalDataAnalysisWithPython.ipynb) giving hands-on for various tools including Cubical Ripser
- How Deep-learning and Persistent homology can be combined is demonstrated at https://github.com/shizuo-kaji/HomologyCNN

## Installation
(Recommended)
Install the Python module only:

    % pip install git+https://github.com/shizuo-kaji/CubicalRipser_3dim

(For those who need commend-line versions as well)
Precompiled Windows/Mac command-line binaries are found under win64/mac directories.
("cripser.cpython-37m-darwin.so" or "cripser.cp37-win_amd64.pyd" are python modules.
You can simply copy them to the same directory as the python script, if you have the right version of Python.)

To build the command-line executable from source:

    % cd build
    % cmake ..
    % make

On Windows, you may get an error like "Python config failure: Python is 64-bit, chosen compiler is 32-bit".
Then, you have to specify your compiler; for example

    % cmake .. -G"Visual Studio 15 2017 Win64"


The executable is "cubicalripser".

If cmake is not available on your system,

    % cd src
    % make all

but perhaps you have to manually modify "Makefile".

To install Python module,

    % pip install .


## How to use
To see the command-line options:

    % ./cubicalripser

Example:

    % ./cubicalripser --print --location birth --maxdim 2 --output out.csv demo/3dimsample.txt

To use from python,

    import cripser
    cripser.computePH(arr,maxdim=2,location="birth")

where arr is a 2D or 3D numpy array of type numpy.float64.

Look at the Jupyter notebook demo/cubicalripser.ipynb and https://github.com/shizuo-kaji/HomologyCNN for practical usage.


## Input file format
CubicalRipser takes three types of input files: NUMPY, TEXT, DIPHA.
Image files (JPEG, PNG, etc.) can be converted to NUMPY array by

    % python demo/img2npy.py input.jpg output.npy

- 1D/2D/3D Numpy array (recommended). The filename should end with ".npy". DType must be float64. see the Jupyter Notebook example found under the demo directory.
- Text file. The filename should end with ".txt"
```
3
max_x
max_y
max_z
val[1,1,1]
val[2,1,1]
...
val[max_x,max_y,max_z]
```
- [DIPHA binary format](https://github.com/DIPHA/dipha#file-formats). The filename should end with ".complex".


## Difference with the original version
I added the following functionality:
- input/output numpy array file (.npy): see the Jupyter Notebook example found under the demo directory.
- input size up to 2^20-1 x 2^20-1 x 2^20-1
- cleaned up/optimised codes (much less memory footprint, much faster for certain data)
- cache control
- option to use the Alexander duality for the highest degree persistent homology
- output birth/death location
- python binding

## TODO
- When using Alexander duality for 2d image, it saves memory if we use the 2d complex (currently, we compute for a thickened 3d complex having the same PH)
