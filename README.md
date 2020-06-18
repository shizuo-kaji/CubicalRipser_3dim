# CubicalRipser : Persistent homology for 2D image and 3D voxel data (and 1D scalar timeseries)

copyright by Takeki Sudo and Kazushi Ahara, Meiji University, 2018

modified by Shizuo Kaji, Kyushu University, 2019

## Description
CubicalRipser is an adaptation of [Ripser](http://ripser.org) by Ulrich Bauer to computation of persistent homology of weighted cubical complexes.

For 2 and 3 dimensional cubical complexes, we believe CubicalRipser is currently the fastest program for computing persistent homology.

For details, please look at our paper
[Cubical Ripser: Software for computing persistent homology of image and volume data](https://arxiv.org/abs/2005.12692)
by Shizuo Kaji, Takeki Sudo, Kazushi Ahara.

## License
CubicalRipser is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

## Get Started
- You can try Cubical Ripser on [Google Colaboratory](https://colab.research.google.com/github/shizuo-kaji/CubicalRipser_3dim/blob/master/demo/cubicalripser.ipynb)
- You may also want to look at [A guide through TDA tools](https://colab.research.google.com/github/shizuo-kaji/TutorialTopologicalDataAnalysis/blob/master/TopologicalDataAnalysisWithPython.ipynb) giving a hands-on tutorial for various tools for topological data analysis including Cubical Ripser
- How Deep-learning and Persistent homology can be combined is demonstrated at https://github.com/shizuo-kaji/HomologyCNN

## Installation
### Recommended: pip
Install the Python module only:

    % pip install git+https://github.com/shizuo-kaji/CubicalRipser_3dim

### Binary
Precompiled Windows/Mac command-line binaries are found under win64/mac directories.
("cripser.cpython-37m-darwin.so" or "cripser.cp37-win_amd64.pyd" are python modules.
You can simply copy them to the same directory as the python script, if you have the right version of Python.)

### Build from source
The command-line executable should be easily build with C++11 compilers such as G++, Clang, or Microsoft C++.
To build the command-line executable from source:

    % cd build
    % cmake ..
    % make

On Windows, you may get an error like "Python config failure: Python is 64-bit, chosen compiler is 32-bit".
Then, you have to specify your compiler; for example

    % cmake .. -G"Visual Studio 15 2017 Win64"


The executable is "cubicalripser".

If cmake is not available on your system, you can also do

    % cd src
    % make all

but perhaps you have to manually modify "Makefile".

To install Python module,

    % pip install .


## How to use
### Python module
To use from python,

    import cripser
    cripser.computePH(arr,maxdim=2,location="birth")

where arr is a 2D or 3D numpy array of type numpy.float64.

Look at the Jupyter notebook demo/cubicalripser.ipynb and https://github.com/shizuo-kaji/HomologyCNN for practical usage.

### Command-line executable
To see the command-line options:

    % ./cubicalripser

Example:

    % ./cubicalripser --print --location birth --maxdim 2 --output out.csv demo/3dimsample.txt


## Input file format
CubicalRipser accepts three types of input files: NUMPY, TEXT, DIPHA.

### 2D Image file
Given a JPEG image **input.jpg**, we can convert it into a 2D Numpy array **input.npy** by

    % python demo/img2npy.py input.jpg input.npy

Then, we can compute its persistent homology by the python module:
```
import numpy as np                                      # import the Numpy module
import cripser                                          # import the Cubical Ripser python module
arr = np.load("input.npy").astype(np.float64)           # load the image in the numpy array format
result = cripser.computePH(arr,maxdim=1,location="birth")   # compute the persistent homology up to degree 1
```
Here, **result** is another 2D Numpy array of shape (M,6), where M is the number of cycles.
The six numbers of each row indicate the dimension of the cycle, birth-time, death-time, and location (x,y,z) of the cell giving birth to the cycle.

Alternatively, one can use the command-line executable to compute the persistent homology
of the 1D/2D/3D Numpy array **input.npy** and obtain results in **result.csv**.

    % ./cubicalripser --location birth --output result.csv input.npy

Each line in the output **result.csv** consists of six numbers indicating
the dimension of the cycle, birth-time, death-time, and location (x,y,z). 

### 3D Volume file
Given a series of DICOM files named **input00.dcm**, **input01.dcm**, **input02.dcm**...,
one can convert them to a single 3D Numpy array **volume.npy**
that is compatible with Cubical Ripser by

    % python demo/img2npy.py input*.dcm volume.npy 

A series of image files such as JPEG and PNG files (as long as the Pillow library can handle them)
can also be made into a volume in a similar way:

    % python demo/img2npy.py input*.jpg volume.npy 

Note that here we rely on the shell's path expansion.
If your shell does not support it,
you can manually specify file names as in the following:

    % python demo/img2npy.py input00.dcm input01.dcm input02.dcm volume.npy 

### 1D time series
A scalar time-series can be considered as a 1D image,
so Cubical Ripser can compute its persistent homology.
Note that other software would be more efficient for this purpose.

An example of regressing the frequency of noisy sine curves
is demonstrated [here](https://github.com/shizuo-kaji/TutorialTopologicalDataAnalysis).


### Deep Learning X Persistent homology

*Lifetime enhanced image* is a way to feed the topological features obtained by persistent homology
into convolutional neural networks (CNNs).

    % ./cubicalripser --location birth --output result.npy input.npy
    % python demo/stackPH.py result.npy -o lifetime_image.npy -i input.npy

In **lifetime_image.npy**, persistent homology is encoded as the extra channels so that it can be used as input for CNNs.

Please look at the example section of [our paper](https://arxiv.org/abs/2005.12692)
and the [demonstration](https://github.com/shizuo-kaji/HomologyCNN) for details.

Similarly, the *persistent histogram image* can be obtained by

    % python demo/stackPH.py result.npy -o lifetime_image.npy -t hist -i input.npy


### Text file (Perseus)
The filename should end with ".txt".
Please look at [Perseus Dense Cubical Grids format](http://people.maths.ox.ac.uk/nanda/perseus/) for specification.

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

### DIPHA file
The filename should end with ".complex".
Look at [DIPHA binary format](https://github.com/DIPHA/dipha#file-formats) for specification.

We can convert input and output files between Cubical Ripser and DIPHA.
- to convert an Numpy array **img.npy** into DIPHA's format **img.complex**

    % python dipha2npy.py img.npy img.complex 

- the other way around

    % python dipha2npy.py img.complex img.npy

- convert DIPHA's output **result.output** into an Numpy array **result.npy**

    % python dipha2npy.py result.output result.npy 


## Difference with the original version
- cleaned up/optimised codes (much less memory footprint, much faster for certain data; sometimes more than 100 times.)
- Python friendly: see the Jupyter Notebook example found under the demo directory.
- virtually infinite input size (compared to 510x510x510)
- cache control
- option to use the Alexander duality for the highest degree persistent homology
- output birth/death location

## TODO
- Do no use bounding cells so that we can reduce memory usage. (but using an accessor for cell weights affects performance, so we have to find another way.)
