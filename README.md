# CubicalRipser : Persistent homology for 2D image and 3D voxel data (and 1D scalar timeseries)

copyright by Takeki Sudo and Kazushi Ahara, Meiji University, 2018

modified by Shizuo Kaji, Kyushu University, 2019

## Description
CubicalRipser is an adaptation of [Ripser](http://ripser.org) by Ulrich Bauer to computation of persistent homology of weighted cubical complexes.

- For 2 and 3 dimensional cubical complexes, CubicalRipser is among the fastest programs for computing persistent homology 
- Cubical Ripser implements both the V- and the T- constructions for the filtration of cubical complexes (see [V and T constructions](#V-and-T-constructions)).
- The coefficients are taken in the field with two elements.
- See [Other software for persistent homology of cubical complexes](#Other-software-for-persistent-homology-of-cubical-complexes)

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

    % pip install -U cripser

### Build from source
The command-line executable should be easily build with C++11 compilers such as G++, Clang, or Microsoft C++.
To build the command-line executable from source:

    % cd build
    % cmake ..
    % make

The executable is "cubicalripser".

If cmake is not available on your system, you can also do

    % cd src
    % make all

but perhaps you have to manually modify "Makefile".

To install Python module,

    % pip install .


### Windows specifics
On Windows, you may get an error like "Python config failure: Python is 64-bit, chosen compiler is 32-bit".
Then, you have to specify your compiler; for example

    % cmake .. -G"Visual Studio 15 2017 Win64"
    % cmake --build . --target ALL_BUILD --config Release

Also, due to the non-standard type used in pybind11, you may encounter an error 
saying "the type ssize_t is undefined". This error may be resolved by adding

    typedef SSIZE_T ssize_t;

right after the first appearance of 

    #if defined(_MSC_VER)

in pybind11/include/pybind11/numpy.h


## How to use
### Python module
To use from python,

    import cripser
    pd = cripser.computePH(arr,maxdim=2)

where arr is a 2D or 3D numpy array of type numpy.float64.
The result is stored in (n,9)-array, where n is the number of cycles.
Each row consists of

    dim birth   death   x1  y1  z1  x2  y2  z2

where (x1,y1,z1) is the location of the creator cell of the cycle and (x2,y2,z2) is the location of the destroyer cell of the cycle.
See also [Creator and Destroyer cells](#Creator-and-Destroyer-cells).

If you want to compute with the T-construction instead of the V-construction,

    import tcripser
    pd = tcripser.computePH(arr,maxdim=2)
    
Look at the Jupyter notebook demo/cubicalripser.ipynb and https://github.com/shizuo-kaji/HomologyCNN for practical usage.

### Command-line executable
(See also [2D Image file](#2D-Image-file) for a Python-based command-line utility.)

To see the command-line options:

    % ./cubicalripser

Example:

    % ./cubicalripser --print --maxdim 2 --output out.csv demo/3dimsample.txt

The results are recorded in **result.csv**.
Each line in the output **result.csv** consists of nine numbers indicating
the dimension of the cycle, birth-time, death-time, the creator location (x,y,z), and the destroyer location (x,y,z). 

Cubical Ripser accepts 1D/2D/3D Numpy arrays

    % ./cubicalripser --output result.csv input.npy

## Input file format
The python version accepts NUMPY arrays as input.
A small utility is included that converts images in various formats into NUMPU arrays.

The command-line version of CubicalRipser accepts three types of input files: NUMPY (.npy), Perseus TEXT (.txt), CSV (.csv), DIPHA (.complex).

### 2D Image file
Given a JPEG image **input.jpg**, we can compute its persistent homology by

    % python demo/cr.py input.jpg -o output.csv

The result is saved in the CSV file **output.csv** whose rows look

    dim birth   death   x1  y1  z1  x2  y2  z2

where (x1,y1,z1) is the location of the creator cell of the cycle and (x2,y2,z2) is the location of the destroyer cell of the cycle.

Alternatively, we can first convert the image into a 2D Numpy array **input.npy** by

    % python demo/img2npy.py input.jpg input.npy

and compute its persistent homology by the python module:
```
import numpy as np                                      # import the Numpy module
import cripser                                          # import the Cubical Ripser python module
arr = np.load("input.npy").astype(np.float64)           # load the image in the numpy array format
result = cripser.computePH(arr,maxdim=1)   # compute the persistent homology up to degree 1
```
Here, **result** is another 2D Numpy array of shape (M,9), where M is the number of cycles.
The none numbers of each row indicate the dimension of the cycle, birth-time, death-time, location (x1,y1,z1) of the cell giving birth to the cycle, and location (x2,y2,z2) of the cell destroying the cycle.

### 3D Volume file
Given a series of DICOM files named **input00.dcm**, **input01.dcm**, **input02.dcm**... under the directory **dicom**,
we can compute its persistent homology by

    % python demo/cr.py dicom  --sort -it dcm -o output.csv

by reading .dcm files from the directry **dicom** in a sorted order.

Alternatively, we can first convert the DICOM files to a single 3D Numpy array **volume.npy**
that is compatible with Cubical Ripser by

    % python demo/img2npy.py dicom/input*.dcm volume.npy 

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

    % ./cubicalripser --output result.npy input.npy
    % python demo/stackPH.py result.npy -o lifetime_image.npy -i input.npy

In **lifetime_image.npy**, persistent homology is encoded as the extra channels so that it can be used as input for CNNs.

Please look at the example section of [our paper](https://arxiv.org/abs/2005.12692)
and the [demonstration](https://github.com/shizuo-kaji/HomologyCNN) for details.

Similarly, the *persistent histogram image* can be obtained by

    % python demo/stackPH.py result.npy -o lifetime_image.npy -t hist -i input.npy


### CSV file (only for 2D image)
The filename should end with ".csv".

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

## V and T constructions
There are two major ways to build a filtred cubical complex from an image (that is, a function over a grid).
- In the V-construction, each pixel in the image corresponds to the 0 cell. 
- In the T-construction, each pixel in the image corresponds to the top cell.

In the 2D setting, the V-construction amounts to considering 4-neighbour pixel connectivity,
whereas the T-construction amounts to considering 8-neighbour pixel connectivity.

Cubical Ripser provides two versions of executables: 
- for the V-construction: cubicalripser, cripser (python module)
- for the T-construction: tcubicalripser (no python module provided)

By the Alexander duality, the following two give essentially the same results:

    ./cubicalripser input.npy
    ./tcubicalripser --embedded input.npy

The difference is in the sign of the filtration and the permanent cycle.
Here, (--embedded) converts the input I to -I^\infty described in the paper below.

Look at the following paper for details:
[Duality in Persistent Homology of Images by
Adélie Garin, Teresa Heiss, Kelly Maggs, Bea Bleile, Vanessa Robins](https://arxiv.org/abs/2005.04597)

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

It computes for the T-construction of the image.
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
