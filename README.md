# CubicalRipser_3dim : Persistent homology for 3D voxel data

copyright by Takeki Sudo and Kazushi Ahara, Meiji University, 2018

modified by Shizuo Kaji, Kyushu University, 2019

## Description
CubicalRipser depends heavily on [Ripser](http://ripser.org), software for calculating Vietoris-Rips 
persistent homology by Ulrich Bauer, 2015-2016.  
CubicalRipser is an adaptation for a filtered cubical complex.
In 2 and 3 dimensional case, we believe that CubicalRipser is much faster than DIPHA.

## License
CubicalRipser is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

## How to use
To build the software:

    % make

To see the command-line options:

    % ./CR3

Example:

    % ./CR3 --format perseus dat/3dimsample1.txt --print --location --maxdim 2 --output out.csv

