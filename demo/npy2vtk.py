#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import argparse
import os
from pyevtk.hl import *

#%%
parser = argparse.ArgumentParser("Convert numpy array to VTK")
parser.add_argument('from_fn', type=str, help="input npy")
args = parser.parse_args()

# %%
im = np.load(args.from_fn).astype(np.float64)
fn,ext = os.path.splitext(args.from_fn)
if(len(im.shape))==2:
    im = im[:,:,np.newaxis]
#w,h,c = im.shape
#x = np.arange(0, w+1)
#y = np.arange(0, h+1)
#z = np.arange(0, c+1)
#gridToVTK(args.to_fn, x, y, z, cellData = {'image': im})
imageToVTK(fn,pointData={'quantity':im})
