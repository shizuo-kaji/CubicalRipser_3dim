######
import struct
import numpy as np
import argparse
import os
from PIL import Image

#%%
parser = argparse.ArgumentParser("Convert image file to Numpy")
parser.add_argument('from_fn')
parser.add_argument('to_fn')
parser.add_argument('--reduce','-r', type=int, default=1)
args = parser.parse_args()

s = args.reduce

# %%
if ".npy" in args.from_fn:
    im = np.load(args.from_fn).astype(np.float64)
else:
    im = np.array(Image.open(args.from_fn).convert('L'),dtype=np.float64)

print("input shape: ",im.shape)

if len(im.shape)==3:
    im = im[::s,::s,::s]
else:
    im = im[::s,::s]

print("output shape: ",im.shape)
np.save(args.to_fn,im)
