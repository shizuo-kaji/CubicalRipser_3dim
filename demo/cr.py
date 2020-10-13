#!/usr/bin/env python
# -*- coding: utf-8 -*-
# test for python module
#%%
import struct
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import cripser
from PIL import Image
import re
num = lambda val : int(re.sub("\\D", "", val))

#%%
parser = argparse.ArgumentParser("")
parser.add_argument('input',type=str, nargs="*", help="numpy array or multiple images or directory containing multiple images")
parser.add_argument('--output', '-o', default=None)
parser.add_argument('--top_dim',action='store_true')
parser.add_argument('--embedded', '-e', action='store_true')
parser.add_argument('--maxdim','-m', default=2,type=int)
parser.add_argument('--sort','-s', action='store_true', help="Sort file names before stacking")
parser.add_argument('--software',type=str,default="cubicalripser")
args = parser.parse_args()

if os.path.isdir(args.input[0]):
    fns = os.listdir(args.input[0])
    args.input = [os.path.join(args.input[0],f) for f in fns]

if args.sort:
    args.input.sort(key=num)

fn,ext = os.path.splitext(args.input[0])
if ext == ".dcm":
    try:
        import pydicom as dicom
    except:
        print("Install pydicom first by: pip install pydicom")
        exit()
        
images = []

for ffn in args.input:
    print("processing {}".format(ffn))
    fn,ext = os.path.splitext(ffn)
    if ext == ".npy":
        im = np.load(ffn).astype(np.float64)
    elif ext == ".dcm":
        ref_dicom_in = dicom.read_file(ffn, force=True)
#        ref_dicom_in.file_meta.TransferSyntaxUID = dicom.uid.ImplicitVRLittleEndian
        im = ref_dicom_in.pixel_array.astype(np.float64) +ref_dicom_in.RescaleIntercept
    elif ext == ".csv":
        im = np.loadtxt(ffn,delimiter=",")
    else:
        im = np.array(Image.open(ffn).convert('L'),dtype=np.float64)

    images.append(im)

img_arr = np.squeeze(np.stack(images,axis=-1))
print("input shape: ",img_arr.shape)

if args.software=="gudhi":
    try:
        import gudhi
    except:
        print("Install pydicom first by: conda install -c conda-forge gudhi")
        exit()
    
    print("Computing PH with GUDHI (T-construction)")
    gd = gudhi.CubicalComplex(top_dimensional_cells=img_arr)
#    gd.compute_persistence()
    res = np.array(gd.persistence(2,0)) # coeff = 2
    print("Betti numbers: ", gd.persistent_betti_numbers(np.inf,-np.inf))

else:
    res = cripser.computePH(img_arr,maxdim=args.maxdim,top_dim=args.top_dim,embedded=args.embedded)
    print("Betti numbers: ", [res[res[:,0]==i].shape[0] for i in range(3)])
#    print(res[:10])

if args.output is not None:
    np.save(args.output,res)
