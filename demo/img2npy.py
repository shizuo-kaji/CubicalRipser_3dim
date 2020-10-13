#!/usr/bin/env python
# -*- coding: utf-8 -*-
import struct
import numpy as np
import argparse
import os
from PIL import Image
import re,struct
import shutil
num = lambda val : int(re.sub("\\D", "", val))
dtype ={"uint8": np.uint8, "uint16": np.uint16, "float": np.float32, "double": np.float64}

#%%
parser = argparse.ArgumentParser("Convert image file to Numpy array")
parser.add_argument('from_fn', nargs="*", help="multiple files of the same dimension would be stacked to a 3D image. If directory is specified, combine all the files in the directory.")
parser.add_argument('to_fn', help="output filename (.npy)")
parser.add_argument('--reduce','-r', type=int, default=1)
parser.add_argument('--tile','-t', type=int, default=1)
#parser.add_argument('--dtype','-d', type=str, default="double", choices=dtype.keys())
parser.add_argument('--sort','-s', action='store_true', help="Sort file names before stacking")
parser.add_argument('--zsplit','-z', action='store_true', help="save one file for each slice along z-axis")
args = parser.parse_args()

# %%
if os.path.isdir(args.from_fn[0]):
    fns = os.listdir(args.from_fn[0])
    args.from_fn = [os.path.join(args.from_fn[0],f) for f in fns]

if args.sort:
    args.from_fn.sort(key=num)

frfn,ext = os.path.splitext(args.from_fn[0])
if ext == ".dcm":
    try:
        import pydicom as dicom
    except:
        print("Install pydicom first by: pip install pydicom")
        exit()
        
images = []
for ffn in args.from_fn:
    print("processing {}".format(ffn))
    fn,ext = os.path.splitext(ffn)
    if ext == ".npy":
        im = np.load(ffn).astype(np.float64)
    elif ext == ".npz":
        imz = np.load(ffn)
        im = imz[imz.files[0]]
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

# size reduction
s = args.reduce
if len(img_arr.shape)==4:
    img_arr = img_arr[::s,::s,::s,::s]
    img_arr = np.tile(img_arr, (args.tile,args.tile,args.tile,args.tile))
elif len(img_arr.shape)==3:
    img_arr = img_arr[::s,::s,::s]
    img_arr = np.tile(img_arr, (args.tile,args.tile,args.tile))
elif len(img_arr.shape)==2:
    img_arr = img_arr[::s,::s]
    img_arr = np.tile(img_arr, (args.tile,args.tile))

# save
print("output ",args.to_fn, " shape: ",img_arr.shape)

tofn,ext = os.path.splitext(args.to_fn)
if ext == ".raw":  # for use with cubicle
    data = img_arr.astype(np.uint16).flatten()
    with open(args.to_fn, 'wb') as out:
        for v in data:
            out.write(struct.pack('H', v))   # B: uchar, H: ushort
elif ext == ".pgm":   # for use with diamorse
    if len(img_arr.shape)==3:
        with open(args.to_fn, 'wb') as outfile:
            for z in range(img_arr.shape[2]):
                fname = tofn+"_{:0>4}.pgm".format(z)
                Image.fromarray(img_arr[:,:,z].astype(np.uint8),mode='L').save(fname)
                with open(fname, 'rb') as infile:
                    shutil.copyfileobj(infile,outfile)
                os.remove(fname)
    else:
        Image.fromarray(img_arr.astype(np.uint8),mode='L').save(args.to_fn)
else:
    if args.zsplit:
        os.makedirs(args.to_fn,exist_ok=True)
        for i in range(img_arr.shape[2]):
            np.save(os.path.join(args.to_fn,frfn)+"_{:0>4}.npy".format(i),img_arr[:,:,i])
    else:
        np.save(args.to_fn,img_arr)

