######
import struct
import numpy as np
import argparse
import os
from PIL import Image

#%%
parser = argparse.ArgumentParser("Convert image file to Numpy array")
parser.add_argument('from_fn', nargs="*", help="multiple files of the same dimension would be stacked to a 3D image")
parser.add_argument('to_fn', help="output filename (.npy)")
parser.add_argument('--reduce','-r', type=int, default=1)
parser.add_argument('--tile','-t', type=int, default=1)
args = parser.parse_args()

s = args.reduce

# %%
fn,ext = os.path.splitext(args.from_fn[0])
if ext == ".dcm":
    import pydicom as dicom

images = []
for ffn in args.from_fn:
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

# size reduction
if len(img_arr.shape)==4:
    img_arr = img_arr[::s,::s,::s,::s]
    img_arr = np.tile(img_arr, (args.tile,args.tile,args.tile,args.tile))
elif len(img_arr.shape)==3:
    img_arr = img_arr[::s,::s,::s]
    img_arr = np.tile(img_arr, (args.tile,args.tile,args.tile))
elif len(img_arr.shape)==2:
    img_arr = img_arr[::s,::s]
    img_arr = np.tile(img_arr, (args.tile,args.tile))

# save to npy
print("output shape: ",img_arr.shape)
np.save(args.to_fn,img_arr)
