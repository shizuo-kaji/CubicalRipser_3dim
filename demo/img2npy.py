#!/usr/bin/env python
# -*- coding: utf-8 -*-
#%%
import struct
import numpy as np
import argparse
import os
from PIL import Image
import re,struct
import shutil
from scipy.ndimage.morphology import distance_transform_edt
from skimage.filters import threshold_otsu
from skimage.transform import rescale
from skimage import io

try:
    import nrrd
except:
    print("Install nrrd first by: pip install pynrrd")
    exit()


num = lambda val : int(re.sub("\\D", "", val+"0"))
dtype ={None: None, "uint8": np.uint8, "uint16": np.uint16, "float": np.float32, "double": np.float64}

#%%
parser = argparse.ArgumentParser("Convert image file to Numpy array")
parser.add_argument('input', nargs="*", help="multiple files of the same dimension would be stacked to a 3D image. If directory is specified, combine all the files in the directory.")
parser.add_argument('output', help="output filename (.npy or .nrrd or .raw or .pgm or .dcm or .jpg)")
parser.add_argument('--scaling_factor','-sf', type=float, default=1)
parser.add_argument('--forceSpacing', '-fs', type=float, help='scale dicom to match the specified spacing (pixel size)')
parser.add_argument('--tile','-tl', type=int, default=1)
parser.add_argument('--transform', '-tr', choices=[None,'distance','signed_distance','distance_inv','signed_distance_inv','radial','radial_inv','geodesic','geodesic_inv','upward','downward'], help="apply transform")
parser.add_argument('--transpose', '-tp', type=int, nargs='*', default=None, help='axis order (argument to numpy.transpose)')
parser.add_argument('--origin', '-o', type=int, nargs='*', default=(0,0), help='origin for the radial transformation (z,y,x)')
parser.add_argument('--threshold', '-th', type=float, default=None, help='binarisation threshold for distance transform')
parser.add_argument('--shift_value', '-sv', default=0, type=int, help='pixel values are shifted before computing PH')
parser.add_argument('--dtype','-d', type=str, default=None, choices=dtype.keys())
parser.add_argument('--sort','-s', action='store_true', help="Sort file names before stacking")
#parser.add_argument('--zsplit','-z', action='store_true', help="save one file for each slice along z-axis")
parser.add_argument('--zrange', type=int, nargs=2, help="select specific z-slices")
args = parser.parse_args()

# %%
if os.path.isdir(args.input[0]):
    is_input_dir = True
    fns = os.listdir(args.input[0])
    args.input = [os.path.join(args.input[0],f) for f in fns]
else:
    is_input_dir = False

if args.sort:
    args.input.sort(key=num)

frfn,ext = os.path.splitext(args.input[0])
ext = ext.lower()
if ext == ".dcm":
    try:
        import pydicom as dicom
    except:
        print("Install pydicom first by: pip install pydicom")
        exit()

try:
    from tqdm import tqdm
except:
    tqdm = lambda x: x

images = []
for ffn in args.input:
    #print("reading {}".format(ffn))
    fn,ext = os.path.splitext(ffn)
    if ext == ".npy":
        im = np.load(ffn)
    elif ext == ".npz":
        imz = np.load(ffn)
        im = imz[imz.files[0]]
    elif ext == ".dcm":
        ref_dicom_in = dicom.read_file(ffn, force=True)
        if not hasattr(ref_dicom_in.file_meta, "TransferSyntaxUID"):
            ref_dicom_in.file_meta.TransferSyntaxUID = dicom.uid.ImplicitVRLittleEndian ## We have to specify the right one here.
        if ref_dicom_in.file_meta.TransferSyntaxUID != dicom.uid.ImplicitVRLittleEndian:
            ref_dicom_in.file_meta.TransferSyntaxUID = dicom.uid.ImplicitVRLittleEndian
            ref_dicom_in.is_little_endian = True
            ref_dicom_in.is_implicit_VR = True
        dt=ref_dicom_in.pixel_array.dtype
        im = ref_dicom_in.pixel_array + ref_dicom_in.RescaleIntercept
        if args.forceSpacing is not None:
            scaling = float(ref_dicom_in.PixelSpacing[0])/args.forceSpacing
            im = rescale(im,scaling,mode="reflect",preserve_range=True)
    elif ext == ".csv":
        im = np.loadtxt(ffn,delimiter=",")
    elif ext == ".nrrd":
        im, header = nrrd.read(ffn, index_order='C')
    elif ext == ".complex": # DIPHA
        dat = open(args.input,'rb').read()
        magic,tp,sz,dim = struct.unpack_from("qqqq",dat,0)
        sp = struct.unpack_from("q"*dim,dat,8*4) # offset = 8byte x 4
        im = np.array(struct.unpack_from("d"*sz,dat,8*(4+dim))).reshape(sp)
    else:
        im = np.array(Image.open(ffn).convert('L'))

    images.append(im)

img_arr = np.squeeze(np.stack(images,axis=0))
print("input shape: ",img_arr.shape, "values {}--{}, dtype: {}".format(im.min(),im.max(),im.dtype))

### processing
if args.transpose is not None:
    img_arr = img_arr.transpose(args.transpose)

# select slices
if args.zrange is not None:
    img_arr = img_arr[slice(*args.zrange)]

# size reduction
if args.scaling_factor != 1:
    #from scipy.ndimage import zoom
    #img_arr = zoom(img_arr, args.scaling_factor)
    if np.issubdtype(img_arr.dtype, np.integer) and np.ptp(img_arr)<2:
        img_arr = img_arr.astype(np.bool_)
    img_arr = rescale(img_arr,args.scaling_factor,mode="reflect",preserve_range=True) #,anti_aliasing=True)

# repeating periodically
if args.tile > 1:
    img_arr = np.tile(img_arr, [args.tile]*len(img_arr.shape))

# thresholding and transform
if args.transform is not None:
    if args.threshold is None:
        img_arr = (img_arr >= threshold_otsu(img_arr))
    else:
        img_arr = (img_arr >= args.threshold)

    if 'distance' in args.transform: # distance from the background
        if '_inv' in args.transform:
            img_arr = ~img_arr
        im = distance_transform_edt(img_arr)
        if 'signed' in args.transform:
            im -= distance_transform_edt(~img_arr)
        img_arr = im
    elif args.transform == 'signed_distance':
        img_arr = distance_transform_edt(img_arr)-distance_transform_edt(~img_arr)
    elif args.transform in ['downward','upward']:
        null_idx = img_arr == 0
        ## height transform
        if len(img_arr.shape) == 3: #(z,y,x)
            h = np.arange(img_arr.shape[0]).reshape(-1,1,1)
        else:
            h = np.arange(img_arr.shape[0]).reshape(-1,1)
        if args.transform=='upward':
            #h = np.max(h) - h
            h = -h
        img_arr = (img_arr * h)
        img_arr[null_idx] = np.max(img_arr)
    elif 'radial' in args.transform:
        null_idx = img_arr == 0
        h = np.linalg.norm(np.stack(np.meshgrid(*map(range,img_arr.shape),indexing='ij'),axis=-1)-np.array(args.origin), axis=-1)
        img_arr = (img_arr * h)
        if args.transform=='radial_inv':
            #img_arr = np.max(img_arr) - img_arr
            img_arr *= -1  ## background pixels are 0
        else:
            img_arr[null_idx] = np.max(h)
    elif 'geodesic' in args.transform:
        try:
            import skfmm
        except:
            print("install skfmm by 'pip install scikit-fmm'")
            exit()
        roi = np.ones(img_arr.shape)
        if len(roi.shape)==3:
            roi[args.origin[0],args.origin[1],args.origin[2]] = 0
        else:
            roi[args.origin[0],args.origin[1]] = 0
        img_arr = skfmm.distance(np.ma.MaskedArray(roi,~img_arr))  # mask=True specifies the obstacle
        if '_inv' in args.transform:
            img_arr *= -1
        img_arr = img_arr.filled(fill_value=img_arr.max())
        #print(roi.min(),roi.max())


if args.shift_value:
    img_arr += args.shift_value

if args.dtype is not None:
    img_arr = img_arr.astype(dtype[args.dtype])

### save
print("output ",args.output, " shape: ",img_arr.shape, " dtype: ",img_arr.dtype, " values: {} -- {}".format(img_arr.min(),img_arr.max()))

tofn,ext = os.path.splitext(args.output)
ext = ext.lower()
if ext == ".raw":  # for use with cubicle
    data = img_arr.astype(np.uint16).flatten()
    with open(args.output, 'wb') as out:
        for v in data:
            out.write(struct.pack('H', v))   # B: uchar, H: ushort
elif ext == ".pgm":   # for use with diamorse
    if len(img_arr.shape)==3:
        with open(args.output, 'wb') as outfile:
            for z in range(img_arr.shape[2]):
                fname = tofn+"_{:0>4}.pgm".format(z)
                Image.fromarray(img_arr[:,:,z].astype(np.uint8),mode='L').save(fname)
                with open(fname, 'rb') as infile:
                    shutil.copyfileobj(infile,outfile)
                os.remove(fname)
    else:
        Image.fromarray(img_arr.astype(np.uint8),mode='L').save(args.output)
elif ext == ".nrrd":
    nrrd.write(args.output, img_arr, index_order='C')
elif ext == ".npy":
    np.save(args.output,img_arr)
elif ext == ".dcm": ## TODO: correct header info
    if len(img_arr.shape) == 3:
        os.makedirs(tofn,exist_ok=True)
        for i in range(img_arr.shape[0]):
            ref_dicom_in.PixelData = (img_arr[i]-ref_dicom_in.RescaleIntercept).astype(dt).tobytes()
            ref_dicom_in.Rows, ref_dicom_in.Columns = img_arr[i].shape
            ref_dicom_in.save_as(os.path.join(tofn,"{:0>4}.dcm".format(i)), write_like_original=False)
    else:
        ref_dicom_in.PixelData = (img_arr-ref_dicom_in.RescaleIntercept).astype(dt).tobytes()
        ref_dicom_in.Rows, ref_dicom_in.Columns = img_arr.shape
        ref_dicom_in.save_as(args.output, write_like_original=False)

else: # a file for each slice
    if len(img_arr.shape) == 3:
        os.makedirs(tofn,exist_ok=True)
        for i in range(img_arr.shape[0]):
            io.imsave(os.path.join(tofn,"_{:0>4}.jpg".format(i)),img_arr[i])
    else:
        io.imsave(args.output,img_arr)
