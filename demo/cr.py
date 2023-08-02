#!/usr/bin/env python
# -*- coding: utf-8 -*-
# test for python module
#%%
import struct
import numpy as np
import argparse
import os,time
import matplotlib.pyplot as plt
import cripser,tcripser
from PIL import Image
import re, glob
from scipy.ndimage import distance_transform_edt
from skimage.filters import threshold_otsu,threshold_niblack,threshold_sauvola
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor
from skimage.transform import rescale

try:
    from tqdm import tqdm
except:
    tqdm = lambda x: x

num = lambda val : int(re.sub("\\D", "", val+"0"))

def comp_PH(fn):
    im = np.array(Image.open(fn).convert('L'),dtype=np.float64)
    #im = dt(im)
    pd = cripser.computePH(im)
    #pd = cripser.computePH(im,maxdim=args.maxdim,top_dim=args.top_dim,embedded=args.embedded)
    np.save(os.path.splitext(fn)[0],pd)
    return

# load a series of files to form a volume
def load_vol(fns,transform=None,shift_value=None, threshold=None,threshold_upper_limit=None,scaling_factor=1,negative=False,sort=False,origin=(0,0,0)):
    print("loading {}".format(fns[0]))
    if sort:
        fns.sort(key=num)
    images = []
    for ffn in fns:
        fn,ext = os.path.splitext(ffn)
        ext = ext.lower()
        dtype = int
        if ext == ".npy":
            im = np.load(ffn)
            dtype = im.dtype
        elif ext == ".dcm":
            ref_dicom_in = dicom.read_file(ffn, force=True)
    #        ref_dicom_in.file_meta.TransferSyntaxUID = dicom.uid.ImplicitVRLittleEndian
            im = ref_dicom_in.pixel_array +ref_dicom_in.RescaleIntercept
        elif ext == ".csv":
            im = np.loadtxt(ffn,delimiter=",")
            dtype = im.dtype
        elif ext == ".nrrd":
            im, header = nrrd.read(ffn, index_order='C')
            dtype = im.dtype
        elif ext == ".complex": # DIPHA
            dat = open(ffn,'rb').read()
            magic,tp,sz,dim = struct.unpack_from("qqqq",dat,0)
            sp = struct.unpack_from("q"*dim,dat,8*4) # offset = 8byte x 4
            im = np.array(struct.unpack_from("d"*sz,dat,8*(4+dim))).reshape(sp)
        else:
            im = np.array(Image.open(ffn).convert('L'))
        images.append(im)

    img_arr = np.squeeze(np.stack(images,axis=0)) #(z,y,x)

    # pre-process
    if transform is not None:
        # binarisation
        if threshold is not None:
            if threshold_upper_limit is not None:
                img_arr = np.logical_and(img_arr >= threshold,img_arr <= threshold_upper_limit)
            else:    
                img_arr = (img_arr >= threshold)
        elif threshold_upper_limit is not None:
            img_arr = (img_arr <= threshold_upper_limit)
        else:        
            img_arr = (img_arr >= threshold_otsu(img_arr))
        
        if 'distance' in transform: # distance from the background
            if '_inv' in transform:
                img_arr = ~img_arr
            im = distance_transform_edt(img_arr)
            if 'signed' in transform:
                im -= distance_transform_edt(~img_arr)
            img_arr = im
        elif transform == 'signed_distance':
            img_arr = distance_transform_edt(img_arr)-distance_transform_edt(~img_arr)
        elif transform in ['downward','upward']:
            null_idx = img_arr == 0
            ## height transform
            if len(img_arr.shape) == 3: #(z,y,x)
                h = np.arange(img_arr.shape[0]).reshape(-1,1,1)
            else:
                h = np.arange(img_arr.shape[0]).reshape(-1,1)
            if transform=='upward':
                #h = np.max(h) - h
                h = -h
            img_arr = (img_arr * h)
            img_arr[null_idx] = np.max(img_arr)
        elif 'radial' in transform:
            null_idx = img_arr == 0
            h = np.linalg.norm(np.stack(np.meshgrid(*map(range,img_arr.shape),indexing='ij'),axis=-1)-np.array(origin), axis=-1)
            img_arr = (img_arr * h)
            if transform=='radial_inv':
                #img_arr = np.max(img_arr) - img_arr
                img_arr *= -1
            else:
                img_arr[null_idx] = np.max(h)
        elif 'geodesic' in transform:
            try:
                import skfmm
            except:
                print("install skfmm by 'pip install scikit-fmm'")
                exit()
            roi = np.ones(img_arr.shape)
            if args.origin_mask is not None:
                skl, header = nrrd.read(args.origin_mask, index_order='C')
                roi[skl>0]=0  # >0 specifies outside the object
                #print(roi.size, roi.sum(), skl.sum(),skl.min(),skl.max())
            else:
                if len(roi.shape)==3:
                    roi[args.origin[0],args.origin[1],args.origin[2]] = 0
                else:
                    roi[args.origin[0],args.origin[1]] = 0
            # compute the signed distance from the 0-contour
            img_arr = skfmm.distance(np.ma.MaskedArray(roi,~img_arr))  # outside the region (~img_arr) is masked
            if '_inv' in transform:
                img_arr *= -1
            img_arr = img_arr.filled(fill_value=img_arr.max())
            #print(roi.min(),roi.max())

    # data type
    img_arr = img_arr.astype(np.float64)

    if scaling_factor != 1:
        img_arr = rescale(img_arr,scaling_factor,order=1, mode="reflect",preserve_range=True)

    if negative:
        img_arr *= -1

    if shift_value:
        img_arr += shift_value

    print("input shape: ",img_arr.shape)
    return(img_arr,dtype)

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser("")
    parser.add_argument('input',type=str, nargs="*", help="numpy array or multiple images or directory containing multiple images")
    parser.add_argument('--output', '-o', default=None)
    parser.add_argument('--filtration', '-f', choices=['V','T'], default='V')
    parser.add_argument('--top_dim',action='store_true')
    parser.add_argument('--embedded', '-e', action='store_true', help='compute for the Alexander dual complex')
    parser.add_argument('--negative', '-n', action='store_true', help='negate the pixel values')
    parser.add_argument('--maxdim','-m', default=2,type=int)
    parser.add_argument('--sort','-s', action='store_true', help="Sort file names before stacking")
    parser.add_argument('--software',type=str,default="cubicalripser")
    parser.add_argument('--batch', '-b', type=int, default=0, help='batch computation (num of threads)')
    parser.add_argument('--imgtype', '-it', type=str, default=None)
    parser.add_argument('--transform', '-tr', choices=[None,'distance','signed_distance','distance_inv','signed_distance_inv','radial','radial_inv','geodesic','geodesic_inv','upward','downward'], help="apply transform")
    parser.add_argument('--origin', type=int, nargs='*', default=(0,0), help='origin for the radial transformation (z,y,x)')
    parser.add_argument('--origin_mask', '-om', type=str, default=None, help='nrrd file of the same dimension with the input specifying the object from which the geodesic distance will be calculated in geodesic mode')
    parser.add_argument('--shift_value', '-sv', default=None, type=float, help='pixel values are shifted before computing PH')
    parser.add_argument('--threshold', '-th', type=float, default=None, help='binarisation threshold for distance transform')
    parser.add_argument('--threshold_upper_limit', '-thu', type=float, default=None, help='binarisation threshold for distance transform')
    parser.add_argument('--scaling_factor','-sf', default=1,type=float, help="scale before computation for saving comutational cost")
    args = parser.parse_args()

    if args.batch>0:
        fns = []
        if args.imgtype is None:
            args.imgtype = "png"
        for dn in args.input:
            fns.extend(glob.glob(os.path.join(dn,"**/*.{}".format(args.imgtype)), recursive=True))
        pool = Pool(args.batch)
        #with ProcessPoolExecutor(args.batch) as executor:
        #    tqdm(executor.map(comp_PH,fns), total=len(fns))
        with tqdm(total=len(fns)) as t:
            for _ in pool.imap_unordered(comp_PH, fns):
                t.update(1)
        exit(0)

    if os.path.isdir(args.input[0]):
        if args.output is None:
            args.output = os.path.basename(os.path.normpath(args.input[0]))
        if args.imgtype is not None:
            fns = glob.glob(os.path.join(args.input[0],"**/*.{}".format(args.imgtype)), recursive=True)
        else:
            fns = os.listdir(args.input[0])
            fns = [os.path.join(args.input[0],f) for f in fns]
    else:
        fns = args.input

    fn,ext = os.path.splitext(fns[0])
    ext = ext.lower()
    if ext == ".dcm":
        try:
            import pydicom as dicom
        except:
            print("Install pydicom first by: pip install pydicom")
            exit()
    if ext == ".nrrd":
        try:
            import nrrd
        except:
            print("Install nrrd first by: pip install pynrrd")
            exit()

    img_arr, dtype = load_vol(fns,args.transform,args.shift_value,args.threshold,args.threshold_upper_limit,args.scaling_factor,args.negative,args.sort,origin=args.origin)

    # compute PH
    print("computing PH for {}..".format(fns[0]))
    start = time.time()
    if args.software=="gudhi":
        try:
            import gudhi
        except:
            print("Install gudhi first by: conda install -c conda-forge gudhi")
            exit()

        print("Computing PH with GUDHI (T-construction)")
        gd = gudhi.CubicalComplex(top_dimensional_cells=img_arr)
    #    gd.compute_persistence()
        res = np.array(gd.persistence(2,0)) # coeff = 2
        print("Betti numbers: ", gd.persistent_betti_numbers(np.inf,-np.inf))

    else: # cripser
        if args.filtration=='V':
            res = cripser.computePH(img_arr,maxdim=args.maxdim,top_dim=args.top_dim,embedded=args.embedded)
        else:
            res = tcripser.computePH(img_arr,maxdim=args.maxdim,top_dim=args.top_dim,embedded=args.embedded)
        print("Betti numbers: ", [res[res[:,0]==i].shape[0] for i in range(3)])
    #    print(res[:10])

    print ("computation took:{} [sec]".format(time.time() - start))

    # save results to file
    if args.output is not None:
        if args.output.endswith(".csv"):
            np.savetxt(args.output,res,delimiter=',',fmt='%d,%18.10f,%18.10f,%d,%d,%d,%d,%d,%d')
        else:
            np.save(args.output,res)
