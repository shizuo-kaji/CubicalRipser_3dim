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
from scipy.ndimage.morphology import distance_transform_edt
from skimage.filters import threshold_otsu
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor

num = lambda val : int(re.sub("\\D", "", val))


def dt(img):
    bw_img = (img >= threshold_otsu(img))
    dt_img = distance_transform_edt(bw_img)-distance_transform_edt(~bw_img)
    return(dt_img)

def comp_PH(fn):
    im = np.array(Image.open(fn).convert('L'),dtype=np.float64)
    im = dt(im)
    pd = cripser.computePH(im)
    #pd = cripser.computePH(im,maxdim=args.maxdim,top_dim=args.top_dim,embedded=args.embedded)
    np.save(os.path.splitext(fn)[0],pd)
    return

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser("")
    parser.add_argument('input',type=str, nargs="*", help="numpy array or multiple images or directory containing multiple images")
    parser.add_argument('--output', '-o', default=None)
    parser.add_argument('--filtration', '-f', choices=['V','T'], default='V')
    parser.add_argument('--top_dim',action='store_true')
    parser.add_argument('--embedded', '-e', action='store_true')
    parser.add_argument('--maxdim','-m', default=2,type=int)
    parser.add_argument('--sort','-s', action='store_true', help="Sort file names before stacking")
    parser.add_argument('--software',type=str,default="cubicalripser")
    parser.add_argument('--batch', '-b', type=int, default=0, help='batch computation (num of threads)')
    parser.add_argument('--imgtype', '-it', type=str, default=None)
    parser.add_argument('--distance_transform', '-dt', action="store_true", help="apply distance transform before computing PH")
    args = parser.parse_args()

    if args.batch>0:
        from tqdm import tqdm
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
    if args.distance_transform:
        img_arr = dt(img_arr)
    print("input shape: ",img_arr.shape)

    start = time.time()
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
        if args.filtration=='V':
            res = cripser.computePH(img_arr,maxdim=args.maxdim,top_dim=args.top_dim,embedded=args.embedded)
        else:
            res = tcripser.computePH(img_arr,maxdim=args.maxdim,top_dim=args.top_dim,embedded=args.embedded)
        print("Betti numbers: ", [res[res[:,0]==i].shape[0] for i in range(3)])
    #    print(res[:10])

    print ("computation took:{} [sec]".format(time.time() - start))

    if args.output is not None:
        np.save(args.output,res)
