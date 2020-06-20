######
# test for python module
#%%
import struct
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import cripser

#%%
parser = argparse.ArgumentParser("")
parser.add_argument('input',type=str, help="numpy array")
parser.add_argument('--output',default=None)
parser.add_argument('--location',default="birth")
parser.add_argument('--top_dim',action='store_true')
parser.add_argument('--embedded',action='store_true')
parser.add_argument('--maxdim',default=2,type=int)
args = parser.parse_args()

a = np.load(args.input)
print("input shape: ", a.shape)
res = cripser.computePH(a,maxdim=args.maxdim,top_dim=args.top_dim,embedded=args.embedded,location=args.location)
print("Betti numbers: ", [res[res[:,0]==i].shape[0] for i in range(3)])
