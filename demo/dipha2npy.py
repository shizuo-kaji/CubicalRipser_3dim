#%%
import struct
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt

#%%
parser = argparse.ArgumentParser("Convert DIPHA format to Numpy and vice versa")
parser.add_argument('from_fn')
parser.add_argument('to_fn')
args = parser.parse_args()

# %%
if ".complex" in args.from_fn:
    dat = open(args.from_fn,'rb').read()
    magic,tp,sz,dim = struct.unpack_from("qqqq",dat,0)
    sp = struct.unpack_from("q"*dim,dat,8*4)
    df = np.array(struct.unpack_from("d"*sz,dat,8*(4+dim)))
    df = df.reshape(sp)
    print(df.shape,df.min(),df.mean(),df.max())
    np.save(args.to_fn,df)
elif ".npy" in args.from_fn:
    dat = np.load(args.from_fn)
    print(dat.shape,dat.min(),dat.mean(),dat.max())
    sz = np.prod(dat.shape)
    dim = len(dat.shape)
    with open(args.to_fn,"wb") as fh:
        fh.write(struct.pack("q"*4,8067171840,1,sz,dim))
        fh.write(struct.pack("q"*dim,*dat.shape))
        fh.write(struct.pack("d"*sz,*dat.transpose((2,1,0)).flatten()))

# %%
#plt.imshow(df[:,:,0])

# %%
