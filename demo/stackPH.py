######
import numpy as np
import argparse
import os
from PIL import Image

#%%

# compute the PH lifetime image
def lifeVect(ph,mx,my,mz,dims=None, max_life=-1):
    if dims is None:
        if mz == 1:
            dims = [0,1]
        else:
            dims = [0,1,2]

    if max_life<0:
        max_life = np.quantile(ph[:,2]-ph[:,1],0.8)
        #print(max_life)
    life_vec = np.zeros((len(dims),mx,my,mz))
    for c in ph:
        d = int(c[0]) # dim
        if d in dims:
            life = min(c[2]-c[1],max_life)
            di = dims.index(d)
            x,y,z=int(c[3]),int(c[4]),int(c[5]) # location
            life_vec[di,x,y,z] = max(life_vec[di,x,y,z],life)

    return(np.squeeze(life_vec))


# compute the PH histogram image
def histVect(ph,mx,my,mz,min_life=-1,max_life=-1,n_life_bins=4,n_birth_bins=4,dims=None):
    if dims is None:
        if mz == 1:
            dims = [0,1]
        else:
            dims = [0,1,2]
    #
    m,M = np.quantile(ph[:,2]-ph[:,1],[0.1,0.9])
    if min_life < 0:
        min_life = m
    if max_life < 0:
        max_life = M

    min_birth, max_birth = np.quantile(ph[:,1],[0.1,0.9])
    #
    life_bins = np.linspace(min_life,max_life,n_life_bins-1)
    birth_bins =np.linspace(min_birth,max_birth,n_birth_bins-1)
    cycle_vec = np.zeros((mx,my,mz,len(dims),n_life_bins,n_birth_bins))
#    print(life_bins,birth_bins)

    # histogram
    for c in ph:
        d = int(c[0]) # dim
        if d in dims:
            b = c[1] # birth
            l = c[2]-c[1] # life
            x,y,z=int(c[3]),int(c[4]),int(c[5]) # location
            i = np.searchsorted(life_bins, l)
            j = np.searchsorted(birth_bins, b)
            cycle_vec[x,y,z,d,i,j] += 1
    
    if mz==1:
        res = cycle_vec.reshape((mx,my,-1))
        res = res.transpose(2,0,1)
    else:
        res = cycle_vec.reshape((mx,my,mz,-1))
        res = res.transpose(3,0,1,2)
    return(res)

##
if __name__ == '__main__':
    parser = argparse.ArgumentParser("create PH stacked image")
    parser.add_argument('PH_file', help="Output of Cubical Ripser in Numpy format")
    parser.add_argument('--image','-i', help="original image (in Numpy format) to be stacked")
    parser.add_argument('--dims','-d', type=int, nargs="*", default=[0,1,2], help="dimensions of cycles to be used")
    parser.add_argument('--width','-x', type=int, default=0, help="Width of the image")
    parser.add_argument('--height','-y', type=int, default=0, help="Height of the image")
    parser.add_argument('--depth','-z', type=int, default=0, help="Depth of the image")
    parser.add_argument('--max_life','-Ml', type=float, default=-1, help="maximum lifetime to be counted")
    parser.add_argument('--min_life','-ml', type=float, default=-1, help="minimum lifetime to be counted")
    parser.add_argument('--n_life_bins','-nl', type=int, default=4, help="number of bins in histogram for lifetime")
    parser.add_argument('--n_birth_bins','-nb', type=int, default=4, help="number of bins in histogram for birthtime")

    parser.add_argument('--out','-o', default="stacked.npy")
    parser.add_argument('--type','-t', default="life", choices=["life","hist"])
    args = parser.parse_args()

    # load computed PH
    ph = np.load(args.PH_file)
    if args.image:
        img = np.load(args.image)
        args.width = img.shape[0]
        args.height = img.shape[1]
        if len(img.shape)==2:
            args.depth = 1
        else:
            args.depth = img.shape[2]

    if args.width ==0 or args.height == 0 or args.depth == 0:
        print("Specify either an input image file or image width/height/depth")
        exit()

    if args.type=="hist":
        PHvec = histVect(ph,args.width,args.height,args.depth,max_life=args.max_life,n_life_bins=args.n_life_bins,
                                n_birth_bins=args.n_birth_bins,dims=args.dims)
    elif args.type=="life":
        PHvec = lifeVect(ph,args.width,args.height,args.depth,max_life=args.max_life,dims=args.dims)

    if args.image:
        out_arr = np.concatenate([np.expand_dims(img, 0),PHvec],axis=0)
    else:
        out_arr = PHvec
    
    np.save(args.out, out_arr)
    print("Output file: {}, shape: {}".format(args.out,out_arr.shape))
