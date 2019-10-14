#!/usr/bin/python3
#%%
import numpy
data = [ [0, 1, 1], [1, 0, 1], [1,1,0]]
a = numpy.array(data, numpy.dtype('f8'))
numpy.save("test.npy", a[:,:,numpy.newaxis])


#%%
