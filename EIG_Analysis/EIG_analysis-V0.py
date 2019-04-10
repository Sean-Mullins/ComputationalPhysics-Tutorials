# File for loading and plotting the EIG file from SIESTA.
#Version0: Sean - Jan23

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.ndimage import gaussian_filter
import time


start_time = time.time()

#Creates a list of files in the current working directory with the .EIG extension
#If there is not one .EIG file, it will raise an Error.

EIG_files = [f for f in os.listdir('.') if f.endswith('.EIG')]
#print(len(EIG_files))
if len(EIG_files) != 1:
    raise ValueError('There should be only one EIG file in the current directory')

filename = EIG_files[0]

with open(filename, 'r') as myfile:
    data = myfile.read()

# Reads the EIG file and loads it as a string, then splits string into list of strings,
# then converts each string in list into a float.

data_str = data.split()
#print(data_str)
data_flt = [float(i) for i in data_str]
eigenvalue_data = data_flt

#print(type(data))
#print(type(data_str), type(data_str[0]))
#print(type(data_flt), type(data_flt[0]))

fermi_energy = eigenvalue_data[0]
unknown1 = eigenvalue_data[1]
spin = eigenvalue_data[2]
unknown2 = eigenvalue_data[3]
unknown3 = eigenvalue_data[4]

# Removing the values that are not eigenvalues
eigenvalues = eigenvalue_data
del eigenvalues[0]; del eigenvalues[0]; del eigenvalues[0]; del eigenvalues[0]; del eigenvalues[0]
#print(eigenvalues)

ndarray_eigenvalues = np.array(eigenvalues)
#print(ndarray_eigenvalues.max())
#print(ndarray_eigenvalues.min())

# The eigenvalues are given for each spin state, need to split these up somehow
# Will do later

#%%
#               Initiating a gaussian
'''
Finding the max number of decimal places in the eigenvalues in order to
calculate the spacing of the gaussian
'''
rev_data_str = []
for str in data_str:
    rev_data_str.append(str[::-1].find('.'))
#print (rev_data_str)
#print(max(rev_data_str))
decimals = max(rev_data_str)
spacing = 10 ** decimals
#print (spacing)

#%%

sigma = 
dx = float(ndarray_eigenvalues.max() - ndarray_eigenvalues.min())*spacing

gx = np.arange(-3*sigma, 3*sigma, dx)
gaussian = np.exp(-(x/sigma)**2/2)

result = np.convolve(ndarray_eigenvalues, gaussian, mode="full")


# Creating convolution with a gaussian

#np.convolve(ndarray_eigenvalues, gaussian, mode='full')

#%%

# Plotting the eigenvalues
# the histogram of the data
#n, bins, patches = plt.hist(eigenvalues, 50, normed=1, facecolor='green', alpha=0.75)

#range = [-10.,-5., 0., 5.]
#
#plt.xlabel('Eigenvalues')
#plt.ylabel('Degeneracy')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([min(eigenvalues), fermi_energy-min(eigenvalues), 0, 0.2])
#plt.grid(True)
#plt.plot(range, ndarray_eigenvalues)
#plt.show()

# Want to determine degeneracies near the fermi energy


elapsed_time = time.time() - start_time
print('Elapsed time: ', elapsed_time)
