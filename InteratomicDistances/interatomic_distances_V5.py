#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 16:11:03 2018
#    The goal is to separate the elements and then determine the interatomic distances
#    between each atom and each other atom and sort the distance by the corresponding elements
#    -> So for atom1 in element.type calc distance and append to list of atom1.type to atom2.type
#    Create list of lists with n(n+1)/2,  n = number of elements.
#    4 elements = 10 lists. {(Au-Au, Au-Mg, Au-C, Au-H),(Mg-Mg, Mg-C, Mg-H),(C-C, C-H),(H-H)}

Pending: Output of graphable data files.
Pending: Arguments allowing execution of python file and supplying .xyz file to use,
         i.e., python interatomic_distances.py -f structure.xyz
Pending: Checks for correct number of atoms (compare number in xyz file to the length of the imported database)
Pending: Other checks? Number of interatomic distances based on total number of atoms,

@author: Sean Mullins   Aug,2018
"""
import time
import pandas as pd
import numpy as np
#from scipy.spatial import distance_matrix
from scipy.spatial.distance import squareform, pdist
import argparse
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel
#from decimal import *

start_time = time.time()


parser = argparse.ArgumentParser(description="This program does stuff",
                                 usage="name_of_script ...",
                                 epilog="this is text added to the end of the help")
# argument: files to process
parser.add_argument("-f", "--file",
                    dest="file",
                    default=None,
                    help="the file to process")

args = parser.parse_args()



# Loading dataset. Need to skip first two lines of xyz file.
#df = pd.read_csv(args.file, delim_whitespace=True, header=None, skiprows=1)
df = pd.read_csv('AtomicCoordinates/Au60Mg12C60H60-symm.xyz', delim_whitespace=True, header=None, skiprows=1)
df.columns = ['Element', 'X', 'Y', 'Z']            
df.set_index('Element', inplace=True)
#print(df)
#print(df.iloc[0].name)
#print(df.loc['Au'])

# Creating list of the elements in structure
elements = df.index.unique(level='Element')
#print(len(elements))
#print(type(elements))
#print(elements)
ndarray_elements = elements.values
#print(ndarray_elements)
#   Creating the complete Distance Matrix
dist = pd.DataFrame(squareform(pdist(df.iloc[:, 0:])), columns=df.index, index=df.index)
#print(dist)

#print(dist.loc[['Au'], ['Au']])


elapsed_time = time.time() - start_time
print('Elapsed time: ', elapsed_time)
#%%
#   Creating distance matrices for each element_element pair
#   If the element pair are of the same element type, i.e. Au_Au, then the matrix is symmetric.
for first_element in range(len(elements)):
#    print(elements[first_element])
    for second_element in range(first_element, len(elements)):
#        print(elements[second_element])
        if first_element == second_element:
            dist2 = dist.loc[[elements[first_element]], [elements[second_element]]]
            dist3 = squareform(dist2)
            dist3.sort()
            test = dist3
            filename  = "./InteratomicDistances/dist_" + str(ndarray_elements[first_element]) + "_" + str(ndarray_elements[second_element]) + ".dat"            
            with open(filename, "wb") as f:
                np.savetxt(f, dist3, delimiter=",")
#            print(elements[first_element] + "_" + elements[second_element])
#            print(y)
        else:
            dist2 = dist.loc[[elements[first_element]], [elements[second_element]]]
            dist3 = dist2.values.ravel()
            dist3.sort()
            filename  = "./InteratomicDistances/dist_" + str(ndarray_elements[first_element]) + "_" + str(ndarray_elements[second_element]) + ".dat"
            with open(filename, "wb") as f:
                np.savetxt(f, dist3, delimiter=",")
#            print(filename)
#            print(len(dist3))
  

#%%
#               Initiating a gaussian
'''
1) Setting the number decimal places in the interatomic distances.

2) Can also find the max number of decimal places in the distances in order to
calculate the spacing of the gaussian.
'''
# 1)
data_flt = []
for item in test:
    data_flt.append(float("%.4f" % round(item,4)))
#print(float(data_flt[0]))
#print(min(data_flt))
    
#  2)
#rev_data_str = []
#for string in data_str:
##    print(string)
#    rev_data_str.append(string[::-1].find('.'))
#decimals = max(rev_data_str)

decimals = 4   # Set above with:   %.4f   
spacing = 10 ** decimals
#print (spacing)

#%%
#   Plotting with histogram

#plt.hist(data_flt, 100)
        
#%%
              
#print(max(test)-min(test))

sigma1 = 0.001
dx = 1/spacing
x = np.arange(min(data_flt), max(test), dx)

def gaussian(x,x0,sigma):
  return np.exp(-((x - x0)/sigma **2)/2.)

data_gaussian = gaussian(x, dx, sigma1)
result0 = np.convolve(test, data_gaussian, mode ="full")
result = np.array(result0)

print(result.max())

#plt.plot(dx, result)
#plt.show()

#%%

gauss_kernel = Gaussian1DKernel(.001)
smoothed_data_gauss = convolve(data_flt, gauss_kernel)

plt.plot(smoothed_data_gauss)
plt.show()

#%%
elapsed_time = time.time() - start_time
print('Elapsed time: ', elapsed_time)



#           Pending: Other analytics of the structure.
        #Distance to center of mass. Bond angles, dihedral angles. Symmetry

#       Testing Variance and Standard Deviation
#test_dist = Au_dist[0:150]
#print(test_dist)
#print(np.var(test_dist))
#print(np.std(test_dist))


