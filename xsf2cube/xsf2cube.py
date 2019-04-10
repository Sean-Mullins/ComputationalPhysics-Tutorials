#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 2019



@author: Sean Mullins   Mar 27, 2019
"""
import time
import ase
from ase.visualize import view
import argparse


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




coords = ase.io.read('AtomicCoordinates/ch3.xyz', index=None, format='xyz')
xsf_file = ase.io.read('sample_xsf_files/ch3.XSF', index=None, format='xsf')
cube_file = ase.io.read('sample_cube_files/ch3.RHO.UP.cube', index=None, format='cube')

#print(xsf_file)
#view(xsf_file)

view(coords)

#view(coords, viewer='avogadro', data = cube_file)

#ase.io.write('ch3.cube', xsf_file, 'cube')


elapsed_time = time.time() - start_time
print('Elapsed time: ', elapsed_time)
