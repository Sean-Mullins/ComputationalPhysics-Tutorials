#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 2019

Need to separate the atom types so we can calculate and write distances between
between each element type.

@author: Sean Mullins   Mar 27,2019
"""
import time
import ase
from ase import Atoms
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





coords = ase.io.read('AtomicCoordinates/Au60Mg12C60H60-symm.xyz', index=None, format='xyz')
#print(coords.get_positions())

distances = Atoms.get_all_distances(coords)
#print(distances)




elapsed_time = time.time() - start_time
print('Elapsed time: ', elapsed_time)