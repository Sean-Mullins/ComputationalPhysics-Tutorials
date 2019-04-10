#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 2019

1) Convert the calculation of the center of mass and translation of atoms
into a function we can call from another program.

2) Maybe round the center of mass position to 4 decimal places.

3) Translate the positions of all the atoms so C.O.M. is a origin.

4) Return the centered coordinates

5) Comment timer once function is done.

6) Create short tutorial describing the use of the function.

@author: Sean Mullins   Mar 27, 2019
"""

import time
import ase



start_time = time.time()



coords = ase.io.read('AtomicCoordinates/Au60Mg12C60H60-symm.xyz', index=None, format='xyz')
print(coords.get_center_of_mass())
#print(coords.get_positions())



elapsed_time = time.time() - start_time
print('Elapsed time: ', elapsed_time)