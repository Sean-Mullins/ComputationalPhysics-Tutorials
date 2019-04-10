# © 2019 TheFlyingKeyboard and released under MIT License
# theflyingkeyboard.net

import math
import numpy as np
import scipy
from matplotlib import pyplot as plt

def discrete_gaussian_kernel(t, n):
    return math.exp(-t) * scipy.special.iv(n, t)
ns = np.arange(-5, 5+1)
y0 = discrete_gaussian_kernel(0.1, ns)
y1 = discrete_gaussian_kernel(1, ns)
y2 = discrete_gaussian_kernel(2, ns)
y3 = discrete_gaussian_kernel(4, ns)
plt.plot(ns, y0, ns, y1, ns, y2, ns, y3)
plt.xlim([-4, 4])
plt.ylim([0, 0.7])