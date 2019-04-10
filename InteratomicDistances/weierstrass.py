# © 2019 TheFlyingKeyboard and released under MIT License
# theflyingkeyboard.net

from math import cos

from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
from pylab import rcParams

step = 0.0001  # the smaller the more precise the plot is
a = 0.5
b = 15
n = 100

intervalBegin = 0
intervalEnd = 1

data = []

for x in range(int((1 / step) * abs(intervalEnd - intervalBegin))):
    output = 0

    for i in range(n):
        output += pow(a, i) * cos(pow(b, i) * i * x * step)

    data.append(output)

rcParams['figure.figsize'] = 16, 9

plt.gca().xaxis.grid(True)
plt.gca().yaxis.grid(True)

plt.gca().get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: format(x * step, ',')))

plt.plot(data)
plt.show()