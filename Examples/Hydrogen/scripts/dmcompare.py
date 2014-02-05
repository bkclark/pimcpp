import sys
from math import *
from numpy import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

a = loadtxt('e-e.dm',skiprows=12).flatten()
b = loadtxt('old.dm',skiprows=13).flatten()

plt.plot(a-b)
plt.savefig('dmcompare.png')
plt.clf()
