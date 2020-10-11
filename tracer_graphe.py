#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

arg = sys.argv
if (len(sys.argv) != 2):
    sys.exit(0)

file = open(arg[1], "r")

X = []
Y = []

lines = file.read().splitlines()

for line in lines:
    line_data = line.split(' ')
    X.append(int(line_data[0]))
    Y.append(float(line_data[1]))

file.close()

plt.plot(X, Y)
plt.title("Performance du produit scalaire my_ddot()")
plt.xlabel("Taille des vecteurs")
plt.ylabel("Performance en Mflop/s")
plt.savefig("graphe/graphe.png")
plt.show()
plt.close()
