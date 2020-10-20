#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

arg = sys.argv
if (len(sys.argv) != 8):
    print("Il n'y a pas le bon nombre d'arguments.\n")
    sys.exit(0)

file_vect = open(arg[1], "r")
file_ijk = open(arg[2], "r")
file_jik = open(arg[3], "r")
file_kij = open(arg[4], "r")
file_kji = open(arg[5], "r")
file_bloc = open(arg[6], "r")
file_unroll = open(arg[7], "r")

def file_to_XY(file):
    X = []
    Y = []
    for line in file.read().splitlines():
        line_data = line.split(' ')
        X.append(int(line_data[0]))
        Y.append(float(line_data[1]))
    return X,Y

X_vect,Y_vect = file_to_XY(file_vect)
X_ijk,Y_ijk = file_to_XY(file_ijk)
X_jik,Y_jik = file_to_XY(file_jik)
X_kij,Y_kij = file_to_XY(file_kij)
X_kji,Y_kji = file_to_XY(file_kji)
X_bloc,Y_bloc = file_to_XY(file_bloc)
X_unroll,Y_unroll = file_to_XY(file_unroll)

for file in [file_ijk, file_jik, file_kij, file_kji, file_bloc, file_vect, file_unroll]:
    file.close()

plt.plot(X_vect, Y_vect)
plt.title("Performance du produit scalaire my_ddot()")
plt.xlabel("Taille des vecteurs")
plt.ylabel("Performance en Mflop/s")
plt.savefig("graphe/graphe.png")
plt.show()
plt.close()
