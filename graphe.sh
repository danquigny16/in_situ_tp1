#!/bin/bash

./test | grep Performance | tr -s ' ' | cut -d ' ' -f 8,10 > donnees_graphe.txt
python3 tracer_graphe.py donnees_graphe.txt
