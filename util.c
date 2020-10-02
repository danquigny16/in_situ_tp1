#include "util.h"

/**
Affiche la matrice a
@param m : hauteur
@param n : largeur
@param a : matrice, vu comme un tableau contenant les colonnes de manière continue
@param lda : leading dimension (si a est une sous-matrice d'une matrice M, lda est la hauteur de M (car on stocke par colonne))
@param flux : flux de sortie
*/
void affiche(int m, int n, int * a, int lda, FILE * flux){
  // ligne
  for (int i = 0; i < m; i++){
    // colonne
    for (int j = 0; j < n; j++){
      fprintf(flux, "%d", a[j * lda + i]);
    }
    // fin de la ligne
    fprintf(flux, "\n");
  }
}

/**
Alloue la mémoire pour une matrice de dim m*n, et renvoie un pointeur vers cette matrice
@param m : hauteur
@param n : largeur
*/
int * matrice(int m, int n){
  int * mat = malloc(m * n * sizeof(int));
  return mat;
}

/**
Alloue la mémoire pour un vecteur de dim n, et renvoie un pointeur vers ce vecteur
@param n : hauteur
*/
int * vecteur(int n){
  int * vec = malloc(n * sizeof(int));
  return vec;
}
