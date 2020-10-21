#include "util.h"

////////////////////////////////////////////////////////////////////////////////
// Affichage

/**
Affiche la matrice A
@param M : hauteur
@param N : largeur
@param A : matrice, vu comme un tableau contenant les colonnes de manière continue
@param lda : Leading dimension (si A est une sous-matrice d'une matrice M, lda est la hauteur de M (car on stocke par colonne))
@param flux : flux de sortie
*/
void affiche(int M, int N, double *A, int lda, FILE * flux){
  // ligne
  for (int i = 0; i < M; i++){
    // colonne
    for (int j = 0; j < N; j++){
      fprintf(flux, "%f  ", A[i + j * lda]);
    }
    // fin de la ligne
    fprintf(flux, "\n");
  }
}

/**
Affiche le vecteur X
@param N : Taille
@param X : vecteur
@param lda : Leading dimension (si X est un vecteur d'une ligne matrice M, lda est la hauteur de M (car on stocke par colonne))
@param flux : flux de sortie
*/
void affiche_vecteur(int N, double *X, int lda, FILE * flux){
  for (int i = 0; i < N; i++){
    fprintf(flux, "%f", X[i * lda]);
    fprintf(flux, "\n");
  }
}


////////////////////////////////////////////////////////////////////////////////
// Allocation et libération des matrices et des vecteurs

/**
Alloue la mémoire pour une matrice de dim M*N, et renvoie un pointeur vers cette matrice
@param M : hauteur
@param N : largeur
@return : Un pointeur sur la matrice
*/
double * matrice(int M, int N){
  double * A = malloc(M * N * sizeof(double));
  return A;
}

/**
Alloue la mémoire pour un vecteur de dim N, et renvoie un pointeur vers ce vecteur
@param N : hauteur
@return : Un pointeur sur le vecteur
*/
double * vecteur(int N){
  double * X = malloc(N * sizeof(double));
  return X;
}

/**
Libère la mémoire pour la matrice A
@param A : La matrice
*/
void free_matrice(double *A){
  free(A);
}

/**
Libère la mémoire pour le vecteur X
@param X : La matrice
*/
void free_vecteur(double *X){
  free(X);
}


////////////////////////////////////////////////////////////////////////////////
// Initialisation et modification des matrices et des vecteurs

/**
Fixe la valeur de l'élément en position (M, N) à val
@param M : ligne
@param N : colonne
@param A : matrice
@param val : nouvel valeur de l'élément en position (M, N)
*/
void fixer_elt_matrice(int M, int N, int lda, double *A, double val){
  A[M + N * lda] = val;
}

/**
Fixe la valeur de l'élément en position N à val
@param N : numéro de l'élément à modifier
@param X : matrice
@param val : nouvel valeur de l'élément en position (M, N)
*/
void fixer_elt_vecteur(int N, int lda, double *X, double val){
  X[N * lda] = val;
}

/**
Renvoie la valeur de l'élément en position (M, N) de la matrice A
@param M : ligne
@param N : colonne
@param A : matrice
@param val : nouvel valeur de l'élément en position (M, N)
@return : L'élément en position (M, N) de la matrice A
*/
double elt_matrice(int M, int N, int lda, double *A){
  return A[M + N * lda];
}

/**
Renvoie la valeur de l'élément en position N du vecteur V
@param N : numéro de l'élément à modifier
@param X : matrice
@param val : nouvel valeur de l'élément en position (M, N)
@return : L'élément en position N du vecteur V
*/
double elt_vecteur(int N, int lda, double *X){
  return X[N * lda];
}

/**
Transforme la matrice en une matrice nulle
@param M : hauteur
@param N : largeur
@param A : matrice
*/
void matrice_zero(int M, int N, int lda, double *A){
  // ligne
  for (int i = 0; i < M; i++){
    // colonne
    for (int j = 0; j < N; j++){
      A[i + j * lda] = 0.0;
    }
  }
}

/**
Transforme le vecteur en un vecteur nulle
@param N : Taille
@param X : vecteur
*/
void vecteur_zero(int N, int lda, double *X){
  for (int i = 0; i < N; i++){
    X[i * lda] = 0.0;
  }
}

/**
initialise la matrice en allant de 1.0 à M*N, de gauche à droite puis de haut en bas
@param M : hauteur
@param N : largeur
@param A : matrice
*/
void init_matrice(int M, int N, int lda, double *A){
  double val = 1.0;
  // ligne
  for (int i = 0; i < M; i++){
    // colonne
    for (int j = 0; j < N; j++){
      A[i + j * lda] = val;
      val++;
    }
  }
}

/**
initialise le vecteur en allant de 1.0 à N, de haut en bas
@param N : Taille
@param X : vecteur
*/
void init_vecteur(int N, int lda, double *X){
  double val = 1.0;
  for (int i = 0; i < N; i++){
    X[i * lda] = val;
    val++;
  }
}

/**
Initialise la matrice pour le test de résolution AX = B, l'initialisation 1 engendre des divisions par 0
@param M : hauteur
@param N : largeur
@param A : matrice
*/
void init_2_matrice(int M, int N, int lda, double *A){
  // Initialisation habituelle
  init_matrice(M, N, lda, A);

  // Augmentation des valeur de la diagonale pour éviter les valeurs nulles
  for (int i = 0; i < M; i++){
    A[i + i * lda] *= 20.0;
  }
}
