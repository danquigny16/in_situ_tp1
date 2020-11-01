#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

#include "cblas.h"


////////////////////////////////////////////////////////////////////////////////
// Affichage

void affiche(int M, int N, double *A, int lda, FILE * flux);
void affiche_vecteur(int N, double *X, int lda, FILE * flux);

////////////////////////////////////////////////////////////////////////////////
// Allocation et libération des matrices et des vecteurs

double * matrice(int M, int N);
double * vecteur(int N);
void free_matrice(double *A);
void free_vecteur(double *X);

////////////////////////////////////////////////////////////////////////////////
// Initialisation et modification des matrices et des vecteurs

void fixer_elt_matrice(int M, int N, int lda, double *A, double val);
void fixer_elt_vecteur(int N, int lda, double *X, double val);
double elt_matrice(int M, int N, int lda, double *A);
double elt_vecteur(int N, int lda, double *X);
void matrice_zero(int M, int N, int lda, double *A);
void vecteur_zero(int N, int lda, double *X);
void init_matrice(int M, int N, int lda, double *A);
void init_vecteur(int N, int lda, double *X);
void init_2_matrice(int M, int N, int lda, double *A);
