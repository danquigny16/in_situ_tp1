#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"


////////////////////////////////////////////////////////////////////////////////

void affiche(int m, int n, double * a, int lda, FILE * flux);
void affiche_vecteur(int n, double *v, int lda, FILE * flux);

////////////////////////////////////////////////////////////////////////////////

double * matrice(int m, int n);
double * vecteur(int n);
void free_matrice(double *a);
void free_vecteur(double *v);

////////////////////////////////////////////////////////////////////////////////

void fixer_elt_matrice(int m, int n, int lda, double *a, double val);
void fixer_elt_vecteur(int n, int lda, double *v, double val);
void matrice_zero(int m, int n, int lda, double *a);
void vecteur_zero(int n, int lda, double *v);
void matrice_diag_2(int m, int n, int lda, double *a);
void init_matrice(int m, int n, int lda, double *a);
void init_vecteur(int n, int lda, double *v);
