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

////////////////////////////////////////////////////////////////////////////////

double my_ddot(const int N, const double *X, const int incX, const double *Y, const int incY);

////////////////////////////////////////////////////////////////////////////////

void my_dgemm_scalaire_kij(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                           const int M, const int N, const int K,
                           const double alpha, const double *A, const int lda,
                           const double *B, const int ldb,
                           const double beta, double *C, const int ldc);
void my_dgemm_scalaire_ijk(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                           const int M, const int N, const int K,
                           const double alpha, const double *A, const int lda,
                           const double *B, const int ldb,
                           const double beta, double *C, const int ldc);
void my_dgemm_scalaire_jik(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                           const int M, const int N, const int K,
                           const double alpha, const double *A, const int lda,
                           const double *B, const int ldb,
                           const double beta, double *C, const int ldc);
void my_dgemm_scalaire_jik_unroll(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                           const int M, const int N, const int K,
                           const double alpha, const double *A, const int lda,
                           const double *B, const int ldb,
                           const double beta, double *C, const int ldc);
void my_dgemm_scalaire(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                       const int M, const int N, const int K,
                       const double alpha, const double *A, const int lda,
                       const double *B, const int ldb,
                       const double beta, double *C, const int ldc);
void my_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
              const int M, const int N, const int K,
              const double alpha, const double *A, const int lda,
              const double *B, const int ldb,
              const double beta, double *C, const int ldc);

////////////////////////////////////////////////////////////////////////////////

void my_cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
void my_cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);
void my_cblas_dger(const enum CBLAS_ORDER order, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda);

////////////////////////////////////////////////////////////////////////////////

void my_dgetf2(const enum CBLAS_ORDER Order, int m, int n, double* a, int lda);
void my_cblas_dtrsm(const enum CBLAS_ORDER Order, const int M, const int N, const double *A, const int lda, double *B, const int ldb);
void my_cblas_dgesv(const enum CBLAS_ORDER Order, const int N, double *A, const int lda, double *B, const int ldb);
