#include "util.h"
#include "myblas.h"


////////////////////////////////////////////////////////////////////////////////
// Factorisation LU et résolution de système triangulaire

void my_dgetf2(const enum CBLAS_ORDER Order, int M, int N, double *A, int lda, int *ipiv);
void my_dgetrf(const enum CBLAS_ORDER Order, int M, int N, double *A, int lda, int *ipiv);
void my_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
              const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
              const enum CBLAS_DIAG Diag, const int M, const int N,
              const double alpha, const double *A, const int lda,
              double *B, const int ldb);
void my_dgesv(const enum CBLAS_ORDER Order, const int N, int nrhs, double *A, int lda, int *ipiv, double *B, int ldb);
