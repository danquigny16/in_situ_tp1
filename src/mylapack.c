#include "util.h"
///////////////////////////////////////////////////////////////////////////////
// Factorisation LU

///////////////////////////////////////////////////////////////////////////////
// Factorisation LU de la matrice a
void my_dgetf2(const enum CBLAS_ORDER Order, int m, int n, double* a, int lda){
  if (Order != CblasColMajor){
    printf("erreur dans \"my_cblas_dgetf2\" : condition de l'énoncé non respecté\n");
    exit(0);
  }
  for(int i=0; i<m; i++){
    for(int j=i+1; j<n; j++){
      a[j+i*lda] /= a[i+i*lda];
      for(int k=i+1; k<n; k++){
        a[j+k*lda] -= a[j+i*lda] * a[i+k*lda];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// Résolution de système triangulaire LX = B
void my_cblas_dtrsm(const enum CBLAS_ORDER Order, const int M, const int N, const double *A, const int lda, double *B, const int ldb){
  if (Order != CblasColMajor){
    printf("erreur dans \"my_cblas_dtrsm\" : condition de l'énoncé non respecté\n");
    exit(0);
  }
  (void) ldb;
  int sum;
  B[M-1] /= A[(M-1) + (M-1) * lda];
  for(int i = M-2; i > -1; i--){
    sum = 0;
    for(int j = N-1; j > i; j--){
      sum += A[i+j*lda] * B[j];
    }
    B[i] = (B[i] - sum)/A[i+i*lda];
  }
}

//////////////////////////////////////////////////////////////////////////
// résolution de système linéaire A*X = B
void my_cblas_dgesv(const enum CBLAS_ORDER Order, const int N, double *A, const int lda, double *B, const int ldb){
     if (Order != CblasColMajor){
       printf("erreur dans \"my_cblas_dtrsm\" : condition de l'énoncé non respecté\n");
       exit(0);
     }
     my_dgetf2(Order, N, N, A, lda);
     my_cblas_dtrsm(Order, N, N, A, lda, B, ldb);
}
