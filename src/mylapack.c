#include "mylapack.h"


///////////////////////////////////////////////////////////////////////////////
// Factorisation LU et résolution de système triangulaire

/**
On effectue la factorisation LU de la matrice A
Ne gère pas les cas où il y a des division par 0
@param Order : Indique si la matrice A est stocké en CblasRowMajor ou en CblasColMajor
@param M : Nombre de ligne de A
@param N : nombre de colonne de A
@param A : Matrice A
@param lda : Leading dimension de A
@param ipiv : Pivot pour résoudre le problème des divisions par 0
*/
void my_dgetf2(const enum CBLAS_ORDER Order, int M, int N, double *A, int lda, int *ipiv){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor
  if (Order != CblasColMajor || ipiv != NULL){
    printf("erreur dans \"my_dgetf2\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Factorisation LU
  for (int i = 0; i < M; i++){
    for (int j = i + 1; j < N; j++){
      double diag = A[i + i * lda];
      if (diag == 0.0){
        printf("ouille ouille ouille : division par zero durant la decomposition LU (my_dgetf2)\n");
      }
      A[j + i * lda] /= diag;
      for (int k = i + 1; k < N; k++){
        A[j + k * lda] -= A[j + i * lda] * A[i + k * lda];
      }
    }
  }
}

/**
On effectue la résolution du système AX = B avec A une matrice triangulaire, on stocke le résultat dans B
Ne gère pas les cas où il y a des division par 0
@param Order : Indique si la matrice A est stocké en CblasRowMajor ou en CblasColMajor
@param Side : Indique de quel côté est la matrice A = LU (AX = B ou XA = B)
@param Uplo : Indique si A est une matrice triangulaire supérieure ou inférieure
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param Diag : Indique si A est une matrice unitriangulaire (que des 1 sur la diagonale)
@param M : Nombre de ligne de B et X (nombre de ligne de A si Side == CblasLeft, nombre de colonne de A sinon)
@param N : nombre de colonne de B et X (nombre de colonne de A si Side == CblasLeft, nombre de ligne de A sinon)
@param alpha : Scalaire alpha
@param A : Matrice A
@param lda : Leading dimension de A
@param B : Vecteur B
@param lda : Leading dimension de B
*/
void my_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor
  if (Order != CblasColMajor || Side != CblasLeft || TransA != CblasNoTrans || alpha != 1 || M != N ||
      (Uplo == CblasLower && Diag != CblasUnit) || ((Uplo == CblasLower && Diag != CblasUnit))){
    printf("erreur dans \"my_dtrsm\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // L(UX) = LY = B
  // L matrice unitriangulaire inférieure
  // colonne
  if (Uplo == CblasLower && Diag == CblasUnit){
    for (int j = 0; j < M; j++){
      // ligne
      for (int i = j + 1; i < M; i++){
        B[i * ldb] -= A[i + j * lda] * B[j * ldb];
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // UX = Y
  // U matrice triangulaire supérieure
  else if (Uplo == CblasUpper && Diag == CblasNonUnit){
    // colonne
    for (int j = N - 1; j >= 0; j--){
      double diag = A[j + j * lda];
      if (diag == 0.0){
        printf("aie aie aie : division par zero durant la résolution Ux = y (my_dtrsm)\n");
      }
      B[j * ldb] /= diag;
      // ligne
      for (int i = j - 1; i >= 0; i--){
        B[i * ldb] -= A[i + j * lda] * B[j * ldb];
      }
    }
  }
}

/**
On effectue la résolution du système AX = B (ne gère pas les cas où il y a des division par 0)
@param Order : Indique si la matrice A est stocké en CblasRowMajor ou en CblasColMajor
@param N : Taille de la matrice A et du vecteur B (ou nombre de ligne de la matrice B)
@param nrhs : Nombre de colonne de B
@param A : Matrice A
@param lda : Leading dimension de A
@param ipiv : Pivot pour résoudre le problème des divisions par 0
@param B : Vecteur B
@param lda : Leading dimension de B
*/
void my_dgesv(const enum CBLAS_ORDER Order, const int N, int nrhs, double *A, int lda, int *ipiv, double *B, int ldb){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor
  if (Order != CblasColMajor || nrhs != 1 || ipiv != NULL){
    printf("erreur dans \"my_dtrsm\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Résolution de AX = B
  my_dgetf2(Order, N, N, A, lda, NULL);
  my_dtrsm(Order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, N, N, 1, A, lda, B, ldb);
  my_dtrsm(Order, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, N, 1, A, lda, B, ldb);
}
