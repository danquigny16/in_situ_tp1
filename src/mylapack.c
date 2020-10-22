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
  if (Order != CblasColMajor || ipiv != NULL || M != N){
    printf("erreur dans \"my_dgetf2\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Factorisation LU
  for (int i = 0; i < M; i++){
    for (int j = i + 1; j < N; j++){
      double diag = A[i + i * lda];
      if (diag == 0.0){
        printf("Division par zero durant la decomposition LU (my_dgetf2)\n");
      }
      A[j + i * lda] /= diag;
      for (int k = i + 1; k < N; k++){
        A[j + k * lda] -= A[j + i * lda] * A[i + k * lda];
      }
    }
  }
}

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
void my_dgetrf(const enum CBLAS_ORDER Order, int M, int N, double *A, int lda, int *ipiv){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor
  if (Order != CblasColMajor || ipiv != NULL || M != N){
    printf("erreur dans \"my_dgetf2\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Factorisation LU par block

  // On définit une taille de bloc
  int bloc_size = 5;
  // On suppose N multiple de bloc_size
  int nb_iter = N / bloc_size;

  for (int k = 0; k < nb_iter; k++){
    ////////////////////////////////////////////////////////////////////////////
    // Factorisation AKK = LKK * UKK
    int K = k * bloc_size;
    double *AKK = A + (K + K * lda);
    my_dgetf2(Order, bloc_size, bloc_size, AKK, lda, NULL);

    ////////////////////////////////////////////////////////////////////////////
    // Résolution AIK = LIK * UKK
    for (int i = k + 1; i < nb_iter; i++){
      int I = i * bloc_size;
      // UKK dans AKK
      double *AIK = A + (I + K * lda);
      my_dtrsm(Order, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, bloc_size, bloc_size, 1, AKK, lda, AIK, lda);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Résolution AKJ = LKK * UKJ
    for (int j = k + 1; j < nb_iter; j++){
      int J = j * bloc_size;
      // UKK dans AKK
      double *AKJ = A + (K + J * lda);
      my_dtrsm(Order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, bloc_size, bloc_size, 1, AKK, lda, AKJ, lda);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Faire AIJ -= LIK * UKJ
    for (int i = k + 1; i < nb_iter; i++){
      int I = i * bloc_size;
      for (int j = k + 1; j < nb_iter; j++){
        int J = j * bloc_size;
        // Sous-matrices
        double *LIK = A + (I + K * lda);
        double *UKJ = A + (K + J * lda);
        double *AIJ = A + (I + J * lda);
        my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, bloc_size, bloc_size, bloc_size, -1, LIK, lda, UKJ, lda, 0, AIJ, lda);
      }
    }
  }
}

/**
On effectue la résolution du système AX = alpha * B avec A une matrice triangulaire, on stocke le résultat dans B
Ne gère pas les cas où il y a des division par 0
@param Order : Indique si la matrice A est stocké en CblasRowMajor ou en CblasColMajor
@param Side : Indique de quel côté est la matrice A = LU (AX = B ou XA = B)
@param Uplo : Indique si A est une matrice triangulaire supérieure ou inférieure
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param Diag : Indique si A est une matrice unitriangulaire (que des 1 sur la diagonale)
@param M : Nombre de ligne de B et X (si Side == CblasLeft, A est M * M)
@param N : nombre de colonne de B et X (si Side == CblasRight, A est N * N)
@param alpha : Scalaire alpha
@param A : Matrice A
@param lda : Leading dimension de A
@param B : Matrice B
@param lda : Leading dimension de B
*/
void my_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
              const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
              const enum CBLAS_DIAG Diag, const int M, const int N,
              const double alpha, const double *A, const int lda,
              double *B, const int ldb){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor
  if (Order != CblasColMajor || TransA != CblasNoTrans || alpha != 1 ||
      (Uplo == CblasLower && Diag != CblasUnit) || ((Uplo == CblasLower && Diag != CblasUnit))){
    printf("erreur dans \"my_dtrsm\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // L(UX) = LY = B
  // L matrice unitriangulaire inférieure
  if (Side == CblasLeft && Uplo == CblasLower && Diag == CblasUnit){
    // Résolution de L*Y1 = B1, L*Y2 = B2, ...
    for (int k = 0; k < N; k++){
      // colonne
      for (int j = 0; j < M; j++){
        // ligne
        for (int i = j + 1; i < M; i++){
          B[i + k * ldb] -= A[i + j * lda] * B[j + k * ldb];
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // UX = Y
  // U matrice triangulaire supérieure
  else if (Uplo == CblasUpper && Diag == CblasNonUnit){
    // Résolution de U*X1 = Y1, U*X2 = Y2, ...
    for (int k = 0; k < N; k++){
      // colonne
      for (int j = M - 1; j >= 0; j--){
        double diag = A[j + j * lda];
        if (diag == 0.0){
          printf("Division par zero durant la résolution Ux = y (my_dtrsm)\n");
        }
        B[j + k * ldb] /= diag;
        // ligne
        for (int i = j - 1; i >= 0; i--){
          B[i + k * ldb] -= A[i + j * lda] * B[j + k * ldb];
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // (XL)U = YU = B       (ici X, Y et B vecteur en ligne)
  // U matrice triangulaire supérieure
  else if (Side == CblasRight && Uplo == CblasUpper && Diag == CblasNonUnit){
    // Résolution de Y1*U = B1, Y2*U = B2, ...
    for (int k = 0; k < M; k++){
      // ligne
      for (int i = 0; i < N; i++){
        // colonne
        for (int j = 0; j < i; j++){
          B[k + i * ldb] -= A[j + i * lda] * B[k + j * ldb];
        }
        double diag = A[i + i * lda];
        if (diag == 0.0){
          printf("Division par zero durant la résolution Ux = y (my_dtrsm)\n");
        }
        B[k + i * ldb] /= diag;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // XL = Y       (ici X et Y vecteur en ligne)
  // L matrice unitriangulaire inférieure
  if (Side == CblasRight && Uplo == CblasLower && Diag == CblasUnit){
    // Résolution de X1*L = Y1, X2*L = Y2, ...
    for (int k = 0; k < M; k++){
      // ligne
      for (int i = N - 1; i >= 0; i--){
        // colonne
        for (int j = N - 1; j > i; j--){
          B[k + i * ldb] -= A[j + i * lda] * B[k + j * ldb];
        }
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
  //my_dgetrf(Order, N, N, A, lda, NULL);
  my_dtrsm(Order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, N, 1, 1, A, lda, B, ldb);
  my_dtrsm(Order, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, 1, 1, A, lda, B, ldb);
}
