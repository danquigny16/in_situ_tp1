#include "util.h"
#include "myblas.h"

////////////////////////////////////////////////////////////////////////////////
// Produit scalaire

/**
Effectue le produit scalaire entre les vecteurs "X" et "Y"
@param N Nombre d'éléments des vecteurs "X" et "Y"
@param X Vecteur "X"
@param incX Incrément entre chaque élément de "X"
@param Y Vecteur "Y"
@param incX Incrément entre chaque élément de "Y"
@return Le résultat du produit scalaire
*/
double my_ddot(const int N, const double *X, const int incX, const double *Y, const int incY){
  int res = 0;

  for (int i = 0; i < N; i++){
    res += X[i*incX] * Y[i*incY];
  }

  return res;
}

double my_ddot_unroll(const int N, const double *X, const int incX, const double *Y, const int incY){
  int res = 0;

  for (int i = 0; i < N; i+=2){
    res += X[i*incX] * Y[i*incY];
    res += X[(i+1)*incX] * Y[(i+1)*incY];
  }

  return res;
}
////////////////////////////////////////////////////////////////////////////////
// Produit de matrices

/**
Effectue : C <- alpha * (t)A * (t)B + beta * C
Avec les tailles : -A : m*k
                   -B : k*n
                   -C : m*n
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
Ici on fait les 3 boucle for dans l'ordre (k,i,j)
@param Order : Indique si les matrices A, B et C sont stocké en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param alpha : voir formule
@param A : voir formule
@param lda : leading dimension de A
@param B : voir formule
@param ldb : leading dimension de B
@param beta : voir formule
@param C : voir formule
@param ldc : leading dimension de C
*/
void my_dgemm_scalaire_kij(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                       const int M, const int N, const int K,
                       const double alpha, const double *A, const int lda,
                       const double *B, const int ldb,
                       const double beta, double *C, const int ldc){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor, qu'on prend la transposé de A,
  // qu'on laisse B tel quel, que l'on manipule des matrices carrés m*m, que alpha vaut 1 et beta 0,
  // qu
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || alpha != 1 || beta != 0){
    printf("erreur dans \"my_dgemm_scalaire\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // On effectue le produit et la somme de matrices (alpha * tA) * B + (beta * C)
  // i parcourt les lignes de tA (donc les colonnes de A)
  // j les colonnes de B
  // k les lignes de B (et les colonnes de tA donc les lignes de A)
  for (int k = 0; k < M; k++){
    for (int i = 0; i < M; i++){
      for (int j = 0; j < M; j++){
        // C = A * B ---> C[i + j * ldc] += A[i + k * lda] * B[k + j * ldb];
        // tA(i, j) = A(j, i)
        // C = tA * b ---> C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
        C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
      }
    }
  }
}

/**
Effectue : C <- alpha * (t)A * (t)B + beta * C
Avec les tailles : -A : m*k
                   -B : k*n
                   -C : m*n
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
Ici on fait les 3 boucle for dans l'ordre (i,j,k)
@param Order : Indique si les matrices A, B et C sont stocké en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param alpha : voir formule
@param A : voir formule
@param lda : leading dimension de A
@param B : voir formule
@param ldb : leading dimension de B
@param beta : voir formule
@param C : voir formule
@param ldc : leading dimension de C
*/
void my_dgemm_scalaire_ijk(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                       const int M, const int N, const int K,
                       const double alpha, const double *A, const int lda,
                       const double *B, const int ldb,
                       const double beta, double *C, const int ldc){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor, qu'on prend la transposé de A,
  // qu'on laisse B tel quel, que l'on manipule des matrices carrés m*m, que alpha vaut 1 et beta 0,
  // qu
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || alpha != 1 || beta != 0){
    printf("erreur dans \"my_dgemm_scalaire\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // On effectue le produit et la somme de matrices (alpha * tA) * B + (beta * C)
  // i parcourt les lignes de tA (donc les colonnes de A)
  // j les colonnes de B
  // k les lignes de B (et les colonnes de tA donc les lignes de A)
  for (int i = 0; i < M; i++){
    for (int j = 0; j < M; j++){
      for (int k = 0; k < M; k++){
        // C = A * B ---> C[i + j * ldc] += A[i + k * lda] * B[k + j * ldb];
        // tA(i, j) = A(j, i)
        // C = tA * b ---> C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
        C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
      }
    }
  }
}

/**
Effectue : C <- alpha * (t)A * (t)B + beta * C
Avec les tailles : -A : m*k
                   -B : k*n
                   -C : m*n
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
Ici on fait les 3 boucle for dans l'ordre (j,i,k)
@param Order : Indique si les matrices A, B et C sont stocké en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param alpha : voir formule
@param A : voir formule
@param lda : leading dimension de A
@param B : voir formule
@param ldb : leading dimension de B
@param beta : voir formule
@param C : voir formule
@param ldc : leading dimension de C
*/
void my_dgemm_scalaire_jik(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                       const int M, const int N, const int K,
                       const double alpha, const double *A, const int lda,
                       const double *B, const int ldb,
                       const double beta, double *C, const int ldc){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor, qu'on prend la transposé de A,
  // qu'on laisse B tel quel, que l'on manipule des matrices carrés m*m, que alpha vaut 1 et beta 0,
  // qu
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || alpha != 1 || beta != 0){
    printf("erreur dans \"my_dgemm_scalaire\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // On effectue le produit et la somme de matrices (alpha * tA) * B + (beta * C)
  // i parcourt les lignes de tA (donc les colonnes de A)
  // j les colonnes de B
  // k les lignes de B (et les colonnes de tA donc les lignes de A)
  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      for (int k = 0; k < M; k++){
        // C = A * B ---> C[i + j * ldc] += A[i + k * lda] * B[k + j * ldb];
        // tA(i, j) = A(j, i)
        // C = tA * b ---> C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
        C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
      }
    }
  }
}

void my_dgemm_scalaire_kji(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                       const int M, const int N, const int K,
                       const double alpha, const double *A, const int lda,
                       const double *B, const int ldb,
                       const double beta, double *C, const int ldc){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor, qu'on prend la transposé de A,
  // qu'on laisse B tel quel, que l'on manipule des matrices carrés m*m, que alpha vaut 1 et beta 0,
  // qu
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || alpha != 1 || beta != 0){
    printf("erreur dans \"my_dgemm_scalaire\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // On effectue le produit et la somme de matrices (alpha * tA) * B + (beta * C)
  // i parcourt les lignes de tA (donc les colonnes de A)
  // j les colonnes de B
  // k les lignes de B (et les colonnes de tA donc les lignes de A)
  for (int k = 0; k < M; k++){
    for (int j = 0; j < M; j++){
      for (int i = 0; i < M; i++){
        // C = A * B ---> C[i + j * ldc] += A[i + k * lda] * B[k + j * ldb];
        // tA(i, j) = A(j, i)
        // C = tA * b ---> C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
        C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
      }
    }
  }
}

void my_dgemm_scalaire_jik_unroll(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                       const int M, const int N, const int K,
                       const double alpha, const double *A, const int lda,
                       const double *B, const int ldb,
                       const double beta, double *C, const int ldc){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor, qu'on prend la transposé de A,
  // qu'on laisse B tel quel, que l'on manipule des matrices carrés m*m, que alpha vaut 1 et beta 0,
  // qu
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || alpha != 1 || beta != 0){
    printf("erreur dans \"my_dgemm_scalaire\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // On effectue le produit et la somme de matrices (alpha * tA) * B + (beta * C)
  // i parcourt les lignes de tA (donc les colonnes de A)
  // j les colonnes de B
  // k les lignes de B (et les colonnes de tA donc les lignes de A)
  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      for (int k = 0; k < M; k+=4){
        // C = A * B ---> C[i + j * ldc] += A[i + k * lda] * B[k + j * ldb];
        // tA(i, j) = A(j, i)
        // C = tA * b ---> C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
        C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
        C[i + j * ldc] += A[(k+1) + i * lda] * B[(k+1) + j * ldb];
        C[i + j * ldc] += A[(k+2) + i * lda] * B[(k+2) + j * ldb];
        C[i + j * ldc] += A[(k+3) + i * lda] * B[(k+3) + j * ldb];
      }
    }
  }
}



/**
Effectue : C <- alpha * (t)A * (t)B + beta * C
Avec les tailles : -A : m*k
                   -B : k*n
                   -C : m*n
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
@param Order : Indique si les matrices A, B et C sont stocké en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param alpha : voir formule
@param A : voir formule
@param lda : leading dimension de A
@param B : voir formule
@param ldb : leading dimension de B
@param beta : voir formule
@param C : voir formule
@param ldc : leading dimension de C
*/
void my_dgemm_scalaire(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                       const int M, const int N, const int K,
                       const double alpha, const double *A, const int lda,
                       const double *B, const int ldb,
                       const double beta, double *C, const int ldc){
  // On choisit l'ordre avec les meilleurs performances
  my_dgemm_scalaire_kji(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

void my_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
              const int M, const int N, const int K,
              const double alpha, const double *A, const int lda,
              const double *B, const int ldb,
              const double beta, double *C, const int ldc){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor, qu'on prend la transposé de A,
  // qu'on laisse B tel quel, que l'on manipule des matrices carrés m*m, que alpha vaut 1 et beta 0,
  // qu
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || alpha != 1 || beta != 0){
    printf("erreur dans \"my_dgemm_scalaire\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  // on définit une taille de bloc pour le produit de matrice par bloc, en supposant M multiple de bloc_size
  // sinon pb avec my_dgemm_scalaire, qui doit prendre des matrices carrés
  int bloc_size = 10;

  //////////////////////////////////////////////////////////////////////////////
  // On effectue le produit et la somme de matrices (alpha * tA) * B + (beta * C)
  // i parcourt les lignes de tA (donc les colonnes de A)
  // j les colonnes de B
  // k les lignes de B (et les colonnes de tA donc les lignes de A)
  // C = A * B ---> C[i + j * ldc] += A[i + k * lda] * B[k + j * ldb];
  // tA(i, j) = A(j, i)
  // C = tA * b ---> C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
  for (int k = 0; k < M; k += bloc_size){
    for (int j = 0; j < M; j += bloc_size){
      for (int i = 0; i < M; i += bloc_size){
        const double *NEW_A = A + (k + i * lda);
        const double *NEW_B = B + (k + j * ldb);
        double *NEW_C = C + (i + j * ldc);
        my_dgemm_scalaire(Order, TransA, TransB, bloc_size, bloc_size, bloc_size, alpha, NEW_A, lda, NEW_B, ldb, beta, NEW_C, ldc);
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
// Complétion de la bibliothèque blas

///////////////////////////////////////////////////////////////////
// On effectue l'opération Y = Y + alpha*X avec X et Y des vecteurs
void my_cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY){
  for(int i = 0; i < N; ++i){
    Y[i*incY] += alpha * X[i*incX];
  }
}

/////////////////////////////////////////////////////////////////////
// On effectue l'opération Y = Y + A * X avec Y et X des vecteurs et A une matrice
// ne tiens pas compte pour l'instant de alpha et beta
void my_cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY){
  if (order != CblasColMajor || TransA != CblasNoTrans || alpha != 1 || beta != 1 || M != N){
    printf("erreur dans \"my_cblas_dgemv\" : condition de l'énoncé non respecté\n");
    exit(0);
  }
  for(int i=0; i<M; i++){
    for(int j=0; j<M; j++){
      Y[i*incY] += A[i+j*lda] * X[j*incX];
    }
  }
}

////////////////////////////////////////////////////////////////////////
// On effectue l'opération A = A + alpha * X * tY avec A une matrice et X et Y des vecteurs
void my_cblas_dger(const enum CBLAS_ORDER order, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda){
  if (order != CblasColMajor){
    printf("erreur dans \"my_cblas_dger\" : condition de l'énoncé non respecté\n");
    exit(0);
  }
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      A[i+j*lda] += alpha * X[i*incX] * Y[j*incY];
    }
  }
}