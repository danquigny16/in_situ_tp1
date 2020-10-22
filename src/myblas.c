#include "myblas.h"

////////////////////////////////////////////////////////////////////////////////
// Produit scalaire

/**
Effectue le produit scalaire entre les vecteurs X et Y
@param N Nombre d'éléments des vecteurs X et Y
@param X Vecteur X
@param incX leading dimension de  X
@param Y Vecteur Y
@param incY leading dimension de  Y
@return : Le résultat du produit scalaire
*/
double my_ddot(const int N, const double *X, const int incX, const double *Y, const int incY){
  int res = 0;

  for (int i = 0; i < N; i++){
    res += X[i * incX] * Y[i * incY];
  }

  return res;
}

/**
Effectue le produit scalaire entre les vecteurs X et Y
@param N Nombre d'éléments des vecteurs X et Y
@param X Vecteur X
@param incX leading dimension de  X
@param Y Vecteur Y
@param incY leading dimension de  Y
@return : Le résultat du produit scalaire
*/
double my_ddot_unroll(const int N, const double *X, const int incX, const double *Y, const int incY){
  int res = 0;

  // Unroll sur 2 instructions
  for (int i = 0; i < N; i+=2){
    res += X[i * incX] * Y[i * incY];
    res += X[(i + 1) * incX] * Y[(i + 1) * incY];
  }

  // On n'oublie le dernier élément si N est impaire
  if (N % 2 != 0){
    res += X[(N - 1) * incX] * Y[(N - 1) * incY];
  }

  return res;
}

////////////////////////////////////////////////////////////////////////////////
// Produit de matrices

/**
Effectue : C <- alpha * (t)A * (t)B + beta * C
Avec les tailles : -A : M*K
                   -B : K*N
                   -C : M*N
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
Ici on fait les 3 boucle for dans l'ordre (k,i,j)
@param Order : Indique si les matrices A, B et C sont stockées en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : Nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param Alpha : Scalaire alpha
@param A : Matrice A
@param lda : Leading dimension de A
@param B : Matrice B
@param ldb : Leading dimension de B
@param beta : Scalaire beta
@param C : Matrice C
@param ldc : Leading dimension de C
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
Avec les tailles : -A : M*K
                   -B : K*N
                   -C : M*N
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
Ici on fait les 3 boucle for dans l'ordre (i,j,k)
@param Order : Indique si les matrices A, B et C sont stockées en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : Nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param Alpha : Scalaire alpha
@param A : Matrice A
@param lda : Leading dimension de A
@param B : Matrice B
@param ldb : Leading dimension de B
@param beta : Scalaire beta
@param C : Matrice C
@param ldc : Leading dimension de C
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
Avec les tailles : -A : M*K
                   -B : K*N
                   -C : M*N
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
Ici on fait les 3 boucle for dans l'ordre (j,i,k)
@param Order : Indique si les matrices A, B et C sont stockées en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : Nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param Alpha : Scalaire alpha
@param A : Matrice A
@param lda : Leading dimension de A
@param B : Matrice B
@param ldb : Leading dimension de B
@param beta : Scalaire beta
@param C : Matrice C
@param ldc : Leading dimension de C
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

/**
Effectue : C <- alpha * (t)A * (t)B + beta * C
Avec les tailles : -A : M*K
                   -B : K*N
                   -C : M*N
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
Ici on fait les 3 boucle for dans l'ordre (k,j,i)
@param Order : Indique si les matrices A, B et C sont stockées en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : Nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param Alpha : Scalaire alpha
@param A : Matrice A
@param lda : Leading dimension de A
@param B : Matrice B
@param ldb : Leading dimension de B
@param beta : Scalaire beta
@param C : Matrice C
@param ldc : Leading dimension de C
*/
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

/**
Effectue : C <- alpha * (t)A * (t)B + beta * C
Avec les tailles : -A : M*K
                   -B : K*N
                   -C : M*N
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
Ici on fait les 3 boucle for dans l'ordre (j,i,k), et on ajoute un unroll au niveau de la triple boucle for pour tenter d'opitimiser le code
@param Order : Indique si les matrices A, B et C sont stockées en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : Nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param Alpha : Scalaire alpha
@param A : Matrice A
@param lda : Leading dimension de A
@param B : Matrice B
@param ldb : Leading dimension de B
@param beta : Scalaire beta
@param C : Matrice C
@param ldc : Leading dimension de C
*/
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
      for (int k = 0; k < M; k += 4){
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
Avec les tailles : -A : M*K
                   -B : K*N
                   -C : M*N
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
On se contente de prendre le meilleur des algo précédents (my_dgemm_scalaire_*)
@param Order : Indique si les matrices A, B et C sont stockées en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : Nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param Alpha : Scalaire alpha
@param A : Matrice A
@param lda : Leading dimension de A
@param B : Matrice B
@param ldb : Leading dimension de B
@param beta : Scalaire beta
@param C : Matrice C
@param ldc : Leading dimension de C
*/
void my_dgemm_scalaire(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                       const int M, const int N, const int K,
                       const double alpha, const double *A, const int lda,
                       const double *B, const int ldb,
                       const double beta, double *C, const int ldc){
  // On choisit l'ordre avec les meilleurs performances
  my_dgemm_scalaire_kji(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

/**
Effectue : C <- alpha * (t)A * (t)B + beta * C
Avec les tailles : -A : M*K
                   -B : K*N
                   -C : M*N
Ici l'énoncé nous demande :
                   -de ne prendre en compte que le cas CblasColMajor
                   -de ne prendre en compte que le cas m == k == n
                   -de ne faire que : C <- tA * B
On procède par bloc, et on se sert du "my_dgemm_scalaire" pour faire les bloc
@param Order : Indique si les matrices A, B et C sont stockées en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param TransB : Indique si on doit prendre la matrice B tel quel ou sa transposé
@param M : Nombre de ligne de A / nombre de ligne de C
@param N : nombre de colonne de B / nombre de colonne de C
@param K : nombre de colonne de A / nombre de ligne de B
@param Alpha : Scalaire alpha
@param A : Matrice A
@param lda : Leading dimension de A
@param B : Matrice B
@param ldb : Leading dimension de B
@param beta : Scalaire beta
@param C : Matrice C
@param ldc : Leading dimension de C
*/
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

/**
On effectue l'opération Y <- alpha * X + Y
@param N : Taille des vecteurs X et Y et de la matrice carré A
@param Alpha : Scalaire alpha
@param X : Vecteur X
@param incX : Leading dimension de X
@param Y : Vecteur Y
@param incY : Leading dimension de Y
@return : Le résultat du produit scalaire
*/
void my_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY){
  for (int i = 0; i < N; i++){
    Y[i * incY] += alpha * X[i * incX];
  }
}

/**
On effectue l'opération Y <- alpha * (t)A * X + beta * Y
@param Order : Indique si la matrice A est stocké en CblasRowMajor ou en CblasColMajor
@param TransA : Indique si on doit prendre la matrice A tel quel ou sa transposé
@param M : Nombre de ligne de A (taille de Y si A est CblasNoTrans, taille de X sinon)
@param N : nombre de colonne de A (taille de X si A est CblasNoTrans, taille de Y sinon)
@param Alpha : Scalaire alpha
@param A : Matrice A
@param lda : Leading dimension de A
@param X : Vecteur X
@param incX : Leading dimension de X
@param beta : Scalaire beta
@param Y : Vecteur Y
@param incY : Leading dimension de Y
*/
void my_dgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
              const double alpha, const double *A, const int lda,
              const double *X, const int incX, const double beta,
              double *Y, const int incY){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor
  if (order != CblasColMajor || TransA == CblasConjTrans){
    printf("erreur dans \"my_dgemv\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Cas matrice A pris tel quel
  if (TransA == CblasNoTrans){
    for (int i = 0; i < M; i++){
      Y[i * incY] *= beta;
      double sum = 0.0;
      for (int j = 0; j < N; j++){
        sum += A[i + j * lda] * X[j * incX];
      }
      Y[i * incY] += alpha * sum;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Cas matrice A transposée
  else if (TransA == CblasTrans){
    for (int i = 0; i < N; i++){
      Y[i * incY] *= beta;
      double sum = 0.0;
      for (int j = 0; j < M; j++){
        sum += A[j + i * lda] * X[j * incX];
      }
      Y[i * incY] += alpha * sum;
    }
  }
}

/**
On effectue l'opération A <- alpha * X * tY + A
@param Order : Indique si la matrice A est stocké en CblasRowMajor ou en CblasColMajor
@param M : Nombre de ligne de A et taille de X
@param N : nombre de colonne de A et taille de Y
@param Alpha : Scalaire alpha
@param X : Vecteur X
@param incX : Leading dimension de X
@param Y : Vecteur Y
@param incY : Leading dimension de Y
@param A : Matrice A
@param lda : Leading dimension de A
*/
void my_dger(const enum CBLAS_ORDER order, const int M, const int N,
             const double alpha, const double *X, const int incX,
             const double *Y, const int incY,
             double *A, const int lda){
   //////////////////////////////////////////////////////////////////////////////
   // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor
  if (order != CblasColMajor){
    printf("erreur dans \"my_dger\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Bonne ordre de boucle pour un CblasColMajor
  for (int j = 0; j < N; j++){
    for (int i = 0; i < M; i++){
      A[i + j * lda] += alpha * X[i * incX] * Y[j * incY];
    }
  }
}
