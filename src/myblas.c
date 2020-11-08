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
Ici on rajoute un unroll de 2
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
  for (int i = 0; i < N - 1; i += 2){
    res += X[i * incX] * Y[i * incY];
    res += X[(i + 1) * incX] * Y[(i + 1) * incY];
  }

  // On n'oublie le dernier élément si N est impaire
  if (N % 2 != 0){
    res += X[(N - 1) * incX] * Y[(N - 1) * incY];
  }

  return res;
}

/**
Effectue le produit scalaire entre les vecteurs X et Y
Ici on se sert de la parallelisation de openmp
@param N Nombre d'éléments des vecteurs X et Y
@param X Vecteur X
@param incX leading dimension de  X
@param Y Vecteur Y
@param incY leading dimension de  Y
@return : Le résultat du produit scalaire
*/
double my_ddot_openmp(const int N, const double *X, const int incX, const double *Y, const int incY){
  int res = 0;

  #pragma omp parallel for
  for (int i = 0; i < N; i++){
    res += X[i * incX] * Y[i * incY];
  }

  return res;
}

/**
Effectue le produit scalaire entre les vecteurs X et Y
Ici on se sert des opérations SIMD avec openmp
@param N Nombre d'éléments des vecteurs X et Y
@param X Vecteur X
@param incX leading dimension de  X
@param Y Vecteur Y
@param incY leading dimension de  Y
@return : Le résultat du produit scalaire
*/
double my_ddot_openmp_simd(const int N, const double *X, const int incX, const double *Y, const int incY){
  int res = 0;

  #pragma omp simd
  for (int i = 0; i < N; i++){
    res += X[i * incX] * Y[i * incY];
  }

  return res;
}

/**
Effectue le produit scalaire entre les vecteurs X et Y
Ici on se sert des opérations SIMD avx2
@param N Nombre d'éléments des vecteurs X et Y
@param X Vecteur X
@param incX leading dimension de  X
@param Y Vecteur Y
@param incY leading dimension de  Y
@return : Le résultat du produit scalaire
*/
double my_ddot_avx2(const int N, const double *X, const int incX, const double *Y, const int incY){
  int res = 0;

  //////////////////////////////////////////////////////////////////////////////

  if (incX == 1 && incY == 1){
    // res_tab est un tableau pour contenir notre petit vecteur de double (maximum 4 double avec avx2-fma), qui contiendra les resultats
    __m256d res_tab = _mm256_set1_pd(0.0);
    // res_tab_tmp servira à stocker les résultat des multiplication, on additionnera ensuite son contenu avec res_tab pour avoir le résultat
    __m256d res_tab_tmp = _mm256_set1_pd(0.0);

    // On fait nos opération par paquet de 4 double, si N n'est pas un multiple de 4, on fait le reste à la fin
    int size_vector = N / 4;
    int remain = N % 4;

    // On interprete notre tableau de double en tableau de __mm256d
    __m256d *X_vector = (__m256d *) X;
    __m256d *Y_vector = (__m256d *) Y;

    // On fait nos calculs de  "X[i * incX] * Y[i * incY]" par paquet de 4
    for (int i = 0; i < size_vector; i++, X_vector++, Y_vector++){
      res_tab_tmp = _mm256_mul_pd(*X_vector, *Y_vector);
      res_tab = _mm256_add_pd(res_tab, res_tab_tmp);
    }

    // On compacte notre petit vecteur de double pour récupérer le résultat final
    res += ((double *) &res_tab)[0];
    res += ((double *) &res_tab)[1];
    res += ((double *) &res_tab)[2];
    res += ((double *) &res_tab)[3];

    // On n'oublie pas de faire le reste si N n'était pas un multiple de 4
    for (int i = N - remain; i < N; i++){
      res += X[i * incX] * Y[i * incY];
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  else{
    // res_tab est un tableau pour contenir notre petit vecteur de double (maximum 4 double avec avx2-fma), qui contiendra les resultats
    __m256d res_tab = _mm256_set1_pd(0.0);
    // res_tab_tmp servira à stocker les résultat des multiplication, on additionnera ensuite son contenu avec res_tab pour avoir le résultat
    __m256d res_tab_tmp = _mm256_set1_pd(0.0);

    // On fait nos opération par paquet de 4 double, si N n'est pas un multiple de 4, on fait le reste à la fin
    int size_vector = N / 4;
    int remain = N % 4;

    // On prendra nos variables par paquet de 4
    __m256d X_vector;
    __m256d Y_vector;

    // On fait nos multiplication "X[i * incX] * Y[i * incY]" par paquet de 4, mais nos additions ne peuvent pas être
    // faite avec SIMD, on fait donc un unroll de 4
    for (int i = 0; i < size_vector; i++){
      X_vector = _mm256_set_pd(X[i], X[i + incX], X[i + 2 * incX], X[i + 3 * incX]);
      Y_vector = _mm256_set_pd(Y[i], Y[i + incX], Y[i + 2 * incX], Y[i + 3 * incX]);

      res_tab_tmp = _mm256_mul_pd(X_vector, Y_vector);
      res_tab = _mm256_add_pd(res_tab, res_tab_tmp);
    }

    // On compacte notre petit vecteur de double pour récupérer le résultat final
    res += ((double *) &res_tab)[0];
    res += ((double *) &res_tab)[1];
    res += ((double *) &res_tab)[2];
    res += ((double *) &res_tab)[3];

    // On n'oublie pas de faire le reste si N n'était pas un multiple de 4
    for (int i = N - remain; i < N; i++){
      res += X[i * incX] * Y[i * incY];
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  return res;
}

/**
Effectue le produit scalaire entre les vecteurs X et Y
Ici on se sert des opérations SIMD avx2 et fma
@param N Nombre d'éléments des vecteurs X et Y
@param X Vecteur X
@param incX leading dimension de  X
@param Y Vecteur Y
@param incY leading dimension de  Y
@return : Le résultat du produit scalaire
*/
double my_ddot_avx2_fma(const int N, const double *X, const int incX, const double *Y, const int incY){
  int res = 0;

  //////////////////////////////////////////////////////////////////////////////

  if (incX == 1 && incY == 1){
    // res_tab est un tableau pour contenir notre petit vecteur de double (maximum 4 double avec avx2-fma), qui contiendra les resultats
    __m256d res_tab = _mm256_set1_pd(0.0);

    // On fait nos opération par paquet de 4 double, si N n'est pas un multiple de 4, on fait le reste à la fin
    int size_vector = N / 4;
    int remain = N % 4;

    // On interprete notre tableau de double en tableau de __mm256d
    __m256d *X_vector = (__m256d *) X;
    __m256d *Y_vector = (__m256d *) Y;

    // On fait nos calculs de  "X[i * incX] * Y[i * incY]" par paquet de 4
    for (int i = 0; i < size_vector; i++, X_vector++, Y_vector++){
      res_tab = _mm256_fmadd_pd(*X_vector, *Y_vector, res_tab);
    }

    // On compacte notre petit vecteur de double pour récupérer le résultat final
    res += ((double *) &res_tab)[0];
    res += ((double *) &res_tab)[1];
    res += ((double *) &res_tab)[2];
    res += ((double *) &res_tab)[3];

    // On n'oublie pas de faire le reste si N n'était pas un multiple de 4
    for (int i = N - remain; i < N; i++){
      res += X[i * incX] * Y[i * incY];
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  else{
    // res_tab est un tableau pour contenir notre petit vecteur de double (maximum 4 double avec avx2-fma), qui contiendra les resultats
    __m256d res_tab = _mm256_set1_pd(0.0);

    // On fait nos opération par paquet de 4 double, si N n'est pas un multiple de 4, on fait le reste à la fin
    int size_vector = N / 4;
    int remain = N % 4;

    // On prendra nos variables par paquet de 4
    __m256d X_vector;
    __m256d Y_vector;

    // On fait nos multiplication "X[i * incX] * Y[i * incY]" par paquet de 4, mais nos additions ne peuvent pas être
    // faite avec SIMD, on fait donc un unroll de 4
    for (int i = 0; i < size_vector; i++){
      X_vector = _mm256_set_pd(X[i], X[i + incX], X[i + 2 * incX], X[i + 3 * incX]);
      Y_vector = _mm256_set_pd(Y[i], Y[i + incX], Y[i + 2 * incX], Y[i + 3 * incX]);

      res_tab = _mm256_fmadd_pd(X_vector, Y_vector, res_tab);
    }

    // On compacte notre petit vecteur de double pour récupérer le résultat final
    res += ((double *) &res_tab)[0];
    res += ((double *) &res_tab)[1];
    res += ((double *) &res_tab)[2];
    res += ((double *) &res_tab)[3];

    // On n'oublie pas de faire le reste si N n'était pas un multiple de 4
    for (int i = N - remain; i < N; i++){
      res += X[i * incX] * Y[i * incY];
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  return res;
}

/**
Effectue le produit scalaire entre les vecteurs X et Y
Ici on se sert des opérations SIMD avx2 et fma ainsi que de la parallelisation openmp
@param N Nombre d'éléments des vecteurs X et Y
@param X Vecteur X
@param incX leading dimension de  X
@param Y Vecteur Y
@param incY leading dimension de  Y
@return : Le résultat du produit scalaire
*/
double my_ddot_avx2_fma_openmp(const int N, const double *X, const int incX, const double *Y, const int incY){
  int res = 0;

  //////////////////////////////////////////////////////////////////////////////

  if (incX == 1 && incY == 1){
    // res_tab est un tableau pour contenir notre petit vecteur de double (maximum 4 double avec avx2-fma), qui contiendra les resultats
    __m256d res_tab = _mm256_set1_pd(0.0);

    // On fait nos opération par paquet de 4 double, si N n'est pas un multiple de 4, on fait le reste à la fin
    int size_vector = N / 4;
    int remain = N % 4;

    // On interprete notre tableau de double en tableau de __mm256d
    __m256d *X_vector = (__m256d *) X;
    __m256d *Y_vector = (__m256d *) Y;

    // On fait nos calculs de  "X[i * incX] * Y[i * incY]" par paquet de 4
    #pragma omp parallel for
    for (int i = 0; i < size_vector; i++){
      res_tab = _mm256_fmadd_pd(*X_vector, *Y_vector, res_tab);
      X_vector++;
      Y_vector++;
    }

    // On compacte notre petit vecteur de double pour récupérer le résultat final
    res += ((double *) &res_tab)[0];
    res += ((double *) &res_tab)[1];
    res += ((double *) &res_tab)[2];
    res += ((double *) &res_tab)[3];

    // On n'oublie pas de faire le reste si N n'était pas un multiple de 4
    for (int i = N - remain; i < N; i++){
      res += X[i * incX] * Y[i * incY];
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  else{
    // res_tab est un tableau pour contenir notre petit vecteur de double (maximum 4 double avec avx2-fma), qui contiendra les resultats
    __m256d res_tab = _mm256_set1_pd(0.0);

    // On fait nos opération par paquet de 4 double, si N n'est pas un multiple de 4, on fait le reste à la fin
    int size_vector = N / 4;
    int remain = N % 4;

    // On prendra nos variables par paquet de 4
    __m256d X_vector;
    __m256d Y_vector;

    // On fait nos multiplication "X[i * incX] * Y[i * incY]" par paquet de 4, mais nos additions ne peuvent pas être
    // faite avec SIMD, on fait donc un unroll de 4
    #pragma omp parallel for
    for (int i = 0; i < size_vector; i++){
      X_vector = _mm256_set_pd(X[i], X[i + incX], X[i + 2 * incX], X[i + 3 * incX]);
      Y_vector = _mm256_set_pd(Y[i], Y[i + incX], Y[i + 2 * incX], Y[i + 3 * incX]);

      res_tab = _mm256_fmadd_pd(X_vector, Y_vector, res_tab);
    }

    // On compacte notre petit vecteur de double pour récupérer le résultat final
    res += ((double *) &res_tab)[0];
    res += ((double *) &res_tab)[1];
    res += ((double *) &res_tab)[2];
    res += ((double *) &res_tab)[3];

    // On n'oublie pas de faire le reste si N n'était pas un multiple de 4
    for (int i = N - remain; i < N; i++){
      res += X[i * incX] * Y[i * incY];
    }
  }

  //////////////////////////////////////////////////////////////////////////////

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
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || beta != 0){
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
        C[i + j * ldc] += alpha * A[k + i * lda] * B[k + j * ldb];
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
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || beta != 0){
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
        C[i + j * ldc] += alpha * A[k + i * lda] * B[k + j * ldb];
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
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || beta != 0){
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
        C[i + j * ldc] += alpha * A[k + i * lda] * B[k + j * ldb];
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
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || beta != 0){
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
        C[i + j * ldc] += alpha * A[k + i * lda] * B[k + j * ldb];
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
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || beta != 0){
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
        C[i + j * ldc] += alpha * A[k + i * lda] * B[k + j * ldb];
        C[i + j * ldc] += alpha * A[(k+1) + i * lda] * B[(k+1) + j * ldb];
        C[i + j * ldc] += alpha * A[(k+2) + i * lda] * B[(k+2) + j * ldb];
        C[i + j * ldc] += alpha * A[(k+3) + i * lda] * B[(k+3) + j * ldb];
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
void my_dgemm_scalaire_kji_no_transpose(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                           const int M, const int N, const int K,
                           const double alpha, const double *A, const int lda,
                           const double *B, const int ldb,
                           const double beta, double *C, const int ldc){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor, qu'on prend la transposé de A,
  // qu'on laisse B tel quel, que l'on manipule des matrices carrés m*m, que alpha vaut 1 et beta 0,
  // qu
  if (Order != CblasColMajor || TransA != CblasNoTrans || TransB != CblasNoTrans || M != N || N != K || beta != 0){
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
        C[i + j * ldc] += alpha * A[i + k * lda] * B[k + j * ldb];
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
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor, qu'on prend la transposé de A,
  // qu'on laisse B tel quel, que l'on manipule des matrices carrés m*m, que alpha vaut 1 et beta 0,
  // qu
  if (Order != CblasColMajor || TransB != CblasNoTrans || M != N || N != K || beta != 0){
  printf("erreur dans \"my_dgemm_scalaire\" : condition de l'énoncé non respecté\n");
  exit(0);
  }

  // On choisit l'ordre avec les meilleurs performances
  if (TransA == CblasTrans){
    my_dgemm_scalaire_kji(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }
  else if (TransA == CblasNoTrans){
    my_dgemm_scalaire_kji_no_transpose(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
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
  if (Order != CblasColMajor || TransB != CblasNoTrans || M != N || N != K || beta != 0){
    printf("erreur dans \"my_dgemm_scalaire\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  // on définit une taille de bloc pour le produit de matrice par bloc, en supposant M multiple de bloc_size
  // sinon pb avec my_dgemm_scalaire, qui doit prendre des matrices carrés
  int bloc_size = 5;

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
        // Calcul des adresse des sous-matrices
        const double *NEW_A = A + (k + i * lda);
        if (TransA == CblasNoTrans){
          NEW_A = A + (i + k * lda);
        }
        const double *NEW_B = B + (k + j * ldb);
        double *NEW_C = C + (i + j * ldc);
        // dgemm sur un block
        my_dgemm_scalaire(Order, TransA, TransB, bloc_size, bloc_size, bloc_size, alpha, NEW_A, lda, NEW_B, ldb, beta, NEW_C, ldc);
      }
    }
  }
}

/**
Effectue : C <- alpha * (t)A * (t)B + beta * C
Avec les tailles : -A : M*K
                   -B : K*N
                   -C : M*N
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
void dgemm_seq_opti(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE transA, const enum CBLAS_TRANSPOSE transB,
                    const int M, const int N, const int K,
                    const double alpha, const double *A, const int lda,
                    const double *B, const int ldb,
                    const double beta, double *C, const int ldc){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor, qu'on prend la transposé de A,
  // qu'on laisse B tel quel, que l'on manipule des matrices carrés m*m, que alpha vaut 1 et beta 0,
  // qu
  if (Order != CblasColMajor){
    printf("erreur dans \"dgemm_seq_opti\" : condition de l'énoncé non respecté\n");
    exit(0);
  }

  int m, n, k;

  if (transA == CblasNoTrans){
    if (transB == CblasNoTrans){
      for (m = 0; m < M; m++){
        for (n = 0; n < N; n++){
          C[ldc * n + m] *= beta;
          for (k = 0; k < K; k++){
            C[ldc * n + m] += alpha * A[lda * k + m] * B[ldb * n + k];
          }
        }
      }
    }
    else {
      for (m = 0; m < M; m++){
        for (n = 0; n < N; n++){
          C[ldc * n + m] *= beta;
          for (k = 0; k < K; k++){
            C[ldc * n + m] += alpha * A[lda * k + m] * B[ldb * k + n];
          }
        }
      }
    }
  }
  else {
    if (transB == CblasNoTrans){
      for (m = 0; m < M; m++){
        for (n = 0; n < N; n++){
          C[ldc * n + m] *= beta;
          for (k = 0; k < K; k++){
            C[ldc * n + m] += alpha * A[lda * m + k] * B[ldb * n + k];
          }
        }
      }
    }
    else {
      for (m = 0; m < M; m++){
        for (n = 0; n < N; n++){
          C[ldc * n + m] *= beta;
          for (k = 0; k < K; k++){
            C[ldc * n + m] += alpha * A[lda * m + k] * B[ldb * k + n];
          }
        }
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
