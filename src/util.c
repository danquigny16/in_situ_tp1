#include "util.h"


////////////////////////////////////////////////////////////////////////////////
// Affichage

/**
Affiche la matrice a
@param m : hauteur
@param n : largeur
@param a : matrice, vu comme un tableau contenant les colonnes de manière continue
@param lda : leading dimension (si a est une sous-matrice d'une matrice M, lda est la hauteur de M (car on stocke par colonne))
@param flux : flux de sortie
*/
void affiche(int m, int n, double *a, int lda, FILE * flux){
  // ligne
  for (int i = 0; i < m; i++){
    // colonne
    for (int j = 0; j < n; j++){
      fprintf(flux, "%f  ", a[i + j * lda]);
    }
    // fin de la ligne
    fprintf(flux, "\n");
  }
}

/**
Affiche le vecteur v
@param n : taille
@param v : vecteur
@param lda : leading dimension (si v est un vecteur d'une ligne matrice M, lda est la hauteur de M (car on stocke par colonne))
@param flux : flux de sortie
*/
void affiche_vecteur(int n, double *v, int lda, FILE * flux){
  for (int i = 0; i < n; i++){
    fprintf(flux, "%f", v[i * lda]);
    fprintf(flux, "\n");
  }
}


////////////////////////////////////////////////////////////////////////////////
// Allocation et libération des matrices et des vecteurs

/**
Alloue la mémoire pour une matrice de dim m*n, et renvoie un pointeur vers cette matrice
@param m : hauteur
@param n : largeur
@return Un pointeur sur la matrice
*/
double * matrice(int m, int n){
  double * mat = malloc(m * n * sizeof(double));
  return mat;
}

/**
Alloue la mémoire pour un vecteur de dim n, et renvoie un pointeur vers ce vecteur
@param m : hauteur
@return Un pointeur sur le vecteur
*/
double * vecteur(int n){
  double * vec = malloc(n * sizeof(double));
  return vec;
}

/**
Libère la mémoire pour la matrice "a"
@param a La matrice
*/
void free_matrice(double *a){
  free(a);
}

/**
Libère la mémoire pour le vecteur "a"
@param a La matrice
*/
void free_vecteur(double *v){
  free(v);
}


////////////////////////////////////////////////////////////////////////////////
// Initialisation et modification des matrices et des vecteurs

/**
Fixe la valeur de l'élément en position (m, n) à "val"
@param m : ligne
@param n : colonne
@param a : matrice
@param val : nouvel valeur de l'élément en position (m, n)
*/
void fixer_elt_matrice(int m, int n, int lda, double *a, double val){
  a[m + n * lda] = val;
}

/**
Fixe la valeur de l'élément en position "n" à "val"
@param n : numéro de l'élément à modifier
@param v : matrice
@param val : nouvel valeur de l'élément en position (m, n)
*/
void fixer_elt_vecteur(int n, int lda, double *v, double val){
  v[n * lda] = val;
}

/**
initialise la matrice en allant de 1.0 à m*n, de gauche à droite puis de haut en bas
@param m : hauteur
@param n : largeur
@param a : matrice
*/
void init_matrice(int m, int n, int lda, double *a){
  double val = 1.0;
  // ligne
  for (int i = 0; i < m; i++){
    // colonne
    for (int j = 0; j < n; j++){
      a[i + j * lda] = val;
      val++;
    }
  }
}

/**
initialise le vecteur en allant de 1.0 à n, de haut en bas
@param n : taille
@param v : vecteur
*/
void init_vecteur(int n, int lda, double *v){
  double val = 1.0;
  for (int i = 0; i < n; i++){
    v[i * lda] = val;
    val++;
  }
}


////////////////////////////////////////////////////////////////////////////////
// Opération sur les matrices et les vecteurs

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

/**
Effectue : C <- alpha * A * B + beta * C
Avec les tailles : -A : m*k
                   -B : k*n
                   -C : m*n
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
                       const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc){
  //////////////////////////////////////////////////////////////////////////////
  // Pour cette fonction on suppose dans l'énoncé qu'on est en CblasColMajor, qu'on prend la transposé de A,
  // qu'on laisse B tel quel, que l'on manipule des matrices carrés m*m, que alpha vaut 1 et beta 0,
  // qu
  if (Order != CblasColMajor || TransA != CblasTrans || TransB != CblasNoTrans || M != N || N != K || alpha != 1 || beta != 0){
    printf("erreur dans \"my_dgemm_scalaire\" : condition de l'énoncé non respecté");
    exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // On effectue le produit et la somme de matrices (alpha * tA) * B + (beta * C)
  for (int i = 0; i < M; i++){
    for (int k = 0; k < M; k++){
      for (int j = 0; j < M; j++){
        // C = A * B
        // C[i + j * ldc] += A[i + k * lda] * B[k + j * ldb];
        // tA(i, j) = A(j, i)
        // C = tA * b
        //C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
        C[i + j * ldc] += A[k + i * lda] * B[k + j * ldb];
      }
    }
  }
}

///////////////////////////////////////////////////////////////////
// On effectue l'opération Y = Y + alpha*X avec X et Y des vecteurs
void my_cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY){
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
  if (Order != CblasColMajor || TransA != CblasTrans){
    printf("erreur dans \"my_cblas_dgemv\" : condition de l'énoncé non respecté");
    exit(0);
  }
  for(int i=0; i<M; i++){
    for(int j=0; j<M; j++){
      Y[i*incY] += A[i+j*lda] * X[j*incX];
    }
  }
}

void my_cblas_dger(const enum CBLAS_ORDER order, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda){
  if (Order != CblasColMajor){
    printf("erreur dans \"my_cblas_dger\" : condition de l'énoncé non respecté");
    exit(0);
  }
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      A[i+j*lda] += alpha * X[i*incX] * Y[j*incY];
    }
  }
}
