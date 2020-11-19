#include "util.h"


////////////////////////////////////////////////////////////////////////////////
// Produit scalaire

double my_ddot(const int N, const double *X, const int incX, const double *Y, const int incY);
double my_ddot_unroll(const int N, const double *X, const int incX, const double *Y, const int incY);
double my_ddot_openmp(const int N, const double *X, const int incX, const double *Y, const int incY);
double my_ddot_openmp_simd(const int N, const double *X, const int incX, const double *Y, const int incY);
double my_ddot_avx2(const int N, const double *X, const int incX, const double *Y, const int incY);
double my_ddot_avx2_fma(const int N, const double *X, const int incX, const double *Y, const int incY);
double my_ddot_avx2_fma_openmp(const int N, const double *X, const int incX, const double *Y, const int incY);

////////////////////////////////////////////////////////////////////////////////
// Produit de matrices

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
void my_dgemm_scalaire_kji(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                           const int M, const int N, const int K,
                           const double alpha, const double *A, const int lda,
                           const double *B, const int ldb,
                           const double beta, double *C, const int ldc);
void my_dgemm_scalaire_jik_unroll(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                                  const int M, const int N, const int K,
                                  const double alpha, const double *A, const int lda,
                                  const double *B, const int ldb,
                                  const double beta, double *C, const int ldc);
void my_dgemm_scalaire_kji_no_transpose(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
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
void dgemm_seq_opti(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE transA, const enum CBLAS_TRANSPOSE transB,
                    const int M, const int N, const int K,
                    const double alpha, const double *A, const int lda,
                    const double *B, const int ldb,
                    const double beta, double *C, const int ldc);
void dgemm_seq_simd(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE transA, const enum CBLAS_TRANSPOSE transB,
                    const int M, const int N, const int K,
                    const double alpha, const double *A, const int lda,
                    const double *B, const int ldb,
                    const double beta, double *C, const int ldc);

////////////////////////////////////////////////////////////////////////////////
// Complétion de la bibliothèque blas

void my_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
void my_dgemv(const enum CBLAS_ORDER order,
              const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
              const double alpha, const double *A, const int lda,
              const double *X, const int incX, const double beta,
              double *Y, const int incY);
void my_dger(const enum CBLAS_ORDER order, const int M, const int N,
             const double alpha, const double *X, const int incX,
             const double *Y, const int incY, double *A, const int lda);
