////////////////////////////////////////////////////////////////////////////////

void my_dgetf2(const enum CBLAS_ORDER Order, int m, int n, double* a, int lda);
void my_cblas_dtrsm(const enum CBLAS_ORDER Order, const int M, const int N, const double *A, const int lda, double *B, const int ldb);
void my_cblas_dgesv(const enum CBLAS_ORDER Order, const int N, double *A, const int lda, double *B, const int ldb);
