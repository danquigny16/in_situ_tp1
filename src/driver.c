#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "util.h"
#include "myblas.h"
#include "mylapack.h"


////////////////////////////////////////////////////////////////////////////////
// Tests

/**
Fait les tests relatifs aux allocation et libération mémoire de matrices et de vecteurs
*/
void test_alloc_et_free(){
  // Nom du test
  printf("**************************************************\n");
  printf("*************** TEST ALLOC ET FREE ***************\n");
  printf("**************************************************\n");

  // Matrice
  printf("\n----- Matrice -----\n\n");

  double * a = matrice(3, 4);
  affiche(3, 4, a, 3, stdout);
  free_matrice(a);

  // Vecteur
  printf("\n----- Vecteur -----\n\n");

  double * v = vecteur(5);
  affiche_vecteur(5, v, 1, stdout);
  free_vecteur(v);

  // Fin du test
  printf("\n");
}


////////////////////////////////////////////////////////////////////////////////

/**
Fait les tests relatifs aux initialisations de matrices et de vecteurs
*/
void test_initialisation(){
  // Nom du test
  printf("***************************************************\n");
  printf("*************** TEST INITIALISATION ***************\n");
  printf("***************************************************\n");

  // Matrice
  printf("\n----- Matrice -----\n\n");

  double * a = matrice(3, 4);
  init_matrice(3, 4, 3, a);
  affiche(3, 4, a, 3, stdout);

  printf("\n----- On met 42.0 à la position (1,3) -----\n\n");

  fixer_elt_matrice(1, 3, 3, a, 42.0);
  affiche(3, 4, a, 3, stdout);

  free_matrice(a);

  // Vecteur
  printf("\n----- Vecteur -----\n\n");

  double * v = vecteur(5);
  init_vecteur(5, 1, v);
  affiche_vecteur(5, v, 1, stdout);

  printf("\n----- On met 42.0 à la position 2 -----\n\n");

  fixer_elt_vecteur(2, 1, v, 42.0);
  affiche_vecteur(5, v, 1, stdout);

  free_vecteur(v);

  // Fin du test
  printf("\n");
}


////////////////////////////////////////////////////////////////////////////////

/**
Fait les tests relatifs au produit scalaire my_ddot()
*/
void test_my_ddot(){
  // Nom du test
  printf("*******************************************************************\n");
  printf("*************** TEST PRODUIT SCALAIRE : MY_DDOT() *****************\n");
  printf("*******************************************************************\n\n");

  //////////////////////////////////////////////////////////////////////////////
  // Test du produit scalaire : le résultat est-il correct
  printf("----- Test de resultat (vecteur) -----\n\n");

  // Initialisation des vecteurs
  double * vec1 = vecteur(5);
  double * vec2 = vecteur(5);

  // Initialisation des deux vecteurs
  init_vecteur(5, 1, vec1);
  init_vecteur(5, 1, vec2);

  printf("Vecteur 1\n\n");
  affiche_vecteur(5, vec1, 1, stdout);
  printf("\nVecteur 2\n\n");
  affiche_vecteur(5, vec2, 1, stdout);
  printf("\nResultat du produit scalaire (attendue 55) : %f\n", my_ddot(5, vec1, 1, vec2, 1));
  printf("\nResultat du produit scalaire unroll (attendue 55) : %f\n", my_ddot_unroll(5, vec1, 1, vec2, 1));
  printf("\nResultat du produit scalaire openmp (attendue 55) : %f\n", my_ddot_openmp(5, vec1, 1, vec2, 1));
  printf("\nResultat du produit scalaire openmp simd (attendue 55) : %f\n", my_ddot_openmp_simd(5, vec1, 1, vec2, 1));
  printf("\nResultat du produit scalaire avx2 (attendue 55) : %f\n", my_ddot_avx2(5, vec1, 1, vec2, 1));
  printf("\nResultat du produit scalaire avx2-fma (attendue 55) : %f\n", my_ddot_avx2_fma(5, vec1, 1, vec2, 1));
  printf("\nResultat du produit scalaire avx2-fma openmp (attendue 55) : %f\n\n", my_ddot_avx2_fma_openmp(5, vec1, 1, vec2, 1));

  // Libération mémoire des précédents vecteurs
  free_vecteur(vec1);
  free_vecteur(vec2);

  //////////////////////////////////////////////////////////////////////////////

  printf("\n----- Test de resultat (matrice) -----\n\n");

  // Initialisation des matrice
  double * mat1 = matrice(8, 5);
  double * mat2 = matrice(8, 5);

  // Initialisation des deux vecteurs
  init_matrice(8, 5, 8, mat1);
  init_matrice(8, 5, 8, mat2);

  printf("Matrice 1\n\n");
  affiche(8, 5, mat1, 8, stdout);
  printf("\nMatrice 2\n\n");
  affiche(8, 5, mat2, 8, stdout);
  printf("\nResultat du produit scalaire (attendue 55) : %f\n", my_ddot(5, mat1, 8, mat2, 8));
  printf("\nResultat du produit scalaire unroll (attendue 55) : %f\n", my_ddot_unroll(5, mat1, 8, mat2, 8));
  printf("\nResultat du produit scalaire openmp (attendue 55) : %f\n", my_ddot_openmp(5, mat1, 8, mat2, 8));
  printf("\nResultat du produit scalaire openmp simd (attendue 55) : %f\n", my_ddot_openmp_simd(5, mat1, 8, mat2, 8));
  printf("\nResultat du produit scalaire avx2 (attendue 55) : %f\n", my_ddot_avx2(5, mat1, 8, mat2, 8));
  printf("\nResultat du produit scalaire avx2-fma (attendue 55) : %f\n", my_ddot_avx2_fma(5, mat1, 8, mat2, 8));
  printf("\nResultat du produit scalaire avx2-fma openmp (attendue 55) : %f\n\n", my_ddot_avx2_fma_openmp(5, mat1, 8, mat2, 8));

  // Libération mémoire des précédents vecteurs
  free_matrice(mat1);
  free_matrice(mat2);

  //////////////////////////////////////////////////////////////////////////////
  // Test du produit scalaire : temps d'éxecution
  printf("-----------------------------------------------\n");
  printf("----- Test de temps d'execution (vecteur) -----\n");
  printf("-----------------------------------------------\n\n");

  // Initialisation des variables
  int m = 50;
  double debut, fin;

  // Tests pour des tailles de vecteurs croissant
  while (m < 1000000){
    // Allocation de deux vecteurs de taille m
    vec1 = vecteur(m);
    vec2 = vecteur(m);

    ////////////////////////////////////////////////////////////////////////////

    // Produit scalaire
    debut = omp_get_wtime();
    my_ddot(m, vec1, 1, vec2, 1);
    fin = omp_get_wtime();

    // Affichage des performances
    double temps0 = fin - debut;
    double flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    double Mflop_s = (flop / temps0) / 1000000.0;
    printf("Performance obtenu pour des vecteurs de taille %7d :                         ", m);
    printf("%12.6f Mflop/s pour un temps de %f\n", Mflop_s, temps0);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_vecteur(vec1);
    free_vecteur(vec2);
    vec1 = vecteur(m);
    vec2 = vecteur(m);

    // Produit scalaire unroll
    debut = omp_get_wtime();
    my_ddot_unroll(m, vec1, 1, vec2, 1);
    fin = omp_get_wtime();

    // Affichage des performances
    double temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs de taille %7d avec unroll :             ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_vecteur(vec1);
    free_vecteur(vec2);
    vec1 = vecteur(m);
    vec2 = vecteur(m);

    // Produit scalaire openmp
    debut = omp_get_wtime();
    my_ddot_openmp(m, vec1, 1, vec2, 1);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs de taille %7d avec openmp :             ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_vecteur(vec1);
    free_vecteur(vec2);
    vec1 = vecteur(m);
    vec2 = vecteur(m);

    // Produit scalaire openmp simd
    debut = omp_get_wtime();
    my_ddot_openmp_simd(m, vec1, 1, vec2, 1);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs de taille %7d avec openmp simd :        ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_vecteur(vec1);
    free_vecteur(vec2);
    vec1 = vecteur(m);
    vec2 = vecteur(m);

    // Produit scalaire avx2
    debut = omp_get_wtime();
    my_ddot_avx2(m, vec1, 1, vec2, 1);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs de taille %7d avec avx2 :               ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_vecteur(vec1);
    free_vecteur(vec2);
    vec1 = vecteur(m);
    vec2 = vecteur(m);

    // Produit scalaire avx2-fma
    debut = omp_get_wtime();
    my_ddot_avx2_fma(m, vec1, 1, vec2, 1);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs de taille %7d avec avx2-fma :           ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_vecteur(vec1);
    free_vecteur(vec2);
    vec1 = vecteur(m);
    vec2 = vecteur(m);

    // Produit scalaire avx2-fma et openmp
    debut = omp_get_wtime();
    my_ddot_avx2_fma_openmp(m, vec1, 1, vec2, 1);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs de taille %7d avec avx2-fma et openmp : ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // Libération mémoire des précédents vecteurs
    free_vecteur(vec1);
    free_vecteur(vec2);

    // Increment de m de 25%
    int increment = m / 4;
    m += increment;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Test du produit scalaire : temps d'éxecution
  printf("-----------------------------------------------\n");
  printf("----- Test de temps d'execution (matrice) -----\n");
  printf("-----------------------------------------------\n\n");

  // Initialisation des variables
  m = 50;

  // Tests pour des tailles de vecteurs croissant
  while (m < 1000000){
    // Allocation de deux vecteurs de taille m
    mat1 = matrice(10, m);
    mat2 = matrice(10, m);

    ////////////////////////////////////////////////////////////////////////////

    // Produit scalaire
    debut = omp_get_wtime();
    my_ddot(m, mat1, 10, mat2, 10);
    fin = omp_get_wtime();

    // Affichage des performances
    double temps0 = fin - debut;
    double flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    double Mflop_s = (flop / temps0) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d :                         ", m);
    printf("%12.6f Mflop/s pour un temps de %f\n", Mflop_s, temps0);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(10, m);
    mat2 = matrice(10, m);

    // Produit scalaire unroll
    debut = omp_get_wtime();
    my_ddot_unroll(m, mat1, 10, mat2, 10);
    fin = omp_get_wtime();

    // Affichage des performances
    double temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec unroll :             ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(10, m);
    mat2 = matrice(10, m);

    // Produit scalaire openmp
    debut = omp_get_wtime();
    my_ddot_openmp(m, mat1, 10, mat2, 10);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec openmp :             ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(10, m);
    mat2 = matrice(10, m);

    // Produit scalaire openmp simd
    debut = omp_get_wtime();
    my_ddot_openmp_simd(m, mat1, 10, mat2, 10);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec openmp simd :        ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(10, m);
    mat2 = matrice(10, m);

    // Produit scalaire avx2
    debut = omp_get_wtime();
    my_ddot_avx2(m, mat1, 10, mat2, 10);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec avx2 :               ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(10, m);
    mat2 = matrice(10, m);

    // Produit scalaire avx2-fma
    debut = omp_get_wtime();
    my_ddot_avx2_fma(m, mat1, 10, mat2, 10);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec avx2-fma :           ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(10, m);
    mat2 = matrice(10, m);

    // Produit scalaire avx2-fma et openmp
    debut = omp_get_wtime();
    my_ddot_avx2_fma_openmp(m, mat1, 10, mat2, 10);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec avx2-fma et openmp : ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // Libération mémoire des précédents vecteurs
    free_matrice(mat1);
    free_matrice(mat2);

    // Increment de m de 25%
    int increment = m / 4;
    m += increment;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Test du produit scalaire : temps d'éxecution
  printf("----------------------------------------------------\n");
  printf("----- Test de temps d'execution (matrice hard) -----\n");
  printf("--------------------------------------------------\n\n");

  // Initialisation des variables
  m = 50;

  // Tests pour des tailles de vecteurs croissant
  while (m < 10000){
    // Allocation de deux vecteurs de taille m
    mat1 = matrice(m, m);
    mat2 = matrice(m, m);

    ////////////////////////////////////////////////////////////////////////////

    // Produit scalaire
    debut = omp_get_wtime();
    my_ddot(m, mat1, m, mat2, m);
    fin = omp_get_wtime();

    // Affichage des performances
    double temps0 = fin - debut;
    double flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    double Mflop_s = (flop / temps0) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d :                         ", m);
    printf("%12.6f Mflop/s pour un temps de %f\n", Mflop_s, temps0);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(m, m);
    mat2 = matrice(m, m);

    // Produit scalaire unroll
    debut = omp_get_wtime();
    my_ddot_unroll(m, mat1, m, mat2, m);
    fin = omp_get_wtime();

    // Affichage des performances
    double temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec unroll :             ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(m, m);
    mat2 = matrice(m, m);

    // Produit scalaire openmp
    debut = omp_get_wtime();
    my_ddot_openmp(m, mat1, m, mat2, m);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec openmp :             ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(m, m);
    mat2 = matrice(m, m);

    // Produit scalaire openmp simd
    debut = omp_get_wtime();
    my_ddot_openmp_simd(m, mat1, m, mat2, m);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec openmp simd :        ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(m, m);
    mat2 = matrice(m, m);

    // Produit scalaire avx2
    debut = omp_get_wtime();
    my_ddot_avx2(m, mat1, m, mat2, m);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec avx2 :               ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(m, m);
    mat2 = matrice(m, m);

    // Produit scalaire avx2-fma
    debut = omp_get_wtime();
    my_ddot_avx2_fma(m, mat1, m, mat2, m);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec avx2-fma :           ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // liberation et reallocation
    free_matrice(mat1);
    free_matrice(mat2);
    mat1 = matrice(m, m);
    mat2 = matrice(m, m);

    // Produit scalaire avx2-fma et openmp
    debut = omp_get_wtime();
    my_ddot_avx2_fma_openmp(m, mat1, m, mat2, m);
    fin = omp_get_wtime();

    // Affichage des performances
    temps = fin - debut;
    flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs (1ere ligne matrice) de taille %7d avec avx2-fma et openmp : ", m);
    printf("%12.6f Mflop/s pour un temps de %f", Mflop_s, temps);
    printf(", soit une amelioration de %f\n\n", temps0 / temps);

    ////////////////////////////////////////////////////////////////////////////

    // Libération mémoire des précédents vecteurs
    free_matrice(mat1);
    free_matrice(mat2);

    // Increment de m de 25%
    int increment = m / 4;
    m += increment;
  }

  //////////////////////////////////////////////////////////////////////////////

  // Fin du test
  printf("\n");
}

////////////////////////////////////////////////////////////////////////////////

/**
Fait les tests relatifs au produit de matrices my_dgemm_scalaire()
*/
void test_my_dgemm(){
  // Nom du test
  printf("*********************************************************************\n");
  printf("*************** TEST PRODUIT DE MATRICES : MY_DGEMM() ***************\n");
  printf("*********************************************************************\n\n");

  //////////////////////////////////////////////////////////////////////////////
  // Test du produit de matrices (scalaire) : le résultat est-il correct
  printf("***** Test de resultat (scalaire) *****\n\n");

  // Allocation mémoire des matrices
  double * A = matrice(5, 5);
  double * B = matrice(5, 5);
  double * C = matrice(5, 5);

  // Initialisation des matrice A et B
  // init_matrice(5, 5, 5, A);
  // init_matrice(5, 5, 5, B);
  //
  // // Affichage des matrices
  // printf("----- Matrice A -----\n\n");
  // affiche(5, 5, A, 5, stdout);
  //
  // printf("\n----- Matrice B -----\n\n");
  // affiche(5, 5, B, 5, stdout);
  //
  // printf("\n----- Resultat attendue -----\n\n");
  // printf("[[ 855  910  965 1020 1075]\n");
  // printf(" [ 910  970 1030 1090 1150]\n");
  // printf(" [ 965 1030 1095 1160 1225]\n");
  // printf(" [1020 1090 1160 1230 1300]\n");
  // printf(" [1075 1150 1225 1300 1375]]\n");
  //
  // printf("\n----- Matrice (scalaire) C = A * B -----\n\n");
  // my_dgemm_scalaire(CblasColMajor, CblasTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  // affiche(5, 5, C, 5, stdout);
  // matrice_zero(5, 5, 5, C);
  //
  // printf("\n----- Matrice (ordre kij) C = A * B -----\n\n");
  // my_dgemm_scalaire_kij(CblasColMajor, CblasTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  // affiche(5, 5, C, 5, stdout);
  // matrice_zero(5, 5, 5, C);
  //
  // printf("\n----- Matrice (ordre ijk) C = A * B -----\n\n");
  // my_dgemm_scalaire_ijk(CblasColMajor, CblasTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  // affiche(5, 5, C, 5, stdout);
  // matrice_zero(5, 5, 5, C);
  //
  // printf("\n----- Matrice (ordre jik) C = A * B -----\n\n");
  // my_dgemm_scalaire_jik(CblasColMajor, CblasTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  // affiche(5, 5, C, 5, stdout);
  // matrice_zero(5, 5, 5, C);

  //////////////////////////////////////////////////////////////////////////////

  free_matrice(A);
  free_matrice(B);
  free_matrice(C);

  A = matrice(5, 5);
  B = matrice(5, 5);
  C = matrice(5, 5);

  init_matrice(5, 5, 5, A);
  init_matrice(5, 5, 5, B);
  init_matrice(5, 5, 5, C);

  printf("\n----- Matrice (simd) correct C = A * B -----\n\n");
  dgemm_seq_opti(CblasColMajor, CblasNoTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  affiche(5, 5, C, 5, stdout);
  matrice_zero(5, 5, 5, C);

  printf("\n----- Matrice (simd) C = A * B -----\n\n");
  dgemm_seq_simd(CblasColMajor, CblasNoTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  affiche(5, 5, C, 5, stdout);
  matrice_zero(5, 5, 5, C);

  //////////////////////////////////////////////////////////////////////////////

  // Libération mémoire des matrices
  free_matrice(A);
  free_matrice(B);
  free_matrice(C);

  //////////////////////////////////////////////////////////////////////////////
  // // Test du produit de matrices (scalaire) : le résultat est-il correct
  // printf("\n***** Test de resultat (par blocs) *****\n\n");
  //
  // // Allocation mémoire des matrices
  // A = matrice(10, 10);
  // B = matrice(10, 10);
  // C = matrice(10, 10);
  //
  // // Initialisation des matrice A et B
  // init_matrice(10, 10, 10, A);
  // init_matrice(10, 10, 10, B);
  //
  // // Affichage des matrices
  // printf("----- Matrice A -----\n\n");
  // affiche(10, 10, A, 10, stdout);
  //
  // printf("\n----- Matrice B -----\n\n");
  // affiche(10, 10, B, 10, stdout);
  //
  // printf("\n----- Resultat attendue -----\n\n");
  // printf("[[29410 29870 30330 30790 31250 31710 32170 32630 33090 33550]\n");
  // printf(" [29870 30340 30810 31280 31750 32220 32690 33160 33630 34100]\n");
  // printf(" [30330 30810 31290 31770 32250 32730 33210 33690 34170 34650]\n");
  // printf(" [30790 31280 31770 32260 32750 33240 33730 34220 34710 35200]\n");
  // printf(" [31250 31750 32250 32750 33250 33750 34250 34750 35250 35750]\n");
  // printf(" [31710 32220 32730 33240 33750 34260 34770 35280 35790 36300]\n");
  // printf(" [32170 32690 33210 33730 34250 34770 35290 35810 36330 36850]\n");
  // printf(" [32630 33160 33690 34220 34750 35280 35810 36340 36870 37400]\n");
  // printf(" [33090 33630 34170 34710 35250 35790 36330 36870 37410 37950]\n");
  // printf(" [33550 34100 34650 35200 35750 36300 36850 37400 37950 38500]]\n");
  //
  // printf("\n----- Matrice (par bloc) C = A * B -----\n\n");
  // my_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 10, 10, 10, 1, A, 10, B, 10, 0, C, 10);
  // affiche(10, 10, C, 10, stdout);
  // matrice_zero(10, 10, 10, C);
  //
  // // Libération mémoire des matrices
  // free_matrice(A);
  // free_matrice(B);
  // free_matrice(C);
  //
  // //////////////////////////////////////////////////////////////////////////////
  // // Test du produit scalaire : temps d'éxecution
  // printf("\n***** Test de temps d'execution *****\n");
  //
  // // Variable pour mesure du temps
  // double debut, fin;
  //
  // //////////////////////////////////////////////////////////////////////////////
  // // Tests pour des tailles de matrices croissantes avec l'ordre kij
  // printf("\nTest dgemm avec l'ordre kij :\n\n");
  // for (int size_mat = 1000; size_mat <= 2000; size_mat += 200){
  //   // Allocation mémoire des matrices
  //   A = matrice(size_mat, size_mat);
  //   B = matrice(size_mat, size_mat);
  //   C = matrice(size_mat, size_mat);
  //   matrice_zero(size_mat, size_mat, size_mat, C);
  //
  //   // Produit de matrices
  //   debut = omp_get_wtime();
  //   my_dgemm_scalaire_kij(CblasColMajor, CblasTrans, CblasNoTrans, size_mat, size_mat, size_mat, 1, A, 5, B, 5, 0, C, 5);
  //   fin = omp_get_wtime();
  //
  //   // Affichage des performances
  //   double temps = fin - debut;
  //   // (size_mat ligne) * (size_mat colonne) * (size_mat multiplication et size_mat - 1 addition)
  //   int Mflop = (size_mat/1000 * size_mat/1000 * (2 * size_mat - 1));
  //   double Mflop_s = (Mflop / temps);
  //   printf("Performance obtenu pour des matrices de taille %7d ordre kij : %10.6f Mflop/s pour un temps de %f s\n", size_mat, Mflop_s, temps);
  //
  //   // Libération mémoire des matrices
  //   free_matrice(A);
  //   free_matrice(B);
  //   free_matrice(C);
  // }
  //
  // // Tests pour des tailles de matrices croissantes avec l'ordre ijk
  // printf("\nTest dgemm avec l'ordre ijk :\n\n");
  // for (int size_mat = 1000; size_mat <= 2000; size_mat += 200){
  //   // Allocation mémoire des matrices
  //   A = matrice(size_mat, size_mat);
  //   B = matrice(size_mat, size_mat);
  //   C = matrice(size_mat, size_mat);
  //   matrice_zero(size_mat, size_mat, size_mat, C);
  //
  //   // Produit de matrices
  //   debut = omp_get_wtime();
  //   my_dgemm_scalaire_ijk(CblasColMajor, CblasTrans, CblasNoTrans, size_mat, size_mat, size_mat, 1, A, 5, B, 5, 0, C, 5);
  //   fin = omp_get_wtime();
  //
  //   // Affichage des performances
  //   double temps = fin - debut;
  //   // (size_mat ligne) * (size_mat colonne) * (size_mat multiplication et size_mat - 1 addition)
  //   int Mflop = (size_mat/1000 * size_mat/1000 * (2 * size_mat - 1));
  //   double Mflop_s = (Mflop / temps);
  //
  //   printf("Performance obtenu pour des matrices de taille %7d ordre ijk : %10.6f Mflop/s pour un temps de %f s\n", size_mat, Mflop_s, temps);
  //
  //   // Libération mémoire des matrices
  //   free_matrice(A);
  //   free_matrice(B);
  //   free_matrice(C);
  // }
  //
  // // Tests pour des tailles de matrices croissantes avec l'ordre jik
  // printf("\nTest dgemm avec l'ordre jik :\n\n");
  // for (int size_mat = 1000; size_mat <= 2000; size_mat += 200){
  //   // Allocation mémoire des matrices
  //   A = matrice(size_mat, size_mat);
  //   B = matrice(size_mat, size_mat);
  //   C = matrice(size_mat, size_mat);
  //   matrice_zero(size_mat, size_mat, size_mat, C);
  //
  //   // Produit de matrices
  //   debut = omp_get_wtime();
  //   my_dgemm_scalaire_jik(CblasColMajor, CblasTrans, CblasNoTrans, size_mat, size_mat, size_mat, 1, A, 5, B, 5, 0, C, 5);
  //   fin = omp_get_wtime();
  //
  //   // Affichage des performances
  //   double temps = fin - debut;
  //   // (size_mat ligne) * (size_mat colonne) * (size_mat multiplication et size_mat - 1 addition)
  //   int Mflop = (size_mat/1000 * size_mat/1000 * (2 * size_mat - 1));
  //   double Mflop_s = (Mflop / temps);
  //
  //   printf("Performance obtenu pour des matrices de taille %7d ordre jik : %10.6f Mflop/s pour un temps de %f s\n", size_mat, Mflop_s, temps);
  //
  //   // Libération mémoire des matrices
  //   free_matrice(A);
  //   free_matrice(B);
  //   free_matrice(C);
  // }
  //
  // // Tests pour des tailles de matrices croissantes avec l'ordre kji
  // printf("\nTest dgemm avec l'ordre kji :\n\n");
  // for (int size_mat = 1000; size_mat <= 2000; size_mat += 200){
  //   // Allocation mémoire des matrices
  //   A = matrice(size_mat, size_mat);
  //   B = matrice(size_mat, size_mat);
  //   C = matrice(size_mat, size_mat);
  //   matrice_zero(size_mat, size_mat, size_mat, C);
  //
  //   // Produit de matrices
  //   debut = omp_get_wtime();
  //   my_dgemm_scalaire_kji(CblasColMajor, CblasTrans, CblasNoTrans, size_mat, size_mat, size_mat, 1, A, 5, B, 5, 0, C, 5);
  //   fin = omp_get_wtime();
  //
  //   // Affichage des performances
  //   double temps = fin - debut;
  //   // (size_mat ligne) * (size_mat colonne) * (size_mat multiplication et size_mat - 1 addition)
  //   int Mflop = (size_mat/1000 * size_mat/1000 * (2 * size_mat - 1));
  //   double Mflop_s = (Mflop / temps);
  //   printf("Performance obtenu pour des matrices de taille %7d ordre kji : %10.6f Mflop/s pour un temps de %f s\n", size_mat, Mflop_s, temps);
  //
  //   // Libération mémoire des matrices
  //   free_matrice(A);
  //   free_matrice(B);
  //   free_matrice(C);
  // }
  //
  // // Tests pour des tailles de matrices croissantes avec l'ordre jik
  // printf("\nTest dgemm avec l'ordre kji et boucle for déroulée:\n\n");
  // for (int size_mat = 1000; size_mat <= 2000; size_mat += 200){
  //   // Allocation mémoire des matrices
  //   A = matrice(size_mat, size_mat);
  //   B = matrice(size_mat, size_mat);
  //   C = matrice(size_mat, size_mat);
  //   matrice_zero(size_mat, size_mat, size_mat, C);
  //
  //   // Produit de matrices
  //   debut = omp_get_wtime();
  //   my_dgemm_scalaire_jik_unroll(CblasColMajor, CblasTrans, CblasNoTrans, size_mat, size_mat, size_mat, 1, A, 5, B, 5, 0, C, 5);
  //   fin = omp_get_wtime();
  //
  //   // Affichage des performances
  //   double temps = fin - debut;
  //   // (size_mat ligne) * (size_mat colonne) * (size_mat multiplication et size_mat - 1 addition)
  //   int Mflop = (size_mat/1000 * size_mat/1000 * (2 * size_mat - 1));
  //   double Mflop_s = (Mflop / temps);
  //
  //   printf("Performance obtenu pour des matrices de taille %7d ordre jik unroll: %10.6f Mflop/s pour un temps de %f s\n", size_mat, Mflop_s, temps);
  //
  //   // Libération mémoire des matrices
  //   free_matrice(A);
  //   free_matrice(B);
  //   free_matrice(C);
  // }
  //
  // // Tests pour des tailles de matrices croissantes
  // printf("\nTest dgemm par bloc :\n\n");
  // for (int size_mat = 1000; size_mat <= 2000; size_mat += 200){
  //   // Allocation mémoire des matrices
  //   A = matrice(size_mat, size_mat);
  //   B = matrice(size_mat, size_mat);
  //   C = matrice(size_mat, size_mat);
  //   matrice_zero(size_mat, size_mat, size_mat, C);
  //
  //   // Produit de matrices
  //   debut = omp_get_wtime();
  //   my_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, size_mat, size_mat, size_mat, 1, A, 5, B, 5, 0, C, 5);
  //   fin = omp_get_wtime();
  //
  //   // Affichage des performances
  //   double temps = fin - debut;
  //   // (size_mat ligne) * (size_mat colonne) * (size_mat multiplication et size_mat - 1 addition)
  //   int Mflop = (size_mat/1000 * size_mat/1000 * (2 * size_mat - 1));
  //   double Mflop_s = (Mflop / temps);
  //
  //   printf("Performance obtenu pour des matrices de taille %7d par bloc : %10.6f Mflop/s pour un temps de %f s\n", size_mat, Mflop_s, temps);
  //
  //   // Libération mémoire des matrices
  //   free_matrice(A);
  //   free_matrice(B);
  //   free_matrice(C);
  // }

  // Fin du test
  printf("\n");
}

////////////////////////////////////////////////////////////////////////////////

/**
Fait les tests relatifs aux fonctions blas complémentaires
*/
void test_blas(){
  // Nom du test
  printf("*****************************************************************************************\n");
  printf("*************** TEST FONCTIONS BLAS : MY_DAXPY(), MY_DGEMV(), MY_DGER() *****************\n");
  printf("*****************************************************************************************\n\n");

  //////////////////////////////////////////////////////////////////////////////
  // Test de my_daxpy()
  printf("***** Test de my_daxpy() *****\n\n");

  // Initialisation des vecteurs
  double * vec1 = vecteur(5);
  double * vec2 = vecteur(5);

  // Initialisation des deux vecteurs
  init_vecteur(5, 1, vec1);
  init_vecteur(5, 1, vec2);

  // Affichage des résultats
  printf("----- Vecteur 1 -----\n\n");
  affiche_vecteur(5, vec1, 1, stdout);

  printf("\n----- Vecteur 2 -----\n\n");
  affiche_vecteur(5, vec2, 1, stdout);

  printf("\n----- Resultat attendue -----\n\n");
  printf("[ 3  6  9 12 15]\n");

  my_daxpy(5, 2, vec1, 1, vec2, 1);
  printf("\n----- Resultat obtenue -----\n\n");
  affiche_vecteur(5, vec2, 1, stdout);

  // Libération mémoire des précédents vecteurs
  free_vecteur(vec1);
  free_vecteur(vec2);

  //////////////////////////////////////////////////////////////////////////////
  // Test de my_daxpy()
  printf("\n***** Test de my_dgemv() *****\n\n");

  // Initialisation des vecteurs
  vec1 = vecteur(5);
  vec2 = vecteur(5);
  double * mat = matrice(5, 5);

  // Initialisation des deux vecteurs
  init_vecteur(5, 1, vec1);
  init_vecteur(5, 1, vec2);
  init_matrice(5, 5, 5, mat);

    // Affichage des résultats
  printf("----- Vecteur 1 -----\n\n");
  affiche_vecteur(5, vec1, 1, stdout);

  printf("\n----- Vecteur 2 -----\n\n");
  affiche_vecteur(5, vec2, 1, stdout);

  printf("\n----- Matrice -----\n\n");
  affiche(5, 5, mat, 5, stdout);

  printf("\n----- Resultat attendue -----\n\n");
  printf("[ 56 132 208 284 360]\n");

  my_dgemv(CblasColMajor, CblasNoTrans, 5, 5, 1, mat, 5, vec1, 1, 1, vec2, 1);
  printf("\n----- Resultat obtenue -----\n\n");
  affiche_vecteur(5, vec2, 1, stdout);

  // Libération mémoire des précédents vecteurs
  free_vecteur(vec1);
  free_vecteur(vec2);
  free_matrice(mat);

  //////////////////////////////////////////////////////////////////////////////
  // Test de my_daxpy()
  printf("\n***** Test de my_dger() *****\n\n");

  // Initialisation des vecteurs
  vec1 = vecteur(5);
  vec2 = vecteur(5);
  mat = matrice(5, 5);

  // Initialisation des deux vecteurs
  init_vecteur(5, 1, vec1);
  init_vecteur(5, 1, vec2);
  init_matrice(5, 5, 5, mat);

    // Affichage des résultats
  printf("----- Vecteur 1 -----\n\n");
  affiche_vecteur(5, vec1, 1, stdout);

  printf("\n----- Vecteur 2 -----\n\n");
  affiche_vecteur(5, vec2, 1, stdout);

  printf("\n----- Matrice -----\n\n");
  affiche(5, 5, mat, 5, stdout);

  printf("\n----- Resultat attendue -----\n\n");
  printf("[[ 3  6  9 12 15]\n");
  printf(" [10 15 20 25 30]\n");
  printf(" [17 24 31 38 45]\n");
  printf(" [24 33 42 51 60]\n");
  printf(" [31 42 53 64 75]]\n");

  my_dger(CblasColMajor, 5, 5, 2, vec1, 1, vec2, 1, mat, 5);
  printf("\n----- Resultat obtenue -----\n\n");
  affiche(5, 5, mat, 5, stdout);

  // Libération mémoire des précédents vecteurs
  free_vecteur(vec1);
  free_vecteur(vec2);
  free_matrice(mat);

  //////////////////////////////////////////////////////////////////////////////
  // Fin du test
  printf("\n");
}


////////////////////////////////////////////////////////////////////////////////

/**
Fait les tests relatifs à la factorisation LU
*/
void test_factorisation_LU(){
  // // Nom du test
  // printf("*****************************************************\n");
  // printf("*************** TEST FACTORISATION LU ***************\n");
  // printf("*****************************************************\n");
  //
  // printf("\n\n***** test de resultat en seq *****\n\n");
  //
  // // Initialisation des vecteurs
  // double * vec = vecteur(5);
  // double * mat = matrice(5, 5);
  //
  // // Initialisation des deux vecteurs
  // init_vecteur(5, 1, vec);
  // // Notre init_matrice donne des divisions par zéro, on l'évite
  // init_2_matrice(5, 5, 5, mat);
  //
  // // Affichage des résultats
  // printf("\n----- Matrice -----\n\n");
  // affiche(5, 5, mat, 5, stdout);
  //
  // printf("\n----- Vecteur -----\n\n");
  // affiche_vecteur(5, vec, 1, stdout);
  //
  // printf("\n----- Resultat attendue -----\n\n");
  // printf("[[0.0444439 ]\n");
  // printf(" [0.01093364]\n");
  // printf(" [0.00835592]\n");
  // printf(" [0.00740624]\n");
  // printf(" [0.0069124 ]]\n");
  //
  // my_dgesv(CblasColMajor, 5, 1, mat, 5, NULL, vec, 1);
  // printf("\n----- Resultat obtenue -----\n\n");
  // affiche_vecteur(5, vec, 1, stdout);
  //
  // // Libération mémoire des précédents vecteurs
  // free_vecteur(vec);
  // free_matrice(mat);
  //
  // //////////////////////////////////////////////////////////////////////////////
  // printf("\n\n***** test de resultat par block *****\n\n");
  //
  // int size = 11;
  //
  // // Initialisation des vecteurs
  // vec = vecteur(size);
  // mat = matrice(size, size);
  //
  // // Initialisation des deux vecteurs
  // init_vecteur(size, 1, vec);
  // // Notre init_matrice donne des divisions par zéro, on l'évite
  // init_2_matrice(size, size, size, mat);
  //
  // // Affichage des résultats
  // printf("\n----- Matrice -----\n\n");
  // affiche(size, size, mat, size, stdout);
  //
  // printf("\n----- Vecteur -----\n\n");
  // affiche_vecteur(size, vec, 1, stdout);
  //
  // printf("\n----- Resultat attendue -----\n\n");
  // // printf("0.043652\n");
  // // printf("0.005048\n");
  // // printf("0.003370\n");
  // // printf("0.002778\n");
  // // printf("0.002475\n");
  // // printf("0.002291\n");
  // // printf("0.002168\n");
  // // printf("0.002079\n");
  // // printf("0.002012\n");
  // // printf("0.001960\n");
  //
  // printf("0.043843\n");
  // printf("0.004487\n");
  // printf("0.002913\n");
  // printf("0.002359\n");
  // printf("0.002077\n");
  // printf("0.001906\n");
  // printf("0.001791\n");
  // printf("0.001709\n");
  // printf("0.001647\n");
  // printf("0.001598\n");
  // printf("0.001559\n");
  //
  // my_dgesv(CblasColMajor, size, 1, mat, size, NULL, vec, 1);
  // printf("\n----- Resultat obtenue -----\n\n");
  // affiche_vecteur(size, vec, 1, stdout);
  //
  // // Libération mémoire des précédents vecteurs
  // free_vecteur(vec);
  // free_matrice(mat);

  //////////////////////////////////////////////////////////////////////////////

  printf("************* Test rectangulaire qui pose pb ******************\n");

  double *mat = matrice(57, 3);

  // init_2_matrice(57, 3, 57, mat);
  // my_dgetf2_1(CblasColMajor, 57, 3, mat, 57, NULL);
  // printf("\n----------- matrice factoriser -----------------\n\n");
  // affiche(57, 3, mat, 57, stdout);

  init_2_matrice(57, 3, 57, mat);
  dgetrf_seq_opti(CblasColMajor, 57, 3, mat, 57, NULL);
  printf("\n----------- matrice factoriser -----------------\n\n");
  affiche(57, 3, mat, 57, stdout);

  free_matrice(mat);

  //////////////////////////////////////////////////////////////////////////////

  // printf("\n\n***** test de performance en seq *****\n\n");
  //
  // // Initialisation des variables
  // int m = 100;
  // double debut, fin;
  // double temps;
  //
  // while (m <= 1200){
  //   // Allocation de deux vecteurs de taille m
  //   mat = matrice(m, m);
  //   init_2_matrice(m, m, m, mat);
  //
  //   ////////////////////////////////////////////////////////////////////////////
  //
  //   // dgetrf seq
  //   debut = omp_get_wtime();
  //   my_dgetf2(CblasColMajor, m, m, mat, m, NULL);
  //   fin = omp_get_wtime();
  //
  //   // Affichage des performances
  //   temps = fin - debut;
  //   printf("(dgetf2)        Performance obtenu pour une matrice de taille %7d : %12.6f s\n", m, temps);
  //
  //   ////////////////////////////////////////////////////////////////////////////
  //
  //   // dgetrf par block
  //   debut = omp_get_wtime();
  //   my_dgetrf(CblasColMajor, m, m, mat, m, NULL);
  //   fin = omp_get_wtime();
  //
  //   // Affichage des performances
  //   temps = fin - debut;
  //   printf("(dgetrf)        Performance obtenu pour une matrice de taille %7d : %12.6f s\n", m, temps);
  //
  //   ////////////////////////////////////////////////////////////////////////////
  //
  //   // dgetrf par block avec openmp
  //   debut = omp_get_wtime();
  //   my_dgetrf_openmp(CblasColMajor, m, m, mat, m, NULL);
  //   fin = omp_get_wtime();
  //
  //   // Affichage des performances
  //   temps = fin - debut;
  //   printf("(dgetrf openmp) Performance obtenu pour une matrice de taille %7d : %12.6f s\n", m, temps);
  //
  //   ////////////////////////////////////////////////////////////////////////////
  //
  //   printf("\n");
  //
  //   // Libération mémoire des précédents vecteurs
  //   free_matrice(mat);
  //
  //   // Increment de m de 25%
  //   m += 100;
  // }

  //////////////////////////////////////////////////////////////////////////////

  // Fin du test
  printf("\n");
}


////////////////////////////////////////////////////////////////////////////////
// Main

int main(/*int argc, char ** argv*/){
  // test_alloc_et_free();
  // test_initialisation();
  // test_my_ddot();
  test_my_dgemm();
  // test_blas();
  // test_factorisation_LU();

  return 0;
}
