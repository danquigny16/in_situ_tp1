#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util.h"


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
  printf("----- Test de resultat -----\n\n");

  // Initialisation des vecteurs
  double * vec1 = vecteur(5);
  double * vec2 = vecteur(5);

  // Initialisation des deux vecteurs
  init_vecteur(5, 1, vec1);
  init_vecteur(5, 1, vec2);

  printf("Vecteur 1\n\n");
  affiche_vecteur(5, vec1, 1, stdout);
  printf("\nVecteur 2\n\n");
  affiche_vecteur(5, vec1, 1, stdout);
  printf("\nResultat du produit scalaire (attendue 55) : %f\n\n", my_ddot(5, vec1, 1, vec2, 1));

  // Libération mémoire des précédents vecteurs
  free_vecteur(vec1);
  free_vecteur(vec2);

  //////////////////////////////////////////////////////////////////////////////
  // Test du produit scalaire : temps d'éxecution
  printf("----- Test de temps d'execution -----\n\n");

  // Initialisation des variables
  int m = 50;
  clock_t debut, fin;

  // Tests pour des tailles de vecteurs croissant
  while (m < 1000000){
    // Allocation de deux vecteurs de taille m
    vec1 = vecteur(m);
    vec2 = vecteur(m);

    // Produit scalaire
    debut = clock();
    my_ddot(m, vec1, 1, vec2, 1);
    fin = clock();

    // Affichage des performances
    double temps = ((double) (fin - debut)) / ((double) CLOCKS_PER_SEC);
    double flop = (double) (2 * m - 1); // m multiplication et m-1 addition
    double Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des vecteurs de taille %7d : %10.6f Mflop/s\n", m, Mflop_s);

    // Libération mémoire des précédents vecteurs
    free_vecteur(vec1);
    free_vecteur(vec2);

    // Increment de m de 25%
    int increment = m / 4;
    m += increment;
  }

  // Fin du test
  printf("\n");
}

/**
Fait les tests relatifs au produit de matrices my_dgemm_scalaire()
*/
void test_my_dgemm_scalaire(){
  // Nom du test
  printf("******************************************************************************\n");
  printf("*************** TEST PRODUIT DE MATRICES : MY_DGEMM_SCALAIRE() ***************\n");
  printf("******************************************************************************\n\n");

  //////////////////////////////////////////////////////////////////////////////
  // Test du produit de matrices : le résultat est-il correct
  printf("----- Test de resultat -----\n\n");

  // Allocation mémoire des matrices
  double * A = matrice(5, 5);
  double * B = matrice(5, 5);
  double * C = matrice(5, 5);

  // Initialisation des matrice A et B
  init_matrice(5, 5, 5, A);
  init_matrice(5, 5, 5, B);

  // Affichage des matrices
  printf("----- Matrice A -----\n\n");
  affiche(5, 5, A, 5, stdout);

  printf("\n----- Matrice B -----\n\n");
  affiche(5, 5, B, 5, stdout);

  printf("\n----- Resultat attendue -----\n\n");
  printf("855.000000  910.000000  965.000000  1020.000000  1075.000000\n");
  printf("910.000000  970.000000  1030.000000  1090.000000  1150.000000\n");
  printf("965.000000  1030.000000  1095.000000  1160.000000  1225.000000\n");
  printf("1020.000000  1090.000000  1160.000000  1230.000000  1300.000000\n");
  printf("1075.000000  1150.000000  1225.000000  1300.000000  1375.000000\n");

  printf("\n----- Matrice (scalaire) C = A * B -----\n\n");
  my_dgemm_scalaire(CblasColMajor, CblasTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  affiche(5, 5, C, 5, stdout);
  matrice_zero(5, 5, 5, C);

  printf("\n----- Matrice (ordre kij) C = A * B -----\n\n");
  my_dgemm_scalaire_kij(CblasColMajor, CblasTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  affiche(5, 5, C, 5, stdout);
  matrice_zero(5, 5, 5, C);

  printf("\n----- Matrice (ordre ijk) C = A * B -----\n\n");
  my_dgemm_scalaire_ijk(CblasColMajor, CblasTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  affiche(5, 5, C, 5, stdout);
  matrice_zero(5, 5, 5, C);

  printf("\n----- Matrice (ordre jik) C = A * B -----\n\n");
  my_dgemm_scalaire_jik(CblasColMajor, CblasTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  affiche(5, 5, C, 5, stdout);
  matrice_zero(5, 5, 5, C);

  printf("\n----- Matrice (par bloc) C = A * B -----\n\n");
  my_dgemm_scalaire(CblasColMajor, CblasTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);
  affiche(5, 5, C, 5, stdout);
  matrice_zero(5, 5, 5, C);

  // Libération mémoire des matrices
  free_matrice(A);
  free_matrice(B);
  free_matrice(C);

  //////////////////////////////////////////////////////////////////////////////
  // Test du produit scalaire : temps d'éxecution
  printf("\n----- Test de temps d'execution -----\n");

  // Variable pour mesure du temps
  clock_t debut, fin;

  //////////////////////////////////////////////////////////////////////////////
  // Tests pour des tailles de matrices croissantes avec l'ordre kij
  printf("\nTest dgemm avec l'ordre kij :\n\n");
  for (int size_mat = 100; size_mat <= 1000; size_mat += 100){
    // Allocation mémoire des matrices
    A = matrice(size_mat, size_mat);
    B = matrice(size_mat, size_mat);
    C = matrice(size_mat, size_mat);
    matrice_zero(size_mat, size_mat, size_mat, C);

    // Produit de matrices
    debut = clock();
    my_dgemm_scalaire_kij(CblasColMajor, CblasTrans, CblasNoTrans, size_mat, size_mat, size_mat, 1, A, 5, B, 5, 0, C, 5);
    fin = clock();

    // Affichage des performances
    double temps = ((double) (fin - debut)) / ((double) CLOCKS_PER_SEC);
    // (size_mat ligne) * (size_mat colonne) * (size_mat multiplication et size_mat - 1 addition)
    double flop = (double) (size_mat * size_mat * (2 * size_mat - 1));
    double Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des matrice de taille %7d^2 : %10.6f Mflop/s pour un temps de %f s\n", size_mat, Mflop_s, temps);

    // Libération mémoire des matrices
    free_matrice(A);
    free_matrice(B);
    free_matrice(C);
  }

  // Tests pour des tailles de matrices croissantes avec l'ordre ijk
  printf("\nTest dgemm avec l'ordre ijk :\n\n");
  for (int size_mat = 100; size_mat <= 1000; size_mat += 100){
    // Allocation mémoire des matrices
    A = matrice(size_mat, size_mat);
    B = matrice(size_mat, size_mat);
    C = matrice(size_mat, size_mat);
    matrice_zero(size_mat, size_mat, size_mat, C);

    // Produit de matrices
    debut = clock();
    my_dgemm_scalaire_ijk(CblasColMajor, CblasTrans, CblasNoTrans, size_mat, size_mat, size_mat, 1, A, 5, B, 5, 0, C, 5);
    fin = clock();

    // Affichage des performances
    double temps = ((double) (fin - debut)) / ((double) CLOCKS_PER_SEC);
    // (size_mat ligne) * (size_mat colonne) * (size_mat multiplication et size_mat - 1 addition)
    double flop = (double) (size_mat * size_mat * (2 * size_mat - 1));
    double Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des matrice de taille %7d^2 : %10.6f Mflop/s pour un temps de %f s\n", size_mat, Mflop_s, temps);

    // Libération mémoire des matrices
    free_matrice(A);
    free_matrice(B);
    free_matrice(C);
  }

  // Tests pour des tailles de matrices croissantes avec l'ordre jik
  printf("\nTest dgemm avec l'ordre jik :\n\n");
  for (int size_mat = 100; size_mat <= 1000; size_mat += 100){
    // Allocation mémoire des matrices
    A = matrice(size_mat, size_mat);
    B = matrice(size_mat, size_mat);
    C = matrice(size_mat, size_mat);
    matrice_zero(size_mat, size_mat, size_mat, C);

    // Produit de matrices
    debut = clock();
    my_dgemm_scalaire_jik(CblasColMajor, CblasTrans, CblasNoTrans, size_mat, size_mat, size_mat, 1, A, 5, B, 5, 0, C, 5);
    fin = clock();

    // Affichage des performances
    double temps = ((double) (fin - debut)) / ((double) CLOCKS_PER_SEC);
    // (size_mat ligne) * (size_mat colonne) * (size_mat multiplication et size_mat - 1 addition)
    double flop = (double) (size_mat * size_mat * (2 * size_mat - 1));
    double Mflop_s = (flop / temps) / 1000000.0;
    printf("Performance obtenu pour des matrice de taille %7d^2 : %10.6f Mflop/s pour un temps de %f s\n", size_mat, Mflop_s, temps);

    // Libération mémoire des matrices
    free_matrice(A);
    free_matrice(B);
    free_matrice(C);
  }

  // Tests pour des tailles de matrices croissantes
  // printf("\nTest dgemm par bloc :\n\n");
  // for (int size_mat = 100; size_mat <= 1000; size_mat += 100){
  //   // Allocation mémoire des matrices
  //   A = matrice(size_mat, size_mat);
  //   B = matrice(size_mat, size_mat);
  //   C = matrice(size_mat, size_mat);
  //   matrice_zero(size_mat, size_mat, size_mat, C);
  //
  //   // Produit de matrices
  //   debut = clock();
  //   my_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, size_mat, size_mat, size_mat, 1, A, 5, B, 5, 0, C, 5);
  //   fin = clock();
  //
  //   // Affichage des performances
  //   double temps = ((double) (fin - debut)) / ((double) CLOCKS_PER_SEC);
  //   // (size_mat ligne) * (size_mat colonne) * (size_mat multiplication et size_mat - 1 addition)
  //   double flop = (double) (size_mat * size_mat * (2 * size_mat - 1));
  //   double Mflop_s = (flop / temps) / 1000000.0;
  //   printf("Performance obtenu pour des matrice de taille %7d^2 : %10.6f Mflop/s pour un temps de %f s\n", size_mat, Mflop_s, temps);
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
// Main

int main(/*int argc, char ** argv*/){
  test_alloc_et_free();
  test_initialisation();
  test_my_ddot();
  test_my_dgemm_scalaire();

  return 0;
}
