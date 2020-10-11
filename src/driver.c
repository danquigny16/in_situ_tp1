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
  printf("*************** TEST ALLOC ET FREE ***************\n");

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
  printf("*************** TEST INITIALISATION ***************\n");

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
  printf("*************** TEST PRODUIT SCALAIRE : MY_DDOT() ***************\n\n");

  // Initialisation des variables
  int m = 50;
  clock_t debut, fin;

  while (m < 1000000){
    // Déclaration de deux vecteurs de taille m
    double * vec1 = vecteur(m);
    double * vec2 = vecteur(m);

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
  printf("*************** TEST PRODUIT DE MATRICES : MY_DGEMM_SCALAIRE() ***************\n\n");

  // Allocation matrice
  double * A = matrice(5, 5);
  double * B = matrice(5, 5);
  double * C = matrice(5, 5);

  // Initialisation des matrice A et B
  init_matrice(5, 5, 5, A);
  init_matrice(5, 5, 5, B);

  // Test de dgemm pour la calcul de tA * B
  my_dgemm_scalaire(CblasColMajor, CblasTrans, CblasNoTrans, 5, 5, 5, 1, A, 5, B, 5, 0, C, 5);

  // Affichage des matrices
  printf("----- Matrice A -----\n\n");
  affiche(5, 5, A, 5, stdout);
  printf("\n----- Matrice B -----\n\n");
  affiche(5, 5, B, 5, stdout);
  printf("\n----- Matrice C -----\n\n");
  affiche(5, 5, C, 5, stdout);

  // Libération mémoire
  free_matrice(A);
  free_matrice(B);
  free_matrice(C);

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
