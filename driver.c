#include <stdio.h>
#include <stdlib.h>

#include "util.h"

// gcc -Wall -Wextra -std=c99 util.h util.c driver.c -o test
int main (int argc, char ** argv){
  int * a = matrice(3, 4);
  affiche(3, 4, a, 3, stdout);
  free(a);

  return 0;
}
