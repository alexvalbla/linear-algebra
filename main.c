#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "linalg.h"




int main(){


  matrix_t A;
  A.numC = A.numR = 4;
  //vector_t b;
  //b.len = 4;
  double array1[] = {2,3,1,5,6,13,5,19,2,19,10,23,4,10,11,31};
  //double array2[] = {1,3,3,4};
  A.val = array1;
  //b.val = array2;

  matrix_t A_inv = invert_matrix(A);
  matrix_t Id = matrix_product(A, A_inv);
  print_matrix(Id);
  free_matrix(A_inv);
  free_matrix(Id);


  return 0;
}
