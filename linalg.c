#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "linalg.h"

#define LIN(m, n, i, j) ( (i) * (n) + (j) )
//row i column j on m x n matrix_t



vector_t get_matrix_row(matrix_t A, unsigned int p){
  /*
  return row vector p of A
  */
  if(A.val == NULL){
    printf("NULL matrix_t values pointer (get_matrix_row).\n");
    exit(EXIT_FAILURE);
  }
  unsigned int m = A.numR;
  unsigned int n = A.numC;
  if(p >= m){
    printf("Cannot extract line %u, %u lines in matrix_t (get_matrix_row).\n", p, m);
    exit(EXIT_FAILURE);
  }
  vector_t v = allocate_vector(n);
  memcpy(v.val, &A.val[p * n], sizeof(double) * n);
  return v;
}

vector_t get_matrix_column(matrix_t A, unsigned int p){
  /*
  return column vector p of A
  */
  if(A.val == NULL){
    printf("NULL matrix_t values pointer (get_matrix_column).\n");
    exit(EXIT_FAILURE);
  }
  unsigned int m = A.numR;
  unsigned int n = A.numC;
  if(p >= n){
    printf("Cannot extract line %u, %u lines in matrix_t (get_matrix_column).\n", p, n);
    exit(EXIT_FAILURE);
  }
  vector_t v = allocate_vector(m);
  unsigned int i;
  for(i = 0; i < m; i++){
    v.val[i] = A.val[LIN(m, n, i, p)];
  }
  return v;
}

void set_matrix_row(matrix_t A, vector_t v, unsigned int p){
  /*
  set line p of A to v
  */
  if(v.len != A.numC){
    printf("Incompatible dimensions (set_matrix_row).\n");
    exit(EXIT_FAILURE);
  }
  if(v.val == NULL){
    printf("NULL vector_t values pointer (set_matrix_row).\n");
    exit(EXIT_FAILURE);
  }
  if(A.val == NULL){
    printf("NULL matrix_t values pointer (set_matrix_row).\n");
    exit(EXIT_FAILURE);
  }
  if(p >= A.numR){
    printf("Cannot set line %u, %u lines in matrix_t (set_matrix_row).\n", p, A.numR);
    exit(EXIT_FAILURE);
  }
  memcpy(&A.val[(p-1) * A.numC], v.val, sizeof(double) * v.len);
}

void set_matrix_column(matrix_t A, vector_t v, unsigned int p){
  /*
  set column p of A to v
  */
  if(v.len != A.numR){
    printf("Incompatible dimensions (set_matrix_column).\n");
    exit(EXIT_FAILURE);
  }
  if(v.val == NULL){
    printf("NULL vector_t values pointer (set_matrix_column).\n");
    exit(EXIT_FAILURE);
  }
  if(A.val == NULL){
    printf("NULL matrix_t values pointer (set_matrix_column).\n");
    exit(EXIT_FAILURE);
  }
  if(p >= A.numC){
    printf("Cannot set column %u, %u lines in matrix_t (set_matrix_column).\n", p, A.numC);
    exit(EXIT_FAILURE);
  }
  unsigned int i;
  for(i = 0; i < v.len; i++){
    A.val[i * A.numR + p] = v.val[i];
  }
}

vector_t matrix_vector_product(matrix_t A, vector_t v){
  /*
  matrix-vector product
  returns Av
  */
  if(A.val == NULL){
    printf("NULL matrix_t values pointer (matrix_vector_product).\n");
    exit(EXIT_FAILURE);
  }
  if(v.val == NULL){
    printf("NULL vector_t values pointer (matrix_vector_product).\n");
    exit(EXIT_FAILURE);
  }
  if(A.numC != v.len){
    printf("Incompatible dimensions (matrix_vector_product).\n");
    exit(EXIT_FAILURE);
  }
  unsigned int m = A.numR;
  unsigned int n = A.numC;
  vector_t Av = allocate_vector(m);
  unsigned int i, j;
  for(i = 0; i < m; i++){
    Av.val[i] = 0.0;
    for(j = 0; j < n; j++){
      Av.val[i] += v.val[j] * A.val[LIN(m, n, i, j)];
    }
  }
  return Av;
}

double upper_triangular(matrix_t A, vector_t b){
  /*
  transform Ax = b system into upper triangular system
  destroys modifies both A and b
  returns determinant of A
  */
  int n = A.numR;
  double determinant = 1.;
  double pivot;
  int i, j, k;
  for(i = 0; i < n-1; i++){
    if(A.val[LIN(n,n,i,i)] == 0){
      return 0.; //matrix not invertible
    }
    determinant *= A.val[LIN(n,n,i,i)];
    for(j = i+1; j < n; j++){
      pivot = A.val[LIN(n,n,j,i)] / A.val[LIN(n,n,i,i)];
      A.val[LIN(n,n,j,i)] = 0.;
      b.val[j] -= pivot * b.val[i];
      for(k = i+1; k < n; k++){
        A.val[LIN(n,n,j,k)] -= pivot * A.val[LIN(n,n,i,k)];
      }
    }
  }
  determinant *= A.val[LIN(n,n,n-1,n-1)];
  return determinant;
}

double lower_triangular(matrix_t A, vector_t b){
  /*
  transform Ax = b system into lower triangular system
  destroys modifies both A and b
  returns determinant of A
  */
  int n = A.numR;
  double determinant = 1.;
  double pivot;
  int i, j, k;
  for(i = n-1; i > 0; i--){
    if(A.val[LIN(n,n,i,i)] == 0){
      return 0.; //matrix not invertible
    }
    determinant *= A.val[LIN(n,n,i,i)];
    for(j = i-1; j >= 0; j--){
      pivot = A.val[LIN(n,n,j,i)] / A.val[LIN(n,n,i,i)];
      A.val[LIN(n,n,j,i)] = 0.;
      b.val[j] -= pivot * b.val[i];
      for(k = i-1; k >= 0; k--){
        A.val[LIN(n,n,j,k)] -= pivot * A.val[LIN(n,n,i,k)];
      }
    }
  }
  determinant *= A.val[LIN(n,n,0,0)];
  return determinant;
}

vector_t solve_upper_triangular(matrix_t A, vector_t b){
  /*
  A an n by n upper triangular matrix, b n-vector
  solve Ax = b system, return x
  */
  int n = b.len;
  vector_t x = copy_vector(b);
  int i, j;

  x.val[n-1] /= A.val[LIN(n,n,n-1,n-1)];
  for(i = n-2; i >= 0; i--){
    for(j = n-1; j > i; j--){
      x.val[i] -= x.val[j] * A.val[LIN(n,n,i,j)];
    }
    x.val[i] /= A.val[LIN(n,n,i,i)];
  }
  return x;
}

vector_t solve_lower_triangular(matrix_t A, vector_t b){
  /*
  A an n by n lower triangular matrix, b n-vector
  solve Ax = b system, return x
  */
  int n = b.len;
  vector_t x = copy_vector(b);
  int i, j;

  x.val[0] /= A.val[LIN(n,n,0,0)];
  for(i = 1; i < n; i++){
    for(j = 0; j < i; j++){
      x.val[i] -= x.val[j] * A.val[LIN(n,n,i,j)];
    }
    x.val[i] /= A.val[LIN(n,n,i,i)];
  }
  return x;
}

double triangular_matrix_determinant(matrix_t A){
  /*
  A an n x n triangular matrix
  */
  double determinant = 1.;
  int n = A.numR;
  int i;
  for(i = 0; i < n; i++){
    determinant *= A.val[LIN(n,n,i,i)];
  }
  return determinant;
}

int * LUP_decomposition(matrix_t A){
  /*
  A an n by n matrix
  decomposes A into LU form
  using A's memory space
  */
  if(A.val == NULL || A.numR != A.numC){
    printf("Unsuitable matrix for LU decomposition !\n");
    exit(EXIT_FAILURE);
  }
  int n = A.numR;
  double lines_subtracted;
  double max;
  int max_index;
  int aux;
  double * tmp_vector = alloca(sizeof(double) * n);
  int * P = malloc(sizeof(int) * n);
  int i, j, k;
  for(i = 0; i < n; i++){
    P[i] = i;
  }

  for(i = 0; i < n-1; i++){
    max = 0.;
    max_index = i;
    for(j = i; j < n; j++){
      if(fabs(A.val[LIN(n,n,j,i)]) > max){
        max_index = j;
        max = fabs(A.val[LIN(n,n,j,i)]);
      }
    }
    if(A.val[LIN(n,n,max_index,i)] == 0.) continue;
    if(max_index != i){
      aux = P[i];
      P[i] = P[max_index];
      P[max_index] = aux;
      memcpy(tmp_vector, &A.val[LIN(n,n,i,0)], sizeof(double) * n);
      memcpy(&A.val[LIN(n,n,i,0)], &A.val[LIN(n,n,max_index,0)], sizeof(double) * n);
      memcpy(&A.val[LIN(n,n,max_index,0)], tmp_vector, sizeof(double) * n);
    }

    for(j = i+1; j < n; j++){
      lines_subtracted = A.val[LIN(n,n,j,i)] / A.val[LIN(n,n,i,i)];
      A.val[LIN(n,n,j,i)] = lines_subtracted;
      for(k = i+1; k < n; k++){
        A.val[LIN(n,n,j,k)] -= lines_subtracted * A.val[LIN(n,n,i,k)];
      }
    }
  }
  return P;
}

vector_t solve_LUP_system(matrix_t LU, int * P, vector_t b){
  /*
  solves LUx = Pb
  P permutation matrix represented using an array
  */
  int n = LU.numR;
  int i, j;

  //solve Ly = Pb
  vector_t y;
  y.len = n;
  y.val = alloca(sizeof(double) * n);
  for(i = 0; i < n; i++){
    y.val[i] = b.val[P[i]];
  }
  for(i = 1; i < n; i++){
    for(j = 0; j < i; j++){
      y.val[i] -= y.val[j] * LU.val[LIN(n,n,i,j)];
    }
  }

  //solve Ux = y
  vector_t x = copy_vector(y);
  x.val[n-1] /= LU.val[LIN(n,n,n-1,n-1)];
  for(i = n-2; i >= 0; i--){
    for(j = n-1; j > i; j--){
      x.val[i] -= x.val[j] * LU.val[LIN(n,n,i,j)];
    }
    x.val[i] /= LU.val[LIN(n,n,i,i)];
  }
  return x;
}

matrix_t invert_matrix(matrix_t A){
  /*
  A an n by n matrix
  returns inverse of A
  using LUP decomposition
  */
  matrix_t Ac = copy_matrix(A);
  int * P = LUP_decomposition(Ac);
  double determinant = triangular_matrix_determinant(Ac);
  matrix_t A_inv;
  if(determinant == 0.){
    A_inv.val = NULL;
    return A_inv;
  }
  int n = A.numC;
  A_inv = allocate_matrix(n, n);

  vector_t unit_v = zero_vector(n);
  vector_t w;
  int i;
  for(i = 0; i < n; i++){
    unit_v.val[i] = 1.;
    w = solve_LUP_system(Ac, P, unit_v);
    set_matrix_column(A_inv, w, i);
    free_vector(w);
    unit_v.val[i] = 0.;
  }
  free_vector(unit_v);
  return A_inv;
}
