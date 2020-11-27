#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

#define LIN(m, n, i, j) ( (i) * (n) + (j) )
//row i column j on m x n matrix_t


matrix_t allocate_matrix(int m, int n){
  /*
  allocate m x n matrix
  */
  matrix_t A;
  if(m <= 0 || n <= 0){
    A.val = NULL;
    return A;
  }
  A.val = (double *)malloc(sizeof(double) * m * n);
  A.numR = m;
  A.numC = n;
  return A;
}

void free_matrix(matrix_t A){
  /*
  deallocate matrix
  */
  if(A.val != NULL){
    free(A.val);
    A.val = NULL;
  }
}

matrix_t zero_matrix(int m, int n){
  /*
  m x n matrix with values all zeros
  */
  matrix_t Zero;
  if(m <= 0 || n <= 0){
    Zero.val = NULL;
    return Zero;
  }
  Zero.val = (double *)calloc(m * n, sizeof(double));
  Zero.numR = m;
  Zero.numC = n;
  return Zero;
}

void set_matrix_to_zero(matrix_t A){
  /*
  conserves matrix dimensions
  sets all values to zero
  */
  if(A.val == NULL) return;
  memset(A.val, 0, sizeof(double) * A.numR * A.numC);
}

matrix_t identity_matrix(int n){
  /*
  n x n identity matrix
  */
  matrix_t Id;
  if(n <= 0){
    Id.val = NULL;
    return Id;
  }
  Id.numR = Id.numC = n;
  Id.val = (double *)calloc(n * n, sizeof(double));
  int i;
  for(i = 0; i < n; i++){
    Id.val[LIN(n, n, i, i)] = 1.;
  }
  return Id;
}

matrix_t tranpose_matrix(matrix_t A){
  /*
  returns A^t
  */
  matrix_t A_t;
  if(A.val == NULL){
    A_t.val = NULL;
    return A_t;
  }
  int m = A.numR;
  int n = A.numC;
  A_t = allocate_matrix(n, m); // transposed
  if(A_t.val == NULL) return A_t;
  int i, j;
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      A_t.val[LIN(n, m, j, i)] = A.val[LIN(m, n, i, j)];
    }
  }
  return A_t;
}

int random_matrix(matrix_t A){
  if(A.val == NULL) return -1;
  int i;
  for(i = 0; i < A.numR * A.numC; i++){
    A.val[i] = lrand48() / (double)((1 << 31)-1);
  }
  return 0;
}

matrix_t copy_matrix(matrix_t M){
  /*
  copy matrix M into N
  */
  matrix_t N;
  if(M.val == NULL){
    N.val = NULL;
    return N;
  }
  N = allocate_matrix(M.numR, M.numC);
  if(N.val == NULL) return N;
  memcpy(N.val, M.val, sizeof(double) * N.numR * N.numC);
  return N;
}

void print_matrix(matrix_t A){
  if(A.val == NULL){
    printf("EMPTY MATRIX (%u x %u).\n", A.numR, A.numC);
  }
  else{
    int i, j;
    for(i = 0; i < A.numR; i++){
      for(j = 0; j < A.numC; j++){
        printf("%8.10f\t", A.val[LIN(A.numR, A.numC, i, j)]);
      }
      putchar('\n');
      putchar('\n');
    }
  }
}

int add_matrix(matrix_t A, matrix_t B){
  /*
  adds matrix B to A
  in other words: A <-- A + B
  */
  int m = A.numR;
  int n = A.numC;
  if(A.val == NULL || B.val == NULL || m != B.numR || n != B.numC) return -1;
  int i;
  for(i = 0; i < m * n; i++){
    A.val[i] += B.val[i];
  }
  return 0;
}

matrix_t add_matrix_ext(matrix_t A, matrix_t B){
  /*
  adds matrices A and B into C
  in other words: new matrix C = A + B
  */
  matrix_t C;
  int m = A.numR;
  int n = A.numC;
  if(A.val == NULL || B.val == NULL || m != B.numR || n != B.numC){
    C.val = NULL;
    return C;
  }
  C = allocate_matrix(m, n);
  if(C.val == NULL) return C;
  C.numR = m;
  C.numC = n;
  int i;
  for(i = 0; i < m * n; i++){
    C.val[i] = A.val[i] + B.val[i];
  }
  return C;
}

int subtract_matrix(matrix_t A, matrix_t B){
  /*
  subtract matrix B from A
  in other words: A <-- A - B
  */
  int m = A.numR;
  int n = A.numC;
  if(A.val == NULL || B.val == NULL || m != B.numR || n != B.numC) return -1;
  int i;
  for(i = 0; i < m * n; i++){
    A.val[i] -= B.val[i];
  }
  return 0;
}

matrix_t subtract_matrix_ext(matrix_t A, matrix_t B){
  /*
  subtract matrix B from A into C
  in other words: new matrix_t C = A - B
  */
  matrix_t C;
  int m = A.numR;
  int n = A.numC;
  if(A.val == NULL || B.val == NULL || m != B.numR || n != B.numC){
    C.val = NULL;
    return C;
  }
  C = allocate_matrix(m, n);
  if(C.val == NULL) return C;
  C.numR = m;
  C.numC = n;
  int i;
  for(i = 0; i < m * n; i++){
    C.val[i] = A.val[i] + B.val[i];
  }
  return C;
}

int scale_matrix(matrix_t A, double scalar){
  /*
  multiply A by scalar
  */
  int m = A.numR;
  int n = A.numC;
  if(A.val == NULL) return -1;
  int i;
  for(i = 0; i < m * n; i++){
    A.val[i] *= scalar;
  }
  return 0;
}

int scale_matrix_row(matrix_t A, double scalar, int p){
  /*
  multiply A row p by scalar
  */
  int m = A.numR;
  int n = A.numC;
  if(A.val == NULL || p < 0 || p >= m) return -1;
  int i;
  for(i = 0; i < n; i++){
    A.val[LIN(m,n,p,i)] *= scalar;
  }
  return 0;
}

int scale_matrix_column(matrix_t A, double scalar, int p){
  /*
  multiply A column p by scalar
  */
  int m = A.numR;
  int n = A.numC;
  if(p < 0 || p >= n || A.val == NULL) return -1;
  int i;
  for(i = 0; i < m; i++){
    A.val[LIN(m,n,i,p)] *= scalar;
  }
  return 0;
}

matrix_t matrix_product(matrix_t A, matrix_t B){
  /*
  new matrix C = A x B
  */
  int m = A.numR;
  int n = A.numC;
  int p = B.numC;
  matrix_t C;
  if(A.val == NULL || B.val == NULL || n != B.numR){
    C.val = NULL;
    return C;
  }
  C = allocate_matrix(m, p);
  if(C.val == NULL) return C;
  int i, j, k;
  //temporary vector tmpV
  double * tmpV = alloca(sizeof(double) * n);
  for(j = 0; j < p; j++){
    //store column j of B in tmpV
    for(k = 0; k < n; k++){
      tmpV[k] = B.val[LIN(n, p, k, j)];
    }
    //for each (contiguous) line in A...
    for(i = 0; i < m; i++){
      //multiply line i of A by tmpV
      C.val[LIN(m, p, i, j)] = 0.0;
      for(k = 0; k < n; k++){
        C.val[LIN(m, p, i, j)] += A.val[LIN(m, n, i, k)] * tmpV[k];
      }
    }
  }
  return C;
}

double matrix_absolute_difference(matrix_t A, matrix_t B){
  /*
  returns highest absolute difference
  between values of A and B
  incompatible dimensions returns NaN
  */
  int m = A.numR;
  int n = A.numC;
  if(A.val == NULL || B.val == NULL || m != B.numR || n != B.numC){
    volatile double error_code = 0./0.; // NaN
    return error_code;
  }
  double max = 0.;
  double d;
  int i, j;
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      if((d = fabs(A.val[LIN(m,n,i,j)] - B.val[LIN(m,n,i,j)])) > max) max = d;
    }
  }
  return max;
}






void splitVertical(matrix_t A, matrix_t A1, matrix_t A2){
  /*
  A1 and A2 already allocated
  */
  int m = A.numR;
  int n = A.numC;
  int i;
  for(i = 0; i < m; i++){
    memcpy(&A.val[LIN(m, n, i, 0)], &A1.val[LIN(m, n/2, i, 0)], sizeof(double) * (A.numC/2));
    memcpy(&A.val[LIN(m, n, i, n/2)], &A2.val[LIN(m, n/2, i, 0)], sizeof(double) * (A.numC/2 + (A.numC&1)));
  }
}

static void strassenSplit(matrix_t A, matrix_t A11, matrix_t A12, matrix_t A21, matrix_t A22){
  /*
  A is 2^k x 2^k matrix
  A1, A2, A3, A4 already allocated
  are size 2^(k-1) x 2^(k-1)
  */
  int n = A.numR;// == A.numC
  int i;
  for(i = 0; i < n/2; i++){
    memcpy(&A11.val[LIN(n/2, n/2, i, 0)], &A.val[LIN(n, n, i, 0)], sizeof(double) * (n/2));
    memcpy(&A12.val[LIN(n/2, n/2, i, 0)], &A.val[LIN(n, n, i, n/2)], sizeof(double) * (n/2));
  }
  for(i = 0; i < n/2; i++){
    memcpy(&A21.val[LIN(n/2, n/2, i, 0)], &A.val[LIN(n, n, i+n/2, 0)], sizeof(double) * (n/2));
    memcpy(&A22.val[LIN(n/2, n/2, i, 0)], &A.val[LIN(n, n, i+n/2, n/2)], sizeof(double) * (n/2));
  }
}

static void strassenAssemble(matrix_t C, matrix_t C11, matrix_t C12, matrix_t C21, matrix_t C22){
  /*
  C is 2^k x 2^k matrix
  C1, C2, C3, C4 already allocated
  are size 2^(k-1) x 2^(k-1)
  */
  int n = C.numR;// == C.numC
  int i;
  for(i = 0; i < n/2; i++){
    memcpy(&C.val[LIN(n, n, i, 0)], &C11.val[LIN(n/2, n/2, i, 0)], sizeof(double) * (n/2));
    memcpy(&C.val[LIN(n, n, i, n/2)], &C12.val[LIN(n/2, n/2, i, 0)], sizeof(double) * (n/2));
  }
  for(i = 0; i < n/2; i++){
    memcpy(&C.val[LIN(n, n, i+n/2, 0)], &C21.val[LIN(n/2, n/2, i, 0)], sizeof(double) * (n/2));
    memcpy(&C.val[LIN(n, n, i+n/2, n/2)], &C22.val[LIN(n/2, n/2, i, 0)], sizeof(double) * (n/2));
  }
}

matrix_t strassenRec(matrix_t A, matrix_t B){
  /*
  A and B are size 2^k x 2^k
  n = 2^k
  */
  int n = A.numR;
  if(n < 100){
    matrix_t C = matrix_product(A, B);
    return C;
  }

  matrix_t A11 = allocate_matrix(n/2, n/2);
  matrix_t A12 = allocate_matrix(n/2, n/2);
  matrix_t A21 = allocate_matrix(n/2, n/2);
  matrix_t A22 = allocate_matrix(n/2, n/2);
  strassenSplit(A, A11, A12, A21, A22);

  matrix_t B11 = allocate_matrix(n/2, n/2);
  matrix_t B12 = allocate_matrix(n/2, n/2);
  matrix_t B21 = allocate_matrix(n/2, n/2);
  matrix_t B22 = allocate_matrix(n/2, n/2);
  strassenSplit(B, B11, B12, B21, B22);

  matrix_t tmp1 = add_matrix_ext(A11, A22);
  matrix_t tmp2 = add_matrix_ext(B11, B22);
  matrix_t A1 = strassenRec(tmp1, tmp2);
  free_matrix(tmp1);
  free_matrix(tmp2);

  tmp1 = add_matrix_ext(A21, A22);
  matrix_t A2 = strassenRec(tmp1, B11);
  free_matrix(tmp1);

  tmp1 = subtract_matrix_ext(B12, B22);
  matrix_t A3 = strassenRec(A11, tmp1);
  free_matrix(tmp1);

  tmp1 = subtract_matrix_ext(B21, B11);
  matrix_t A4 = strassenRec(A22, tmp1);
  free_matrix(tmp1);

  tmp1 = add_matrix_ext(A11, A12);
  matrix_t A5 = strassenRec(tmp1, B22);
  free_matrix(tmp1);

  tmp1 = subtract_matrix_ext(A21, A11);
  tmp2 = add_matrix_ext(B11, B12);
  matrix_t A6 = strassenRec(tmp1, tmp2);
  free_matrix(tmp1);
  free_matrix(tmp2);

  tmp1 = subtract_matrix_ext(A12, A22);
  tmp1 = add_matrix_ext(B21, B22);
  matrix_t A7 = strassenRec(tmp1, tmp2);
  free_matrix(tmp1);
  free_matrix(tmp2);


  free_matrix(A11);
  free_matrix(A12);
  free_matrix(A21);
  free_matrix(A22);
  free_matrix(B11);
  free_matrix(B12);
  free_matrix(B21);
  free_matrix(B22);


  matrix_t C11 = add_matrix_ext(A1, A4);
  subtract_matrix(C11, A5);
  add_matrix(C11, A7);
  matrix_t C12 = add_matrix_ext(A3, A5);
  matrix_t C21 = add_matrix_ext(A2, A4);
  matrix_t C22 = subtract_matrix_ext(A1, A2);
  add_matrix(C22, A3);
  add_matrix(C22, A6);


  free_matrix(A1);
  free_matrix(A2);
  free_matrix(A3);
  free_matrix(A4);
  free_matrix(A5);
  free_matrix(A6);
  free_matrix(A7);


  matrix_t C = allocate_matrix(n, n);
  strassenAssemble(C, C11, C12, C21, C22);
  free_matrix(C11);
  free_matrix(C12);
  free_matrix(C21);
  free_matrix(C22);


  return C;
}
