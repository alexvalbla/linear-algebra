#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vector.h"




vector_t allocate_vector(int n){
  /*
  allocate n-length vector
  */
  vector_t v;
  if(n <= 0){
    v.val = NULL;
    return v;
  }
  v.len = n;
  v.val = (double *)malloc(sizeof(double) * n);
  if(v.val == NULL){
    printf("Vector allocation failed (allocate_vector).\n");
    exit(EXIT_FAILURE);
  }
  return v;
}

void free_vector(vector_t v){
  /*
  deallocate vector
  */
  if(v.val != NULL){
    free(v.val);
    v.val = NULL;
  }
}

void print_vector(vector_t v){
  int n = v.len;
  int i;
  for(i = 0; i < n-1; i++){
    printf("%.4f\t", v.val[i]);
  }
  printf("%.4f\n", v.val[n-1]);
}

vector_t copy_vector(vector_t v){
  /*
  copy vector into w
  */
  vector_t w;
  w = allocate_vector(v.len);
  if(v.val == NULL) return w;
  memcpy(w.val, v.val, sizeof(double) * w.len);
  return w;
}

vector_t zero_vector(int n){
  /*
  n length vector with values all zero
  */
  vector_t v;
  if(n <= 0){
    v.val = NULL;
    return v;
  }
  v.len = n;
  v.val = (double *)calloc(n, sizeof(double));
  if(v.val == NULL){
    printf("Vector allocation failed (zero_vector).\n");
    exit(EXIT_FAILURE);
  }
  return v;
}

void set_vector_to_zero(vector_t v){
  if(v.val == NULL){
    printf("Cannot wipe unallocated vector_t (wipe_vector).\n");
    exit(EXIT_FAILURE);
  }
  memset(v.val, 0, sizeof(double) * v.len);
}

void add_vector(vector_t u, vector_t v){
  /*
  adds vector v to u
  in other words: u <-- u + v
  */
  if(u.val == NULL || v.val == NULL){
    printf("NULL vector_t values pointer (add_vector).\n");
    exit(EXIT_FAILURE);
  }
  if(u.len != v.len){
    printf("Incompatible vector_t dimensions (add_vector).\n");
    exit(EXIT_FAILURE);
  }
  int i;
  for(i = 0; i < u.len; i++){
    u.val[i] += v.val[i];
  }
}

vector_t add_vector_ext(vector_t u, vector_t v){
  /*
  adds vectors u and v into w
  in other words: new vector w = u + v
  */
  if(u.val == NULL || v.val == NULL){
    printf("NULL vector_t values pointer (add_vector_ext).\n");
    exit(EXIT_FAILURE);
  }
  if(u.len != v.len){
    printf("Incompatible vector_t dimensions (add_vector_ext).\n");
    exit(EXIT_FAILURE);
  }
  int n = u.len;
  vector_t w = zero_vector(n);
  int i;
  for(i = 0; i < n; i++){
    w.val[i] = u.val[i] + v.val[i];
  }
  return w;
}

void subtract_vector(vector_t u, vector_t v){
  /*
  subtracts vector v from u
  in other words: u <-- u - v
  */
  if(u.val == NULL || v.val == NULL){
    printf("NULL vector_t values pointer (subtract_vector).\n");
    exit(EXIT_FAILURE);
  }
  if(u.len != v.len){
    printf("Incompatible vector_t dimensions (subtract_vector).\n");
    exit(EXIT_FAILURE);
  }
  int i;
  for(i = 0; i < u.len; i++){
    u.val[i] -= v.val[i];
  }
}

vector_t subtract_vector_ext(vector_t u, vector_t v){
  /*
  subtracts vector v from u into w
  in other words: new vector w = u - v
  */
  if(u.val == NULL || v.val == NULL){
    printf("NULL vector_t values pointer (subtract_vector_ext).\n");
    exit(EXIT_FAILURE);
  }
  if(u.len != v.len){
    printf("Incompatible vector_t dimensions (subtract_vector_ext).\n");
    exit(EXIT_FAILURE);
  }
  int n = u.len;
  vector_t w = zero_vector(n);
  int i;
  for(i = 0; i < n; i++){
    w.val[i] = u.val[i] - v.val[i];
  }
  return w;
}

void scale_vector(vector_t v, double scalar){
  /*
  multiply v by scalar
  */
  if(v.val == NULL){
    printf("NULL vector_t values pointer (scale_vector).\n");
    exit(EXIT_FAILURE);
  }
  int i;
  for(i = 0; i < v.len; i++){
    v.val[i] *= scalar;
  }
}

double scalar_product(vector_t u, vector_t v){
  /*
  return u * v scalar product
  */
  if(u.val == NULL || v.val == NULL){
    printf("NULL vector_t values pointer (scalar_product).\n");
    exit(EXIT_FAILURE);
  }
  int n = u.len;
  if(v.len != n){
    printf("Incompatible dimensions (scalar_product).\n");
    exit(EXIT_FAILURE);
  }
  int i;
  double sum = 0.0;
  for(i = 0; i < n; i++){
    sum += u.val[i] * v.val[i];
  }
  return sum;
}

double vector_inf_norm(vector_t v){
  /*
  return supremum norm of v
  */
  int n = v.len;
  double max = 0;
  double d;
  int i;
  for(i = 0; i < n; i++){
    if((d = fabs(v.val[i])) > max) max = d;
  }
  return max;
}

double vector_L1_norm(vector_t v){
  /*
  return Manhattan norm of v
  */
  int n = v.len;
  double sum = 0;
  int i;
  for(i = 0; i < n; i++){
    sum += fabs(v.val[i]);
  }
  return sum;
}

double vector_L2_norm(vector_t v){
  /*
  return Euclidian norm of v
  */
  int n = v.len;
  double sum = 0;
  int i;
  for(i = 0; i < n; i++){
    sum += v.val[i] * v.val[i];
  }
  return sqrt(sum);
}

double vector_absolute_difference(vector_t u, vector_t v){
  /*
  returns highest absolute difference
  between values of u and v
  u and v of equal dimensions
  */
  if(u.len != v.len){
    printf("Incompatible vector_t dimensions (vector_absolute_difference).\n");
    exit(EXIT_FAILURE);
  }
  int n = u.len;
  double max = 0.;
  double d;
  int i;
  for(i = 0; i < n; i++){
    if((d = fabs(u.val[i] - v.val[i])) > max) max = d;
  }
  return max;
}
