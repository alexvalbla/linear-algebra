#include "matrix.h"
#include "vector.h"



vector_t get_matrix_row(matrix_t A, unsigned int p);

vector_t get_matrix_column(matrix_t A, unsigned int p);

void set_matrix_row(matrix_t A, vector_t v, unsigned int p);

void set_matrix_column(matrix_t A, vector_t v, unsigned int p);

vector_t matrix_vector_product(matrix_t A, vector_t v);

double upper_triangular(matrix_t A, vector_t b);

double lower_triangular(matrix_t A, vector_t b);

vector_t solve_upper_triangular(matrix_t A, vector_t b);

vector_t solve_lower_triangular(matrix_t A, vector_t b);

double triangular_matrix_determinant(matrix_t A);

int * LUP_decomposition(matrix_t A);

vector_t solve_LUP_system(matrix_t LU, int * P, vector_t b);

matrix_t invert_matrix(matrix_t A);
