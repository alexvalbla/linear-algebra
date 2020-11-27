
typedef struct {
  int numR;
  int numC;
  double * val;
} matrix_t;


matrix_t allocate_matrix(int m, int n);

void free_matrix(matrix_t A);

matrix_t zero_matrix(int m, int n);

void set_matrix_to_zero(matrix_t A);

matrix_t identity_matrix(int n);

matrix_t tranpose_matrix(matrix_t A);

int random_matrix(matrix_t A);

matrix_t copy_matrix(matrix_t A);

void print_matrix(matrix_t A);

int add_matrix(matrix_t A, matrix_t B);

matrix_t add_matrix_ext(matrix_t A, matrix_t B);

int subtract_matrix(matrix_t A, matrix_t B);

matrix_t subtract_matrix_ext(matrix_t A, matrix_t B);

int scale_matrix(matrix_t A, double scalar);

int scale_matrix_row(matrix_t A, double scalar, int p);

int scale_matrix_column(matrix_t A, double scalar, int p);

matrix_t matrix_product(matrix_t A, matrix_t B);

double matrix_absolute_difference(matrix_t A, matrix_t B);
