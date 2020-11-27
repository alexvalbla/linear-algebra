
typedef struct {
  int len;
  double * val;
} vector_t;


vector_t allocate_vector(int n);

void free_vector(vector_t v);

void print_vector(vector_t v);

vector_t copy_vector(vector_t v);

vector_t zero_vector(int n);

void set_vector_to_zero(vector_t v);

void add_vector(vector_t u, vector_t v);

vector_t add_vector_ext(vector_t u, vector_t v);

void subtract_vector(vector_t u, vector_t v);

vector_t subtract_vector_ext(vector_t u, vector_t v);

void scale_vector(vector_t v, double scalar);

double scalar_product(vector_t u, vector_t v);

double vector_inf_norm(vector_t v);

double vector_L1_norm(vector_t v);

double vector_L2_norm(vector_t v);

double vector_absolute_difference(vector_t u, vector_t v);
