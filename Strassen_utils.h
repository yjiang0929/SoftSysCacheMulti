const int MIN_SIZE;

typedef struct {
    double *arr;
    int size;
} Matrix;

Matrix *make_matrix(int size);
Matrix *to_matrix(double *a, int size);
void free_matrix(Matrix *a);
void print_mat(Matrix *a);
Matrix *sum_matrix(Matrix *a, Matrix *b);
Matrix *subtract_matrix(Matrix *a, Matrix *b);
Matrix *mult_matrix(Matrix *a, Matrix *b);
Matrix *compute_c11(Matrix *a, Matrix *b, Matrix *c, Matrix *d);
Matrix *compute_c22(Matrix *a, Matrix *b, Matrix *c, Matrix *d);
Matrix *subdivide(Matrix *a, int start_row, int start_col);
Matrix *merge(Matrix *a, Matrix *b, Matrix *c, Matrix *d);

