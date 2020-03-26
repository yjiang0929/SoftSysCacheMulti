#include <stdlib.h>
#include <stdio.h>

/* Create macros so that the matrices are stored in column-major order */
#define A(i, j) a->arr[ (i)*a->size + (j) ]
#define B(i, j) b->arr[ (i)*b->size + (j) ]
#define C(i, j) c->arr[ (i)*c->size + (j) ]
#define D(i, j) d->arr[ (i)*d->size + (j) ]
#define R(i, j) r->arr[ (i)*r->size + (j) ]

const int MIN_SIZE = 8;

typedef struct {
    double *arr;
    int size;
} Matrix;

/**
 * Allocate space for a new matrix
 * @param size: size of the matrix
 * @return: a newly allocated matrix
 */
Matrix *make_matrix(int size) {
    Matrix *new = malloc(sizeof(Matrix));
    new->size = size;
    new->arr = (double *) malloc(size * size * sizeof(double));
    return new;
}

/**
 * Convert array to Matrix struct
 * @param a: 1D array that represents a 2D matrix
 * @param size: size of the matrix
 * @return: a newly allocated matrix
 */
Matrix *to_matrix(double *a, int size) {
    Matrix *new = malloc(sizeof(Matrix));
    new->size = size;
    new->arr = a;
    return new;
}

/**
 * Free matrix array and struct
 * @param a: input matrix
 */
void free_matrix(Matrix *a) {
    free(a->arr);
    free(a);
}

/**
 * Print elements of matrix
 * @param a: input matrix
 */
void print_mat(Matrix *a) {
    for (int i = 0; i < a->size; i++) {
        for (int j = 0; j < a->size; j++) {
            printf("%f\t", A(i, j));
        }
        printf("\n");
    }
    printf("\n");
}

/**
 * Element-wise summation of Matrix a and b
 * @param a: input matrix a
 * @param b: input matrix b
 * @return: a newly allocated matrix
 */
Matrix *sum_matrix(Matrix *a, Matrix *b) {
    Matrix *c = make_matrix(a->size);

    for (int i = 0; i < c->size; i++) {
        for (int j = 0; j < c->size; j++) {
            C(i, j) = A(i, j) + B(i, j);
        }
    }

    return c;
}

/**
 * Element-wise subtract Matrix b from a
 * @param a: input matrix a
 * @param b: input matrix b
 * @return: a newly allocated matrix
 */
Matrix *subtract_matrix(Matrix *a, Matrix *b) {
    Matrix *c = make_matrix(a->size);

    for (int i = 0; i < c->size; i++) {
        for (int j = 0; j < c->size; j++) {
            C(i, j) = A(i, j) - B(i, j);
        }
    }

    return c;
}

/**
 * Multiply Matrix a and b
 *  Only used for matrices of small size
 * @param a: input matrix a
 * @param b: input matrix b
 * @return: a newly allocated matrix
 */
Matrix *mult_matrix(Matrix *a, Matrix *b) {
    Matrix *c = make_matrix(a->size);

    for (int i = 0; i < c->size; i++) {
        for (int j = 0; j < c->size; j++) {
            for (int p = 0; p < c->size; p++) {
                C(i, j) += A(i, p) * B(p, j);
            }
        }
    }

    return c;
}

/**
 * According to Strassen algorithm, compute C11 block with Matrix a, b, c, d
 *  C11 = A + B - C + D
 * @param a: input matrix a
 * @param b: input matrix b
 * @param c: input matrix c
 * @param d: input matrix d
 * @return: a newly allocated matrix
 */
Matrix *compute_c11(Matrix *a, Matrix *b, Matrix *c, Matrix *d) {
    Matrix *r = make_matrix(a->size);

    for (int i = 0; i < c->size; i++) {
        for (int j = 0; j < c->size; j++) {
            R(i, j) = A(i, j) + B(i, j) - C(i, j) + D(i, j);
        }
    }

    return r;
}

/**
 * According to Strassen algorithm, compute C22 block with Matrix a, b, c, d
 *  C22 = A - B + C + D
 * @param a: input matrix a
 * @param b: input matrix b
 * @param c: input matrix c
 * @param d: input matrix d
 * @return: a newly allocated matrix
 */
Matrix *compute_c22(Matrix *a, Matrix *b, Matrix *c, Matrix *d) {
    Matrix *r = make_matrix(a->size);

    for (int i = 0; i < c->size; i++) {
        for (int j = 0; j < c->size; j++) {
            R(i, j) = A(i, j) - B(i, j) + C(i, j) + D(i, j);
        }
    }

    return r;
}

/**
 * Divide input matrix by half from start_row and start_col.
 * Returns the divided matrix.
 * @param a: input matrix a
 * @param start_row: the index of the start row
 * @param start_col: the index of the start column
 * @return: a newly allocated matrix
 */
Matrix *subdivide(Matrix *a, int start_row, int start_col) {
    int size = a->size / 2;

    if (size < MIN_SIZE) {
        printf("Trying to divide matrix smaller than MIN_SIZE = %d\n", MIN_SIZE);
        exit(1);
    }

    Matrix *new = malloc(sizeof(Matrix));
    new->size = size;
    new->arr = (double *) malloc(size * size * sizeof(double));

    int end_row = start_row + size;
    int end_col = start_col + size;
    int new_index = 0;
    for (int i = start_row; i < end_row; i++) {
        for (int j = start_col; j < end_col; j++) {
            new->arr[new_index++] = A(i, j);
        }
    }
    return new;
}

/**
 * Merge four small matrices into a full matrix
 * @param a: input matrix a
 * @param b: input matrix b
 * @param c: input matrix c
 * @param d: input matrix d
 * @return: a newly allocated merged matrix
 */
Matrix *merge(Matrix *a, Matrix *b, Matrix *c, Matrix *d) {
    int size = a->size * 2;
    int half_size = size / 2;
    Matrix *r = malloc(sizeof(Matrix));
    r->size = size;
    r->arr = (double *) malloc(size * size * sizeof(double));

    // c11
    int index = 0;
    for (int i = 0; i < half_size; i++) {
        for (int j = 0; j < half_size; j++) {
            R(i, j) = a->arr[index++];
        }
    }

    // c12
    index = 0;
    for (int i = 0; i < half_size; i++) {
        for (int j = half_size; j < size; j++) {
            R(i, j) = b->arr[index++];
        }
    }

    // c21
    index = 0;
    for (int i = half_size; i < size; i++) {
        for (int j = 0; j < half_size; j++) {
            R(i, j) = c->arr[index++];
        }
    }

    // c22
    index = 0;
    for (int i = half_size; i < size; i++) {
        for (int j = half_size; j < size; j++) {
            R(i, j) = d->arr[index++];
        }
    }

    return r;
}