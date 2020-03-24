/* Create macros so that the matrices are stored in column-major order */
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#define A(i, j) a->arr[ (i)*a->size + (j) ]
#define B(i, j) b->arr[ (i)*b->size + (j) ]
#define C(i, j) c->arr[ (i)*c->size + (j) ]
#define SiC(i, j) si->c->arr[ (i)*si->c->size + (j) ]
#define D(i, j) d->arr[ (i)*d->size + (j) ]
#define R(i, j) r->arr[ (i)*r->size + (j) ]

const int MIN_SIZE = 8;

typedef struct {
    double *arr;
    int size;
} Matrix;

typedef struct {
    Matrix *a;
    Matrix *b;
    Matrix *c;
    int isFirst;
} StrassenInput;

/* Allocate space for new matrix */
Matrix *make_matrix(int size) {
    Matrix *new = malloc(sizeof(Matrix));
    new->size = size;
    new->arr = (double *) malloc(size * size * sizeof(double));
    return new;
}

/* Convert array to Matrix struct */
Matrix *to_matrix(double *a, int size) {
    Matrix *new = malloc(sizeof(Matrix));
    new->size = size;
    new->arr = a;
    return new;
}

/* Free matrix array and struct */
void free_matrix(Matrix *a) {
    free(a->arr);
    free(a);
}

StrassenInput *make_strassen_input(Matrix *a, Matrix *b, uint isFirst) {
    StrassenInput *new = malloc(sizeof(StrassenInput));
    new->a = a;
    new->b = b;
    new->isFirst = isFirst;
    return new;
}

/* Free Strassen Input */
void free_strassen_inputs(StrassenInput **si) {
    free(si[0]);
    free(si[1]);
    free(si[2]);
    free(si[3]);
    free(si[4]);
    free(si[5]);
    free(si[6]);
    free(si);
}

/* Print elements of matrix */
void print_mat(Matrix *a) {
    for (int i = 0; i < a->size; i++) {
        for (int j = 0; j < a->size; j++) {
            printf("%f\t", A(i, j));
        }
        printf("\n");
    }
    printf("\n");
}

/* Element-wise add Matrix a and b */
Matrix *sum_matrix(Matrix *a, Matrix *b) {
    Matrix *c = make_matrix(a->size);

    for (int i = 0; i < c->size; i++) {
        for (int j = 0; j < c->size; j++) {
            C(i, j) = A(i, j) + B(i, j);
        }
    }

    return c;
}

/* Element-wise subtract Matrix b from a */
Matrix *subtract_matrix(Matrix *a, Matrix *b) {
    Matrix *c = make_matrix(a->size);

    for (int i = 0; i < c->size; i++) {
        for (int j = 0; j < c->size; j++) {
            C(i, j) = A(i, j) - B(i, j);
        }
    }

    return c;
}

/* Multiply Matrix a and b
   Only used for matrices of small size */
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

/* According to Strassen algorithm, compute C11 block with Matrix a,b,c,d
    C11 = A + B - C + D*/
Matrix *compute_c11(Matrix *a, Matrix *b, Matrix *c, Matrix *d) {
    Matrix *r = make_matrix(a->size);

    for (int i = 0; i < c->size; i++) {
        for (int j = 0; j < c->size; j++) {
            R(i, j) = A(i, j) + B(i, j) - C(i, j) + D(i, j);
        }
    }

    return r;
}

/* According to Strassen algorithm, compute C22 block with Matrix a,b,c,d
    C22 = A - B + C + D*/
Matrix *compute_c22(Matrix *a, Matrix *b, Matrix *c, Matrix *d) {
    Matrix *r = make_matrix(a->size);

    for (int i = 0; i < c->size; i++) {
        for (int j = 0; j < c->size; j++) {
            R(i, j) = A(i, j) - B(i, j) + C(i, j) + D(i, j);
        }
    }

    return r;
}

/* Divide input matrix by half from start_row and start_col
    Returns the divided matrix */
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

/* Merge four small matrices into a full matrix */
Matrix *merge(Matrix *a, Matrix *b, Matrix *c, Matrix *d) {
    int size = a->size * 2;
    Matrix *r = malloc(sizeof(Matrix));
    r->size = size;
    r->arr = (double *) malloc(size * size * sizeof(double));

    // c11
    int index = 0;
    for (int i = 0; i < size / 2; i++) {
        for (int j = 0; j < size / 2; j++) {
            R(i, j) = a->arr[index++];
        }
    }

    // c12
    index = 0;
    for (int i = 0; i < size / 2; i++) {
        for (int j = size / 2; j < size; j++) {
            R(i, j) = b->arr[index++];
        }
    }

    // c21
    index = 0;
    for (int i = size / 2; i < size; i++) {
        for (int j = 0; j < size / 2; j++) {
            R(i, j) = c->arr[index++];
        }
    }

    // c22
    index = 0;
    for (int i = size / 2; i < size; i++) {
        for (int j = size / 2; j < size; j++) {
            R(i, j) = d->arr[index++];
        }
    }

    return r;
}

/* Matrix multiplication with Strassen algorithm
    Multi-threading is used in the first level of recursion for parallelization
    and is avoided further to prevent generating too many threads */
void *Strassen_MMult(void* s) {
    Matrix* matrix_a = ((StrassenInput *)s)->a;
    Matrix* matrix_b = ((StrassenInput *)s)->b;

    int size = matrix_a->size;
    // Base case
    if (size <= MIN_SIZE) {
        ((StrassenInput *)s)->c = mult_matrix(matrix_a, matrix_b);
        return NULL;
    }

    // Sub-divide matrices A and B
    Matrix *a11 = subdivide(matrix_a, 0, 0);
    Matrix *a12 = subdivide(matrix_a, 0, size / 2);
    Matrix *a21 = subdivide(matrix_a, size / 2, 0);
    Matrix *a22 = subdivide(matrix_a, size / 2, size / 2);

    Matrix *b11 = subdivide(matrix_b, 0, 0);
    Matrix *b12 = subdivide(matrix_b, 0, size / 2);
    Matrix *b21 = subdivide(matrix_b, size / 2, 0);
    Matrix *b22 = subdivide(matrix_b, size / 2, size / 2);

    // add and subtract matrix a and b
    Matrix *a11_p_a22 = sum_matrix(a11, a22);
    Matrix *b11_p_b22 = sum_matrix(b11, b22);
    Matrix *a21_p_a22 = sum_matrix(a21, a22);
    Matrix *b12_s_b22 = subtract_matrix(b12, b22);
    Matrix *b21_s_b11 = subtract_matrix(b21, b11);
    Matrix *a11_p_a12 = sum_matrix(a11, a12);
    Matrix *a21_s_a11 = subtract_matrix(a21, a11);
    Matrix *b11_p_b12 = sum_matrix(b11, b12);
    Matrix *a12_s_a22 = subtract_matrix(a12, a22);
    Matrix *b21_p_b22 = sum_matrix(b21, b22);

    // Relation recursion with multi threading
    StrassenInput **si = malloc(7 * sizeof(StrassenInput *));
    si[0] = make_strassen_input(a11_p_a22, b11_p_b22, 0);
    si[1] = make_strassen_input(b11, a21_p_a22, 0);
    si[2] = make_strassen_input(a11, b12_s_b22, 0);
    si[3] = make_strassen_input(a22, b21_s_b11, 0);
    si[4] = make_strassen_input(a11_p_a12, b22, 0);
    si[5] = make_strassen_input(a21_s_a11, b11_p_b12, 0);
    si[6] = make_strassen_input(a12_s_a22, b21_p_b22, 0);

    int i;
    if (((StrassenInput *)s)->isFirst) {
        pthread_t strassen_thread[7];
        for (i=0; i<7; i++) {
            if (pthread_create(&strassen_thread[i], NULL, Strassen_MMult, si[i])) {
                fprintf(stderr, "Error creating thread %d\n", i);
                exit(1);
            }
        }
        for (i=0; i<7; i++) {
            if (pthread_join(strassen_thread[i], NULL)) {
                fprintf(stderr, "Error joining thread %d\n", i);
                exit(2);
            }
        }
    } else {
        for (i=0;i<7;i++) {
            Strassen_MMult(si[i]);
        }
    }

    Matrix *p1 = si[0]->c;
    Matrix *p2 = si[1]->c;
    Matrix *p3 = si[2]->c;
    Matrix *p4 = si[3]->c;
    Matrix *p5 = si[4]->c;
    Matrix *p6 = si[5]->c;
    Matrix *p7 = si[6]->c;

    // free intermediate matrices
    free_matrix(a11);
    free_matrix(a12);
    free_matrix(a21);
    free_matrix(a22);
    free_matrix(b11);
    free_matrix(b12);
    free_matrix(b21);
    free_matrix(b22);

    free_matrix(a11_p_a22);
    free_matrix(b11_p_b22);
    free_matrix(a21_p_a22);
    free_matrix(b12_s_b22);
    free_matrix(b21_s_b11);
    free_matrix(a11_p_a12);
    free_matrix(a21_s_a11);
    free_matrix(b11_p_b12);
    free_matrix(a12_s_a22);
    free_matrix(b21_p_b22);

    free_strassen_inputs(si);

    // Merge
    Matrix *c11 = compute_c11(p1, p4, p5, p7);
    Matrix *c12 = sum_matrix(p3, p5);
    Matrix *c21 = sum_matrix(p2, p4);
    Matrix *c22 = compute_c22(p1, p2, p3, p6);

    free_matrix(p1);
    free_matrix(p2);
    free_matrix(p3);
    free_matrix(p4);
    free_matrix(p5);
    free_matrix(p6);
    free_matrix(p7);

    ((StrassenInput *)s)->c = merge(c11, c12, c21, c22);

    free_matrix(c11);
    free_matrix(c12);
    free_matrix(c21);
    free_matrix(c22);

    return NULL;
}

void MY_MMult(int m, int n, int k, double *a, int lda,
              double *b, int ldb,
              double *c_r, int ldc) {

    Matrix *matrix_a = to_matrix(a, lda);
    Matrix *matrix_b = to_matrix(b, ldb);

    StrassenInput *si = make_strassen_input(matrix_a, matrix_b, 1);
    Strassen_MMult(si);

    for (int i = 0; i < si->c->size; i++) {
        for (int j = 0; j < si->c->size; j++) {
            c_r[i * ldc + j] = SiC(i, j);
        }
    }

    free(si);
}
