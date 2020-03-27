#include <stdlib.h>
#include <stdio.h>

#include "Strassen_utils.h"

#define C(i, j) c->arr[ (i)*c->size + (j) ]

/**
 * Matrix multiplication with Strassen algorithm
 * @param matrix_a: input matrix a
 * @param matrix_b: input matrix c
 * @return: a newly allocated resulting matrix
 */
Matrix *Strassen_MMult(Matrix *matrix_a, Matrix *matrix_b) {
    int size = matrix_a->size;
    int half_size = size / 2;

    // Base case
    if (size <= MIN_SIZE) {
        return mult_matrix(matrix_a, matrix_b);
    }

    // Sub-divide matrices A and B
    Matrix *a11 = subdivide(matrix_a, 0, 0);
    Matrix *a12 = subdivide(matrix_a, 0, half_size);
    Matrix *a21 = subdivide(matrix_a, half_size, 0);
    Matrix *a22 = subdivide(matrix_a, half_size, half_size);

    Matrix *b11 = subdivide(matrix_b, 0, 0);
    Matrix *b12 = subdivide(matrix_b, 0, half_size);
    Matrix *b21 = subdivide(matrix_b, half_size, 0);
    Matrix *b22 = subdivide(matrix_b, half_size, half_size);

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

    // Relation recursion
    Matrix *p1 = Strassen_MMult(a11_p_a22, b11_p_b22);
    Matrix *p2 = Strassen_MMult(b11, a21_p_a22);
    Matrix *p3 = Strassen_MMult(a11, b12_s_b22);
    Matrix *p4 = Strassen_MMult(a22, b21_s_b11);
    Matrix *p5 = Strassen_MMult(a11_p_a12, b22);
    Matrix *p6 = Strassen_MMult(a21_s_a11, b11_p_b12);
    Matrix *p7 = Strassen_MMult(a12_s_a22, b21_p_b22);

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

    Matrix *r = merge(c11, c12, c21, c22);

    free_matrix(c11);
    free_matrix(c12);
    free_matrix(c21);
    free_matrix(c22);

    return r;
}

void MY_MMult(int m, int n, int k, double *a, int lda,
              double *b, int ldb,
              double *c_r, int ldc) {
    Matrix *matrix_a = to_matrix(a, lda);
    Matrix *matrix_b = to_matrix(b, ldb);

    Matrix *c = Strassen_MMult(matrix_a, matrix_b);

    for (int i = 0; i < c->size; i++) {
        for (int j = 0; j < c->size; j++) {
            c_r[i * ldc + j] = C(i, j);
        }
    }
}
