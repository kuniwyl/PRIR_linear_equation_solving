#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gauss_seidel.h"

void gauss_seidel(double** A, double* b, double* x0, double epsilon, int n, int max_iterations) {
    double* x = malloc(sizeof(double) * n);
    int i, j, k;
    double sigma;
    double diff;

    for (k = 0; k < max_iterations; k++) {
        for (i = 0; i < n; i++) {
            sigma = 0.0;
            for (j = 0; j < n; j++) {
                if (j != i) {
                    sigma += A[i][j] * x0[j];
                }
            }
            x[i] = (b[i] - sigma) / A[i][i];
        }

        diff = 0.0;
        for (i = 0; i < n; i++) {
            diff += fabs(x[i] - x0[i]);
            x0[i] = x[i];
        }
        if (diff < epsilon) {
            return;
        }
    }
}
