#include <stdio.h>
#include <math.h>

#include "gaussian_elimination.h"

void gaussian_elimination(double** A, double* b, double* x, int n) {
    int i, j, k;
    double factor;

    for (i = 0; i < n; i++) {
        if (A[i][i] == 0) {
            printf("Macierz nie ma unikalnego rozwiązania.\n");
            return;
        }
        for (j = i + 1; j < n; j++) {
            factor = A[j][i] / A[i][i];
            for (k = i; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    for (i = n - 1; i >= 0; i--) {
        x[i] = b[i] / A[i][i];
        for (j = i - 1; j >= 0; j--) {
            b[j] -= A[j][i] * x[i];
        }
    }
}