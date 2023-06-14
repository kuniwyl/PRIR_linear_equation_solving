#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include "func.h"
#include "reader.h"

void apply_elimination(void (*elimination_func)(double**, double*, double*, double*, int, int, int, int, int), double** A, double* x, double* b, double* x_temp, int cols, int rows, int elems, int rank, int size, char* title) {
    clock_t start, end;
    double cpu_time_used;

    if (rank == 0) {
        printf("%s:\n", title);
        start = clock();
    }

    elimination_func(A, x, b, x_temp, cols, rows, elems, rank, size);

    if (rank == 0) {
        end = clock();
        cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
        printf("Czas wykonania: %f sekund\n", cpu_time_used);
        print_solution(x, rows);
        printf("\n");
    }
    return;
}

void gaussian_elimination_par(double** A, double* x, double* b, double* x_temp, int cols, int rows, int elems, int rank, int size) {
    int i, j;
    double factor;

    for (i = 0; i < rows - 1; i++) {
        for (j = i + 1; j < rows; j++) {
            if (rank == 0) {
                factor = A[j][i] / A[i][i]; 
            }
            MPI_Bcast(&factor, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            for(int k = i + rank; k < rows; k += size) {
                A[j][k] -= factor * A[i][k];
            }

            if (rank == 0) {
                b[j] -= factor * b[i];
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    for (i = rows - 1; i >= 0; i--) {
        if (rank == 0) {
            x[i] = b[i] / A[i][i];
        }
        MPI_Barrier(MPI_COMM_WORLD);

        for (j = i - 1 - rank; j >= 0; j -= size) {
            b[j] -= A[j][i] * x[i];
        }
    }
}

void gaussian_elimination_seq(double** A, double* x, double *b, double* x_temp, int cols, int rows, int elems, int rank, int size) {
    if (rank == 0) {
        int i, j, k;
        double factor;

        for (i = 0; i < rows - 1; i++) {
            for (j = i + 1; j < rows; j++) {
        
                factor = A[j][i] / A[i][i]; 
        
                for (k = i; k < rows; k++) {
                    A[j][k] -= factor * A[i][k];
                }
                b[j] -= factor * b[i];
            }
        }

        for (i = rows - 1; i >= 0; i--) {
            x[i] = b[i] / A[i][i];
            for (j = i - 1; j >= 0; j--) {
                b[j] -= A[j][i] * x[i];
            }
        }

    } 
    return;
}

void gauss_seidel_par(double** A, double* x, double* b, double* x_temp, int cols, int rows, int elems, int rank, int size) {
    double sigma_sum;
    double sigma;
    int i, j;
    double diff;
    double epsilon = 1e-6;

    for (;;) {
        for (i = 0; i < rows; i++) {
            sigma = 0.0;
            for (j = rank; j < rows; j += size) {
                if (j != i) {
                    sigma += A[i][j] * x[j];
                }
            }
            MPI_Allreduce(&sigma, &sigma_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);    
            if(rank == 0)
                x[i] = (b[i] - sigma_sum) / A[i][i];
        }

        diff = 0.0;
        for (i = rank; i < rows; i += size) {
            diff += fabs(x[i] - x_temp[i]);
            x_temp[i] = x[i];
        }
        double global_diff;
        MPI_Allreduce(&diff, &global_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (global_diff < epsilon) {
            return;
        }
    }
}

void gauss_seidel_seq(double** A, double* x, double* b, double* x_temp, int cols, int rows, int elems, int rank, int size) {
    if (rank == 0) {
        int i, j;
        double sigma;
        double diff;
        double epsilon = 1e-6;

        for (int k = 0; k < 10; k++) {
            for (i = 0; i < rows; i++) {
                sigma = 0.0;
                for (j = 0; j < rows; j++) {
                    if (j != i) {
                        sigma += A[i][j] * x[j];
                    }
                }
                x[i] = (b[i] - sigma) / A[i][i];
            }

            diff = 0.0;
            for (i = 0; i < rows; i++) {
                diff += fabs(x[i] - x_temp[i]);
                x_temp[i] = x[i];
            }

            if (diff < epsilon) {
                return;
            }
        }

    }
}
