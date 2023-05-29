#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include "func.h"
#include "reader.h"

#define epsilon 0.000001

void apply_elimination(void (*elimination_func)(double**, double*, int, int, int, int), double** A, double* x, int cols, int rows, int elems, int rank, char* title) {
    clock_t start, end;
    double cpu_time_used;
    double **A_copy;
    double *x_copy;

    if (rank == 0) {
        A_copy = copy_2d_array(A, rows, cols);
        x_copy = copy_vector(x, rows);

        printf("%s:\n", title);
        start = clock();
    }

    elimination_func(A_copy, x_copy, cols, rows, elems, rank);

    if (rank == 0) {
        end = clock();
        cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
        printf("Czas wykonania: %f sekund\n", cpu_time_used);
        print_solution(x_copy, rows);
        printf("\n");
    }
    return;
}

void gaussian_elimination_par(double** A, double* x, int cols, int rows, int elems, int rank) {
    int i, j, k;
    double factor;
    double* up = malloc(sizeof(double) * cols);
    double* down = malloc(sizeof(double) * cols);

    for (i = 0; i < rows - 1; i++) {
        if (A[i][i] == 0) {
            printf("Macierz nie ma unikalnego rozwiązania.\n");
            return;
        }
        for (j = i + 1; j < rows; j++) {
            if (rank == 0) {
                factor = A[j][i] / A[i][i];
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&factor, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            MPI_Scatter(A[i], elems, MPI_DOUBLE, up, elems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Scatter(A[j], elems, MPI_DOUBLE, down, elems, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            for (int q = 0; q < elems; q++) {
                down[q] -= factor * up[q]; 
            }

            MPI_Gather(up, elems, MPI_DOUBLE, A[i], elems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(down, elems, MPI_DOUBLE, A[j], elems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

    if (rank == 0) {
        for (i = rows - 1; i >= 0; i--) {
            x[i] = A[i][rows] / A[i][i];
            for (j = i - 1; j >= 0; j--) {
                A[j][rows] -= A[j][i] * x[i];
            }
        }
    }

    free(up);
    free(down);
}

void gaussian_elimination_seq(double** A, double* x, int cols, int rows, int elems, int rank) {
    if (rank == 0) {
        int i, j, k;
        double factor;
        for (i = 0; i < rows - 1; i++) {
            if (A[i][i] == 0) {
                printf("Macierz nie ma unikalnego rozwiązania.\n");
                return;
            }
            for (j = i + 1; j < rows; j++) {
                factor = A[j][i] / A[i][i]; 
                for (k = i; k < rows + 1; k++) {
                    A[j][k] -= factor * A[i][k];
                }
            }
        }

        for (i = rows - 1; i >= 0; i--) {
            x[i] = A[i][rows] / A[i][i];
            for (j = i - 1; j >= 0; j--) {
                A[j][rows] -= A[j][i] * x[i];
            }
        }
    }
}

void gauss_seidel_par(double** A, double* x, int cols, int rows, int elems, int rank) {

    // double** A_copy = copy_2d_array(A, rows, rows);
    // double* B = copy_from_matrix(A, rows, rows);
    double* x_temp = malloc(sizeof(double) * cols);
    double* x_up = malloc(sizeof(double) * cols);
    double* up = malloc(sizeof(double) * cols);
    double* sigma = malloc(sizeof(double) * cols);
    double* sigma_elemes = malloc(sizeof(double) * cols);
    double* sigma_sums = malloc(sizeof(double) * cols);
    double sigma_sum;
    int i, j, k;
    double diff;

    for (;;) {
        for (i = 0; i < rows; i++) {
            
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Scatter(A[i], elems, MPI_DOUBLE, up, elems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Scatter(x, elems, MPI_DOUBLE, x_up, elems, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            for(int q = 0; q < elems; q++) {
                sigma_elemes[q] = up[q] * x_up[q];
            }
            
            MPI_Gather(sigma_elemes, elems, MPI_DOUBLE, sigma, elems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            
            if (rank == 0) 
            {
                sigma_sum = 0.0;
                for(int j = 0; j < rows + 1; j++) {
                    sigma_sum += sigma[j];
                }
                sigma_sum -= A[i][i] * x[i];
                sigma_sum -= A[i][rows] * x[rows];
                x_temp[i] = (A[i][rows] - sigma_sum) / A[i][i];
            }
        }

        if (rank == 0) 
        {
            diff = 0.0;
            for (i = 0; i < rows; i++) {
                diff += fabs(x_temp[i] - x[i]);
                x[i] = x_temp[i];
            }
            if (diff < epsilon) {
                return;
            }
        }
    }
}

void gauss_seidel_seq(double** A, double* x, int cols, int rows, int elems, int rank) {
    if (rank == 0) {
        double* x_temp = malloc(sizeof(double) * cols);
        int i, j, k;
        double sigma;
        double diff;

        for (;;) {
            for (i = 0; i < rows; i++) {
                sigma = 0.0;
                for (j = 0; j < rows; j++) {
                    if (j != i) {
                        sigma += A[i][j] * x[j];
                    }
                }
                x_temp[i] = (A[i][rows] - sigma) / A[i][i];
            }

            diff = 0.0;
            for (i = 0; i < rows; i++) {
                diff += fabs(x_temp[i] - x[i]);
                x[i] = x_temp[i];
            }
            if (diff < epsilon) {
                return;
            }
        }
    }
}
