#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <mpi.h>

#include "reader.h"
#include "func.h"

int main(int argc, char** argv) {
    FILE* matrix;
    char* filename = "data/matrix.txt";
    char* filename_vec = "data/vector.txt";
    int rows, cols, elems;
    double** A;
    double* b;
    double* x;

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) 
    {
        matrix = fopen(filename, "r");
        read_matrix_sizes(matrix, filename, &rows, &cols, &elems, size);
    }

    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&elems, 1, MPI_INT, 0, MPI_COMM_WORLD);

    create_matrix(&A, &x, &b, rows, cols);

    if (rank == 0) {
        read_matrix_values(matrix, filename, A, rows);

        FILE* vec = fopen(filename_vec, "r");
        read_vector_values(vec, filename_vec, b);
    }

    // printf("rows = %d\ncols = %d\nelems = %d\n", rows, cols, elems);
    apply_elimination(gaussian_elimination_seq, A, x, b, cols, rows, elems, rank, "Gaussian elimination sequential");
    // apply_elimination(gauss_seidel_seq, A, x, b, cols, rows, elems, rank, "Gaussian siedel sequential");
    apply_elimination(gaussian_elimination_par, A, x, b, cols, rows, elems, rank, "Gaussian elimination paraller");
    // apply_elimination(gauss_seidel_par, A, x, b, cols, rows, elems, rank, "Gaussian siedel parallel");

    for (int i = 0; i < rows; i++) {
        free(A[i]);
    }
    free(A);
    free(x);

    MPI_Finalize();

    return 0;
}
