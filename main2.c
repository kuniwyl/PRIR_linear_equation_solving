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
    double* x_temp;
    
    double** A_values;
    double* b_values;

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

    // SHARED MEMORY
    MPI_Offset shared_mem_size = sizeof(double) * (rows * cols + 3 * rows);
    double* data;
    MPI_Win win;
    MPI_Win_allocate_shared(shared_mem_size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &data, &win);
    
    double* base_addr;
    MPI_Aint win_size;
    int disp_unit;
    MPI_Win_shared_query(win, 0, &win_size, &disp_unit, &base_addr);

    A = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        A[i] = (double*)(base_addr + i * cols);
    }

    b = (double*)(base_addr + rows * cols);
    x = (double*)(base_addr + rows * cols + rows);
    x_temp = (double*)(base_addr + rows * cols + rows + rows);

    if (rank == 0) {
        read_matrix_values(matrix, filename, A, rows);

        A_values = copy_2d_array(A, rows, cols);

        FILE* vec = fopen(filename_vec, "r");
        read_vector_values(vec, filename_vec, b);

        b_values = copy_vector(b, rows);

        for(int i = 0; i < rows; i++) {
            x[i] = 0;
            x_temp[i] = 0;
        }
    }

    
    apply_elimination(gaussian_elimination_seq, A, x, b, x_temp, cols, rows, elems, rank, size, "Gaussian elimination sequential");

    if (rank == 0) {
        prepare(A, A_values, rows, b, b_values);
        for(int i = 0; i < rows; i++) {
            x[i] = 0;
            x_temp[i] = 0;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    apply_elimination(gauss_seidel_seq, A, x, b, x_temp, cols, rows, elems, rank, size, "Gaussian siedel sequential");
    
    if (rank == 0) {
        prepare(A, A_values, rows, b, b_values);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    apply_elimination(gaussian_elimination_par, A, x, b, x_temp, cols, rows, elems, rank, size, "Gaussian elimination paraller");
    
    if (rank == 0) {
        prepare(A, A_values, rows, b, b_values);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    apply_elimination(gauss_seidel_par, A, x, b, x_temp, cols, rows, elems, rank, size, "Gaussian siedel parallel");

    MPI_Win_free(&win);
    MPI_Finalize();

    return 0;
}
