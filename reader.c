#include <stdio.h>
#include <stdlib.h>

#include "reader.h"

void print_solution(double* x, int rows) {
    printf("Rozwiązanie:\n");
    for (int i = 0; i < rows; i++) {
        printf("x%d = %f\n", i + 1, x[i]);
    }
}

void print_matrix(double** A, int rows, int cols) {
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
}

double** copy_2d_array(double** original, int rows, int cols) {
    double** copy = (double**)malloc(sizeof(double*) * rows);
    if (copy == NULL) {
        printf("Memory allocation failed.\n");
        return NULL;
    }

    for (int i = 0; i < rows; i++) {
        copy[i] = (double*)malloc(sizeof(double) * cols);
        if (copy[i] == NULL) {
            printf("Memory allocation failed.\n");

            for (int j = 0; j < i; j++) {
                free(copy[j]);
            }
            free(copy);

            return NULL;
        }

        for (int j = 0; j < cols; j++) {
            copy[i][j] = original[i][j];
        }
    }

    return copy;
}

double* copy_from_matrix(double** original, int rows, int index) {
    double* b = malloc(sizeof(double) * rows);
    for(int i = 0; i < rows; i++) {
        b[i] = original[i][index];
    }
    return b;
}

double* copy_vector(double* original, int size) {
    double* copy = (double*)malloc(sizeof(double) * size);
    if (copy == NULL) {
        printf("Memory allocation failed.\n");
        return NULL;
    }

    for (int i = 0; i < size; i++) {
        copy[i] = original[i];
    }

    return copy;
}

void create_matrix(double*** A, double** x, double** b, int rows, int cols){
    *A = (double**)malloc(sizeof(double*) * rows);
    for (int i = 0; i < rows; i++) {
        (*A)[i] = (double*)malloc(sizeof(double) * cols);
    }
    *x = malloc(sizeof(double) * rows);
    *b = malloc(sizeof(double) * rows);
}

int read_vector_values(FILE* file, char* filename, double* b) {
    if (file == NULL) {
        printf("Failed to open file: %s\n", filename);
        return -1;
    }

    int rows;
    if (fscanf(file, "%d", &rows) != 1) {
        printf("Failed to read vector size.\n");
        fclose(file);
        return -1;
    }

    for (int i = 0; i < rows; i++) {
        if (fscanf(file, "%lf", &b[i]) != 1) {
            printf("Failed to read vector element at index %d.\n", i);
            fclose(file);
            return -1;
        }
    }

    fclose(file);
    return 1;
}


int read_matrix_sizes(FILE* file, char* filename, int* rows, int* cols, int* elems, int size) {
    if (file == NULL) {
        printf("Error opening file: %s\n", filename);
        return -1;
    }

    if (fscanf(file, "%d", rows) != 1 || fscanf(file, "%d", cols) != 1) {
        printf("Invalid file format: %s\n", filename);
        fclose(file);
        return -1;
    }

    int row_count = *rows;
    if (row_count % size != 0) {
        *elems = row_count / size + 1;
        *cols = size * *elems;
    } else {
        *elems = row_count / size;
        *cols = row_count;
    }

    return 1;
}


int read_matrix_values(FILE* file, char* filename, double** A, int rows) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            if (fscanf(file, "%lf", &A[i][j]) != 1) {
                printf("Invalid file format: %s\n", filename);
                fclose(file);
                return -1;
            }
        }
    }
    fclose(file);
    return 1;
}

void prepare(double** A, double** A_values, int rows, double* b, double* b_values) {
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < rows; j++) {
            A[i][j] = A_values[i][j];
        }
        b[i] = b_values[i];
    }
}