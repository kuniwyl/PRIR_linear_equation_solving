#ifndef READER_H
#define READER_H

int read_matrix_sizes(FILE* file, char* filename, int *rows, int *cols, int *elems, int *size);
int read_matrix_values(FILE* file, char* filename, double** A, int rows);
int read_vector_values(FILE* file, char* filename, double** A, int size);
void create_matrix(double*** A, double** x, int rows, int cols);
double** copy_2d_array(double** original, int rows, int cols);
double* copy_from_matrix(double** original, int rows, int index);
double* copy_vector(double* original, int size);
void print_solution(double* x, int rows);
void print_matrix(double** A, int rows, int cols);

#endif
