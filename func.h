#ifndef FUNC_H
#define FUNC_H

void apply_elimination(void (*elimination_func)(double**, double*, double*, double*, int, int, int, int, int), double** A, double* x, double* b, double* x_temp, int cols, int rows, int elems, int rank, int size, char* title);
void gaussian_elimination_par(double** A, double* x, double* b, double* x_temp, int cols, int rows, int elems, int rank, int size);
void gaussian_elimination_seq(double** A, double* x, double* b, double* x_temp, int cols, int rows, int elems, int rank, int size);
void gauss_seidel_par(double** A, double* x, double* b, double* x_temp, int cols, int rows, int elems, int rank, int size);
void gauss_seidel_seq(double** A, double* x, double* b, double* x_temp, int cols, int rows, int elems, int rank, int size);

#endif
