#ifndef FUNC_H
#define FUNC_H

void apply_elimination(void (*elimination_func)(double**, double*, int, int, int, int), double** A, double* x, int cols, int rows, int elems, int rank, char* title);
void gaussian_elimination_par(double** A, double* x, int cols, int rows, int elems, int rank);
void gaussian_elimination_seq(double** A, double* x, int cols, int rows, int elems, int rank);
void gauss_seidel_par(double** A, double* x, int cols, int rows, int elems, int rank);
void gauss_seidel_seq(double** A, double* x, int cols, int rows, int elems, int rank);

#endif
