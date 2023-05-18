#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gauss_seidel.h"
#include "gaussian_elimination.h"

#define N 3

void print_solution(double x[N]) {
    printf("Rozwiązanie:\n");
    for (int i = 0; i < N; i++) {
        printf("x%d = %f\n", i + 1, x[i]);
    }
}

int main() {
    double** A;
    double* b;
    double* x;

    int i;
    A = malloc(sizeof(double*) * N);
    for(i = 0; i < N; i++) 
    {
        A[i] = malloc(sizeof(double) * N);
    }
    b = malloc(sizeof(double) * N);
    x = malloc(sizeof(double) * N);

    A[0][0] = 4;
    A[0][1] = 1;
    A[0][2] = -1;
    A[1][0] = 2;
    A[1][1] = 7;
    A[1][2] = 1;
    A[2][0] = 1;
    A[2][1] = 2;
    A[2][2] = 6;

    b[0] = 3;
    b[1] = 19;
    b[2] = 23;

    clock_t start, end;
    double cpu_time_used;


    printf("Metoda Gaussa-Seidla:\n");
    start = clock();
    x[0] = 0;
    x[1] = 0;
    x[2] = 0;
    gauss_seidel(A, b, x, 0.000001, N, 100);
    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Czas wykonania: %f sekund\n", cpu_time_used);
    print_solution(x);


    printf("\nMetoda eliminacji Gaussa:\n");
    start = clock();
    gaussian_elimination(A, b, x, N);
    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Czas wykonania: %f sekund\n", cpu_time_used);
    print_solution(x);

    return 0;
}
