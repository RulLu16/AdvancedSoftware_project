#define _CRT_SECURE_NO_WARNINGS
#include "my_solver.h"

#define index_2d(row, col, dim) row*dim+col

void practice3_4(void) {
    int i, ia, j, n;
    int* l;
    float* a; float* b; 
    float* x; float* s;

    /********************************/
    FILE* fp_r = fopen("linear_system_3-4.txt", "r");
    if (fp_r == NULL) {
        printf("%s file open error...\n", "linear_system_3-4.txt");
        return;
    }

    FILE* fp_w = fopen("solution_3-4.txt", "w");
    if (fp_w == NULL) {
        printf("%s file open error...\n", "solution_3-4.txt");
        return;
    }

    fscanf(fp_r, "%d", &n);

    a = (float*)malloc(sizeof(float) * (n * n));
    l = (int*)malloc(sizeof(int) * n);
    b = (float*)malloc(sizeof(float) * n);
    x = (float*)malloc(sizeof(float) * n);
    s = (float*)malloc(sizeof(float) * n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            fscanf(fp_r, "%f", &a[index_2d(j, i, n)]);
        }
    }

    for (i = 0; i < n; i++) {
        fscanf(fp_r, "%f", &b[i]);
    }

    ia = n;
 
    gespp_(&n, a, &ia, l, s);
    solve_(&n, a, &ia, l, b, x);
    printf("\n***** Solution *****\n");
    fprintf(fp_w, "%d\n", n);
    for (i = 0; i < n; i++) {
        printf("x[%d] = %f\n", i, x[i]);
        fprintf(fp_w, "%.6f\n", x[i]);
    }
    printf("\n");

    float* child = (float*)malloc(sizeof(float) * n);
    for (i = 0; i < n; i++) {
        child[i] = 0;
        for (j = 0; j < n; j++) {
            child[i] += a[index_2d(j, i, n)] * x[j];
        }
        child[i] -= b[i];
    }

    float up = 0; float down = 0;

    for (i = 0; i < n; i++) {
        up += child[i] * child[i];
        down += b[i] * b[i];
    }
    printf("Error rate: %f\n", sqrt(up) / sqrt(down));
    fprintf(fp_w, "%.6f\n", sqrt(up) / sqrt(down));
    /********************************/
}
