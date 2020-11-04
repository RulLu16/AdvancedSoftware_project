#define _CRT_SECURE_NO_WARNINGS
#include "my_solver.h"

#define SOLNUMS 4
#define MATCOLS SOLNUMS
#define MATROWS SOLNUMS
#define TOLERANCE 0.0000001

void test2(double* x, double* func) {
	func[0] = x[0] + 10 * x[1] + 9;
	func[1] = sqrt(5.0) * (x[2] - x[3]) - 2 * sqrt(5.0);
	func[2] = (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) - 9;
	func[3] = sqrt(10.0) * (x[0] - x[3]) * (x[0] - x[3]) - 2 * sqrt(10.0);
}

void print2(FILE* fp_w, double* x, int info) {
	printf("Info value is %d\n", info);
	printf("%lf %lf %lf %lf\n\n", x[0], x[1], x[2], x[3]);
	fprintf(fp_w, "Info value is %d\n", info);
	fprintf(fp_w, "%lf %lf %lf %lf\n\n", x[0], x[1], x[2], x[3]);

	if (info != 1) {
		printf("Error occurs. That solution is not proper.\n\n");
		fprintf(fp_w, "Error occurs. That solution is not proper.\n\n");
	}
	else {
		double func[SOLNUMS] = { 1, };
		test2(x, func);

		printf("|F(X)| : %lf %lf %lf %lf\n\n", fabs(func[0]), fabs(func[1]), fabs(func[2]), fabs(func[3]));
		fprintf(fp_w, "|F(X)| : %lf %lf %lf %lf\n\n", fabs(func[0]), fabs(func[1]), fabs(func[2]), fabs(func[3]));
	}
}

void fcn3_2(int* n, double* x, double* fvec, int* iflag)
{
	fvec[0] = x[0] + 10 * x[1] + 9;
	fvec[1] = sqrt(5.0) * (x[2] - x[3]) - 2 * sqrt(5.0);
	fvec[2] = (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) - 9;
	fvec[3] = sqrt(10.0) * (x[0] - x[3]) * (x[0] - x[3]) - 2 * sqrt(10.0);
}

void practice3_2(void) {
	int n = SOLNUMS;
	double x[SOLNUMS] = { 0.9, -0.9, 1.25, -1.25 };	//need to initilize x0
	double fvec[SOLNUMS];
	double tol = TOLERANCE;
	int info;
	double wa[(SOLNUMS * (3 * SOLNUMS + 13)) / 2];
	int lwa = (SOLNUMS * (3 * SOLNUMS + 13)) / 2;

	FILE* fp_w = fopen("roots_found_3-2.txt", "w");
	if (fp_w == NULL) {
		printf("%s file open error...\n", "roots_found_3-2.txt");
		return;
	}

	/********************************/
	printf("\n***** Problem 3-2 *****\n");
	printf("Initial value of x = (0.9, -0.9, 1.25, -1.25)\n\n");
	fprintf(fp_w, "Initial value of x = (0.9, -0.9, 1.25, -1.25)\n\n");

	hybrd1_(fcn3_2, &n, x, fvec, &tol, &info, wa, &lwa);

	print2(fp_w, x, info);

	/********************************/

	fclose(fp_w);
}