#define _CRT_SECURE_NO_WARNINGS
#include "my_solver.h"

#define SOLNUMS 2
#define MATCOLS SOLNUMS
#define MATROWS SOLNUMS
#define TOLERANCE 0.0000001

void test3(double* x, double* func) {
	func[0] = ((sin(x[0] * x[1] + M_PI / 6) + sqrt(x[0] * x[0] * x[1] * x[1] + 1)) / cos(x[0] - x[1])) + 2.8;
	func[1] = ((x[0] * exp(x[0] * x[1] + M_PI / 6) - sin(x[0] - x[1])) / sqrt(x[0] * x[0] * x[1] * x[1] + 1)) - 1.66;
}

void print3(FILE* fp_w, double* x, int info) {
	printf("Info value is %d\n", info);
	printf("%lf %lf\n\n", x[0], x[1]);
	fprintf(fp_w, "Info value is %d\n", info);
	fprintf(fp_w, "%lf %lf\n\n", x[0], x[1]);

	if (info != 1) {
		printf("Error occurs. That solution is not proper.\n\n");
		fprintf(fp_w, "Error occurs. That solution is not proper.\n\n");
	}
	else {
		double func[SOLNUMS] = { 1, };
		test3(x, func);

		printf("|F(X)| : %lf %lf\n\n", fabs(func[0]), fabs(func[1]));
		fprintf(fp_w, "|F(X)| : %lf %lf\n\n", fabs(func[0]), fabs(func[1]));
	}
}

void fcn3_3(int* n, double* x, double* fvec, int* iflag)
{
	fvec[0] = ((sin(x[0] * x[1] + M_PI / 6) + sqrt(x[0] * x[0] * x[1] * x[1] + 1)) / cos(x[0] - x[1])) + 2.8;
	fvec[1] = ((x[0] * exp(x[0] * x[1] + M_PI / 6) - sin(x[0] - x[1])) / sqrt(x[0] * x[0] * x[1] * x[1] + 1)) - 1.66;
}

void practice3_3(void) {
	int n = SOLNUMS;
	double x[SOLNUMS] = { 20.0, 0.0 };	//need to initilize x0
	double fvec[SOLNUMS];
	double tol = TOLERANCE;
	int info;
	double wa[(SOLNUMS * (3 * SOLNUMS + 13)) / 2];
	int lwa = (SOLNUMS * (3 * SOLNUMS + 13)) / 2;

	FILE* fp_w = fopen("roots_found_3-3.txt", "w");
	if (fp_w == NULL) {
		printf("%s file open error...\n", "roots_found_3-3.txt");
		return;
	}

	/********************************/
	printf("\n***** Problem 3-3 *****\n");
	printf("Initial value of x = (20.0, 0.0)\n\n");
	fprintf(fp_w, "Initial value of x = (20.0, 0.0)\n\n");

	hybrd1_(fcn3_3, &n, x, fvec, &tol, &info, wa, &lwa);

	print3(fp_w, x, info);

	/********************************/

	fclose(fp_w);
}