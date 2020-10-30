#include "my_solver.h"

#define SOLNUMS 3
#define MATCOLS SOLNUMS
#define MATROWS SOLNUMS
#define TOLERANCE 0.0000001

void test5(double* x, double* func) {
	func[0] = 10 * x[0] - 2 * x[1] * x[1] + x[1] - 2 * x[2] - 5;
	func[1] = 8 * x[1] * x[1] + 4 * x[2] * x[2] - 9;
	func[2] = 8 * x[1] * x[2] + 4;
}

void print5(FILE* fp_w, double* x, int info) {
	printf("Info value is %d\n", info);
	printf("%lf %lf %lf\n\n", x[0], x[1], x[2]);
	fprintf(fp_w, "Info value is %d\n", info);
	fprintf(fp_w, "%lf %lf %lf\n\n", x[0], x[1], x[2]);

	if (info != 1) {
		printf("Error occurs. That solution is not proper.\n\n");
		fprintf(fp_w, "Error occurs. That solution is not proper.\n\n");
	}
	else {
		double func[SOLNUMS] = { 1, };
		test5(x, func);

		printf("|F(X)| : %lf %lf %lf\n\n", fabs(func[0]), fabs(func[1]), fabs(func[2]));
		fprintf(fp_w, "|F(X)| : %lf %lf %lf\n\n", fabs(func[0]), fabs(func[1]), fabs(func[2]));
	}
}

void fcn3_5(int* n, double* x, double* fvec, int* iflag)
{
	/********************************/

	fvec[0] = 10 * x[0] - 2 * x[1] * x[1] + x[1] - 2 * x[2] - 5;
	fvec[1] = 8 * x[1] * x[1] + 4 * x[2] * x[2] - 9;
	fvec[2] = 8 * x[1] * x[2] + 4;

	/********************************/
}

void practice3_5(void)
{
	int n = SOLNUMS;
	double x[SOLNUMS] = { 1.0, -1.0, 1.0 };	//need to initilize x0
	double fvec[SOLNUMS];
	double tol = TOLERANCE;
	int info;
	double wa[(SOLNUMS * (3 * SOLNUMS + 13)) / 2];
	int lwa = (SOLNUMS * (3 * SOLNUMS + 13)) / 2;

	FILE* fp_w = fopen("roots_3-5.txt", "w");
	if (fp_w == NULL) {
		printf("%s file open error...\n", "roots_3-5.txt");
		return;
	}

	/********************************/
	printf("***** Problem 3-5 *****\n");
	printf("Initial value of x = (0.1, 1.2, 2.5)\n\n");
	fprintf(fp_w, "Initial value of x = (0.1, 1.2, 2.5)\n\n");

	hybrd1_(fcn3_5, &n, x, fvec, &tol, &info, wa, &lwa);

	print5(fp_w, x, info);

	x[0] = 1.0; x[1] = 0.0; x[2] = 1.0;
	printf("Initial value of x = (1.0, 0.0, 1.0)\n\n");
	fprintf(fp_w, "Initial value of x = (1.0, 0.0, 1.0)\n\n");

	hybrd1_(fcn3_5, &n, x, fvec, &tol, &info, wa, &lwa);

	print5(fp_w, x, info);

	/********************************/

	fclose(fp_w);
}