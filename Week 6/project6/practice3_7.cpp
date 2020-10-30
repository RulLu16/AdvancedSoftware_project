#include "my_solver.h"

#define SOLNUMS 2
#define MATCOLS SOLNUMS
#define MATROWS SOLNUMS
#define TOLERANCE 0.0000001

void test7(double* x, double* func) {
	func[0] = 2 * x[0] * x[0] * x[0] - 12 * x[0] - x[1] - 1;
	func[1] = 3 * x[1] * x[1] - 6 * x[1] - x[0] - 3;
}

void print7(FILE* fp_w, double* x, int info) {
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
		test7(x, func);

		printf("|F(X)| : %lf %lf\n\n", fabs(func[0]), fabs(func[1]));
		fprintf(fp_w, "|F(X)| : %lf %lf\n\n", fabs(func[0]), fabs(func[1]));
	}
}

void fcn3_7(int* n, double* x, double* fvec, int* iflag)
{
	/********************************/

	fvec[0] = 2 * x[0] * x[0] * x[0] - 12 * x[0] - x[1] - 1;
	fvec[1] = 3 * x[1] * x[1] - 6 * x[1] - x[0] - 3;

	/********************************/
}

void practice3_7(void)
{
	int n = SOLNUMS;
	double x[SOLNUMS] = { -4.0, 5.0 };	//need to initilize x0
	double fvec[SOLNUMS];
	double tol = TOLERANCE;
	int info;
	double wa[(SOLNUMS * (3 * SOLNUMS + 13)) / 2];
	int lwa = (SOLNUMS * (3 * SOLNUMS + 13)) / 2;

	FILE* fp_w = fopen("roots_3-7.txt", "w");
	if (fp_w == NULL) {
		printf("%s file open error...\n", "roots_3-7.txt");
		return;
	}

	/********************************/

	printf("***** Problem 3-7 *****\n");
	printf("Initial value of x, y = (-4.0, 5.0)\n\n");
	fprintf(fp_w, "Initial value of x = (-4.0, 5.0)\n\n");

	hybrd1_(fcn3_7, &n, x, fvec, &tol, &info, wa, &lwa);

	print7(fp_w, x, info);

	x[0] = 2.0; x[1] = 5.0;
	printf("Initial value of x, y = (2.0, 5.0)\n\n");
	fprintf(fp_w, "Initial value of x = (2.0, 5.0)\n\n");

	hybrd1_(fcn3_7, &n, x, fvec, &tol, &info, wa, &lwa);

	print7(fp_w, x, info);

	x[0] = -4.0; x[1] = -5.0;
	printf("Initial value of x, y = (-4.0, -5.0)\n\n");
	fprintf(fp_w, "Initial value of x = (-4.0, -5.0)\n\n");

	hybrd1_(fcn3_7, &n, x, fvec, &tol, &info, wa, &lwa);

	print7(fp_w, x, info);

	x[0] = 4.0; x[1] = -3.0;
	printf("Initial value of x, y = (4.0, -3.0)\n\n");
	fprintf(fp_w, "Initial value of x = (4.0, -3.0)\n\n");

	hybrd1_(fcn3_7, &n, x, fvec, &tol, &info, wa, &lwa);

	print7(fp_w, x, info);

	x[0] = 1.0; x[1] = -5.0;
	printf("Initial value of x, y = (1.0, -5.0)\n\n");
	fprintf(fp_w, "Initial value of x = (1.0, -5.0)\n\n");

	hybrd1_(fcn3_7, &n, x, fvec, &tol, &info, wa, &lwa);

	print7(fp_w, x, info);

	x[0] = 0.0; x[1] = 1.0;
	printf("Initial value of x, y = (0.0, 1.0)\n\n");
	fprintf(fp_w, "Initial value of x = (0.0, 1.0)\n\n");

	hybrd1_(fcn3_7, &n, x, fvec, &tol, &info, wa, &lwa);

	print7(fp_w, x, info);

	/********************************/

	fclose(fp_w);
}