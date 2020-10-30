#include "my_solver.h"

#define SOLNUMS 2
#define MATCOLS SOLNUMS
#define MATROWS SOLNUMS
#define TOLERANCE 0.0000001

void test6(double* x, double* func) {
	func[0] = 3 * x[0] * x[0] - 2 * x[1] * x[1] - 1;
	func[1] = x[0] * x[0] - 2 * x[0] + x[1] * x[1] + 2 * x[1] - 8;
}

void print6(FILE* fp_w, double* x, int info) {
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
		test6(x, func);

		printf("|F(X)| : %lf %lf\n\n", fabs(func[0]), fabs(func[1]));
		fprintf(fp_w, "|F(X)| : %lf %lf\n\n", fabs(func[0]), fabs(func[1]));
	}
}

void fcn3_6(int* n, double* x, double* fvec, double* fjac, int* ldfjac, int* iflag) {
	// origin function F(x)
	if (*iflag == 1) {
		/********************************/

		fvec[0] = 3 * x[0] * x[0] - 2 * x[1] * x[1] - 1;
		fvec[1] = x[0] * x[0] - 2 * x[0] + x[1] * x[1] + 2 * x[1] - 8;
	
		/********************************/
	}
	// Jacobi matrix J(x)
	else if (*iflag == 2) {
		/********************************/
		fjac[0] = 6 * x[0];						fjac[2] = -4 * x[1];
		fjac[1] = 2 * x[0] - 2;					fjac[3] = 2 * x[1] + 2;
		
		/********************************/
	}
}

void practice3_6(void)
{
	int n = SOLNUMS;
	double x[SOLNUMS] = { -4.0, 5.0 };	//need to initilize x0
	double fvec[SOLNUMS], fjac[MATCOLS * MATROWS];
	int ldfjac = SOLNUMS;
	double tol = TOLERANCE;
	int info;
	double wa[(SOLNUMS * (SOLNUMS + 13)) / 2];
	int lwa = (SOLNUMS * (SOLNUMS + 13)) / 2;

	FILE* fp_w = fopen("roots_3-6.txt", "w");
	if (fp_w == NULL) {
		printf("%s file open error...\n", "roots_3-6.txt");
		return;
	}

	/********************************/

	printf("***** Problem 3-6 *****\n");
	printf("Initial value of x, y = (-4.0, 5.0)\n\n");
	fprintf(fp_w, "Initial value of x = (-4.0, 5.0)\n\n");

	hybrj1_(fcn3_6, &n, x, fvec, fjac, &ldfjac, &tol, &info, wa, &lwa);

	print6(fp_w, x, info);

	x[0] = 2.0; x[1] = 5.0;
	printf("Initial value of x, y = (2.0, 5.0)\n\n");
	fprintf(fp_w, "Initial value of x = (2.0, 5.0)\n\n");

	hybrj1_(fcn3_6, &n, x, fvec, fjac, &ldfjac, &tol, &info, wa, &lwa);

	print6(fp_w, x, info);

	x[0] = -4.0; x[1] = -5.0;
	printf("Initial value of x, y = (-4.0, -5.0)\n\n");
	fprintf(fp_w, "Initial value of x = (-4.0, -5.0)\n\n");

	hybrj1_(fcn3_6, &n, x, fvec, fjac, &ldfjac, &tol, &info, wa, &lwa);

	print6(fp_w, x, info);

	x[0] = 4.0; x[1] =-3.0;
	printf("Initial value of x, y = (4.0, -3.0)\n\n");
	fprintf(fp_w, "Initial value of x = (4.0, -3.0)\n\n");

	hybrj1_(fcn3_6, &n, x, fvec, fjac, &ldfjac, &tol, &info, wa, &lwa);

	print6(fp_w, x, info);

	/********************************/

	fclose(fp_w);
}