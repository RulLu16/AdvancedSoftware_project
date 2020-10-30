#include "my_solver.h"

#define SOLNUMS 3
#define MATCOLS SOLNUMS
#define MATROWS SOLNUMS
#define TOLERANCE 0.0000001

void test3(double* x, double* func) {

	func[0] = exp(2 * x[0]) - x[1] + 4;
	func[1] = x[1] - x[2] * x[2] - 1;
	func[2] = x[2] - sin(x[0]);
}

void print3(FILE* fp_w, double* x, int info) {
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
		test3(x, func);

		printf("|F(X)| : %lf %lf %lf\n\n", fabs(func[0]), fabs(func[1]), fabs(func[2]));
		fprintf(fp_w, "|F(X)| : %lf %lf %lf\n\n", fabs(func[0]), fabs(func[1]), fabs(func[2]));
	}
}

void fcn3_3(int* n, double* x, double* fvec, double* fjac, int* ldfjac, int* iflag)
{
	// origin function F(x)
	if (*iflag == 1) {
		/********************************/
		fvec[0] = exp(2 * x[0]) - x[1] + 4;
		fvec[1] = x[1] - x[2] * x[2] - 1;
		fvec[2] = x[2] - sin(x[0]);

		/********************************/
	}
	// Jacobi matrix J(x)
	else if (*iflag == 2) {
		/********************************/

		fjac[0] = 2 * exp(1 * x[0]);	fjac[3] = -1;	fjac[6] = 0;
		fjac[1] = 0;					fjac[4] = 1;	fjac[7] = -2 * x[2];
		fjac[2] = -cos(x[0]);			fjac[5] = 0;	fjac[8] = 1;

		/********************************/
	}
}

void practice3_3(void) {
	int n = SOLNUMS;
	double x[SOLNUMS] = { 0.0, 0.0, 0.0 };	//need to initilize x0
	double fvec[SOLNUMS], fjac[MATCOLS * MATROWS];
	int ldfjac = SOLNUMS;
	double tol = TOLERANCE;
	int info;
	double wa[(SOLNUMS * (SOLNUMS + 13)) / 2];
	int lwa = (SOLNUMS * (SOLNUMS + 13)) / 2;

	FILE* fp_w = fopen("roots_3-3.txt", "w");
	if (fp_w == NULL) {
		printf("%s file open error...\n", "roots_3-3.txt");
		return;
	}

	/********************************/
	printf("***** Problem 3-3 *****\n");
	printf("Initial value of x = (0.0, 0.0, 0.0)\n\n");
	fprintf(fp_w, "Initial value of x = (0.0, 0.0, 0.0)\n\n");

	hybrj1_(fcn3_3, &n, x, fvec, fjac, &ldfjac, &tol, &info, wa, &lwa);

	print3(fp_w, x, info);

	x[0] = 1.55; x[1] = 1.39; x[2] = 1.10;
	printf("Initial value of x = (1.55, 1.39, 1.10)\n\n");
	fprintf(fp_w, "Initial value of x = (1.55, 1.39, 1.10)\n\n");

	hybrj1_(fcn3_3, &n, x, fvec, fjac, &ldfjac, &tol, &info, wa, &lwa);

	print3(fp_w, x, info);

	/********************************/

	fclose(fp_w);
}
