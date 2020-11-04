#define _CRT_SECURE_NO_WARNINGS
#include "my_solver.h"

#define SOLNUMS 4
#define MATCOLS SOLNUMS
#define MATROWS SOLNUMS
#define TOLERANCE 0.0000001

double C, b;
double p[4][3] = { 0, };
double t[4] = { 0, };
double tr[4] = { 0, };

void test1(double* x, double* func) {

	func[0] = (x[0] - p[0][0]) * (x[0] - p[0][0]) + (x[1] - p[0][1]) * (x[1] - p[0][1]) + (x[2] - p[0][2]) * (x[2] - p[0][2])
		- (C * (tr[0] + x[3] - t[0])) * (C * (tr[0] + x[3] - t[0]));
	func[1] = (x[0] - p[1][0]) * (x[0] - p[1][0]) + (x[1] - p[1][1]) * (x[1] - p[1][1]) + (x[2] - p[1][2]) * (x[2] - p[1][2])
		- (C * (tr[1] + x[3] - t[1])) * (C * (tr[1] + x[3] - t[1]));
	func[2] = (x[0] - p[2][0]) * (x[0] - p[2][0]) + (x[1] - p[2][1]) * (x[1] - p[2][1]) + (x[2] - p[2][2]) * (x[2] - p[2][2])
		- (C * (tr[2] + x[3] - t[2])) * (C * (tr[2] + x[3] - t[2]));
	func[3] = (x[0] - p[3][0]) * (x[0] - p[3][0]) + (x[1] - p[3][1]) * (x[1] - p[3][1]) + (x[2] - p[3][2]) * (x[2] - p[3][2])
		- (C * (tr[3] + x[3] - t[3])) * (C * (tr[3] + x[3] - t[3]));
}

void print1(FILE* fp_w, double* x, int info) {
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
		test1(x, func);

		printf("|F(X)| : %lf %lf %lf %lf\n\n", fabs(func[0]), fabs(func[1]), fabs(func[2]), fabs(func[3]));
		fprintf(fp_w, "|F(X)| : %lf %lf %lf %lf\n\n", fabs(func[0]), fabs(func[1]), fabs(func[2]), fabs(func[3]));
	}
}

void fcn3_1_hybrj1(int* n, double* x, double* fvec, double* fjac, int* ldfjac, int* iflag)
{
	// origin function F(x)
	if (*iflag == 1) {
		/********************************/
		fvec[0] = (x[0] - p[0][0]) * (x[0] - p[0][0]) + (x[1] - p[0][1]) * (x[1] - p[0][1]) + (x[2] - p[0][2]) * (x[2] - p[0][2])
			- (C * (tr[0] + x[3] - t[0])) * (C * (tr[0] + x[3] - t[0]));
		fvec[1] = (x[0] - p[1][0]) * (x[0] - p[1][0]) + (x[1] - p[1][1]) * (x[1] - p[1][1]) + (x[2] - p[1][2]) * (x[2] - p[1][2])
			- (C * (tr[1] + x[3] - t[1])) * (C * (tr[1] + x[3] - t[1]));
		fvec[2] = (x[0] - p[2][0]) * (x[0] - p[2][0]) + (x[1] - p[2][1]) * (x[1] - p[2][1]) + (x[2] - p[2][2]) * (x[2] - p[2][2])
			- (C * (tr[2] + x[3] - t[2])) * (C * (tr[2] + x[3] - t[2]));
		fvec[3] = (x[0] - p[3][0]) * (x[0] - p[3][0]) + (x[1] - p[3][1]) * (x[1] - p[3][1]) + (x[2] - p[3][2]) * (x[2] - p[3][2])
			- (C * (tr[3] + x[3] - t[3])) * (C * (tr[3] + x[3] - t[3]));

		/********************************/
	}
	// Jacobi matrix J(x)
	else if (*iflag == 2) {
		/********************************/

		fjac[0] = 2 * (x[0] - p[0][0]); 	fjac[4] = 2 * (x[1] - p[0][1]);		fjac[8] = 2 * (x[2] - p[0][2]);		fjac[12] = -2 * C * (C * (tr[0] + x[3] - t[0]));
		fjac[1] = 2 * (x[0] - p[1][0]);		fjac[5] = 2 * (x[1] - p[1][1]);		fjac[9] = 2 * (x[2] - p[1][2]);		fjac[13] = -2 * C * (C * (tr[1] + x[3] - t[1]));
		fjac[2] = 2 * (x[0] - p[2][0]);		fjac[6] = 2 * (x[1] - p[2][1]);		fjac[10] = 2 * (x[2] - p[2][2]);	fjac[14] = -2 * C * (C * (tr[2] + x[3] - t[2]));
		fjac[3] = 2 * (x[0] - p[3][0]);		fjac[7] = 2 * (x[1] - p[3][1]);		fjac[11] = 2 * (x[2] - p[3][2]);	fjac[15] = -2 * C * (C * (tr[3] + x[3] - t[3]));

		/********************************/
	}
}

void fcn3_1_hybrd1(int* n, double* x, double* fvec, int* iflag)
{
	// origin function F(x)
		
	fvec[0] = (x[0] - p[0][0]) * (x[0] - p[0][0]) + (x[1] - p[0][1]) * (x[1] - p[0][1]) + (x[2] - p[0][2]) * (x[2] - p[0][2])
		- (C * (tr[0] + x[3] - t[0])) * (C * (tr[0] + x[3] - t[0]));
	fvec[1] = (x[0] - p[1][0]) * (x[0] - p[1][0]) + (x[1] - p[1][1]) * (x[1] - p[1][1]) + (x[2] - p[1][2]) * (x[2] - p[1][2])
		- (C * (tr[1] + x[3] - t[1])) * (C * (tr[1] + x[3] - t[1]));
	fvec[2] = (x[0] - p[2][0]) * (x[0] - p[2][0]) + (x[1] - p[2][1]) * (x[1] - p[2][1]) + (x[2] - p[2][2]) * (x[2] - p[2][2])
		- (C * (tr[2] + x[3] - t[2])) * (C * (tr[2] + x[3] - t[2]));
	fvec[3] = (x[0] - p[3][0]) * (x[0] - p[3][0]) + (x[1] - p[3][1]) * (x[1] - p[3][1]) + (x[2] - p[3][2]) * (x[2] - p[3][2])
		- (C * (tr[3] + x[3] - t[3])) * (C * (tr[3] + x[3] - t[3]));

}

void practice3_1() {
	int n = SOLNUMS;
	double x[SOLNUMS] = { 0.0, 0.0, 0.0, 0.0 };	//need to initilize x0
	double fvec[SOLNUMS], fjac[MATCOLS * MATROWS];
	int ldfjac = SOLNUMS;
	double tol = TOLERANCE;
	int info;
	double j_wa[(SOLNUMS * (SOLNUMS + 13)) / 2];
	int j_lwa = (SOLNUMS * (SOLNUMS + 13)) / 2;
	double d_wa[(SOLNUMS * (3 * SOLNUMS + 13)) / 2];
	int d_lwa = (SOLNUMS * (3 * SOLNUMS + 13)) / 2;
	char readfile[256];
	char writefile1[256];
	char writefile2[256];

	printf("***** Problem 3-1 *****\n");

	for (int iter = 0; iter < 3; iter++) {
		sprintf(readfile, "GPS_signal_%d.txt", iter);
		sprintf(writefile1, "GPS_position_3-1_%d.txt", iter);
		sprintf(writefile2, "GPS_position_3-2_%d.txt", iter);

		FILE* fp_r = fopen(readfile, "r");
		if (fp_r == NULL) {
			printf("%s file open error...\n", readfile);
			return;
		}

		FILE* fp_w1 = fopen(writefile1, "w");
		if (fp_w1 == NULL) {
			printf("%s file open error...\n", writefile1);
			return;
		}

		FILE* fp_w2 = fopen(writefile2, "w");
		if (fp_w2 == NULL) {
			printf("%s file open error...\n", writefile2);
			return;
		}
		printf("GPS_SIGNAL_%d\n", iter);
		fscanf(fp_r, "%lf %lf", &C, &b);
		for (int i = 0; i < 4; i++) {
			fscanf(fp_r, "%lf %lf %lf", &p[i][0], &p[i][1], &p[i][2]);
			fscanf(fp_r, "%lf %lf", &t[i], &tr[i]);
		}

		printf("Enter initial value: ");
		scanf("%lf %lf %lf", &x[0], &x[1], &x[2]);
		x[3] = b;

		hybrj1_(fcn3_1_hybrj1, &n, x, fvec, fjac, &ldfjac, &tol, &info, j_wa, &j_lwa);

		printf("Hybrj1 method\n\n");
		print1(fp_w1, x, info);

		hybrd1_(fcn3_1_hybrd1, &n, x, fvec, &tol, &info, d_wa, &d_lwa);

		printf("Hybrd1 method\n\n");
		print1(fp_w2, x, info);
		
		if (fp_r != NULL) fclose(fp_r);
		if (fp_w1 != NULL) fclose(fp_w1);
		if (fp_w2 != NULL) fclose(fp_w2);
	}

}
