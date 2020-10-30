#include "my_solver.h"

void complexNum(double* a, double* b, double c, double d) {
	double real = *a * c - *b * d;
	double comp = *a * d + *b * c;

	*a = real;
	*b = comp;
}

void practice3_2() {
	char readfile[256];
	char writefile[256];

	int NCOEF, DEGREE;

	double* poly;
	double* zeror, * zeroi;
	long int fail;
	int i;

	printf("\n***** Problem 3-2 *****\n");
	for (int t = 1; t <= 6; t++) {
		printf("\npolynomial_3-2_%d solution\n", t);
		sprintf(readfile, "polynomial_3-2_%d.txt", t);
		sprintf(writefile, "roots_3-2_%d.txt", t);

		FILE* fp_r = fopen(readfile, "r");
		if (fp_r == NULL) {
			printf("%s file open error...\n", readfile);
			return;
		}

		FILE* fp_w = fopen(writefile, "w");
		if (fp_w == NULL) {
			printf("%s file open error...\n", writefile);
			return;
		}

		fscanf(fp_r, "%d", &DEGREE);
		NCOEF = DEGREE + 1;

		poly = (double*)malloc(sizeof(double) * NCOEF);
		zeror = (double*)malloc(sizeof(double) * DEGREE);
		zeroi = (double*)malloc(sizeof(double) * DEGREE);


		/********************************/
		for (int i = 0; i < NCOEF; i++) {
			fscanf(fp_r, "%lf", &poly[i]);
		}
		rpoly_(poly, &DEGREE, zeror, zeroi, &fail);

		if (fail) {
			printf("Wrong!!\n\n");
			fprintf(fp_w, "Wrong!!\n");
			continue;
		}

		for (int i = 0; i < DEGREE; i++) {
			printf("root%d = %10f + %10f * i\n",i, zeror[i], zeroi[i]);
			fprintf(fp_w, "root%d = %10f + %10f * i\n", i, zeror[i], zeroi[i]);
		}
		printf("\n");
		fprintf(fp_w, "\n");

		for (int i = 0; i < DEGREE; i++) {
			printf("|f(x)| for root%d = ", i);
			fprintf(fp_w, "|f(x)| for root%d = ", i);

			double realResult = 0;
			double compResult = 0;
			for (int j = 0; j < DEGREE; j++) {
				double realmid = zeror[i];
				double compmid = zeroi[i];

				for (int k = DEGREE - j; k > 1; k--) {
					complexNum(&realmid, &compmid, zeror[i], zeroi[i]);
				}
				realmid *= poly[j];
				compmid *= poly[j];
				//printf("%10f %10f \n\n", realmid, compmid);

				realResult += realmid;
				compResult += compmid;
			}
			//realResult += poly[DEGREE - 1] * zeror[i];
			realResult += poly[NCOEF - 1];
			//compResult += poly[DEGREE - 1] * zeroi[i];
			printf("%10f + %10f * i\n", fabs(realResult), fabs(compResult));
			fprintf(fp_w, "%10f + %10f * i\n", fabs(realResult), fabs(compResult));
		}

		/********************************/


		free(zeroi);
		free(zeror);
		free(poly);

		if (fp_r != NULL) fclose(fp_r);
		if (fp_w != NULL) fclose(fp_w);
	}
}