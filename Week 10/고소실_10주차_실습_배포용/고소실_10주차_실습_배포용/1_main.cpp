#include <stdio.h>
#include <string.h>
#include <random>
#include <time.h>

#include <math.h>
#include <time.h>
#include <Windows.h>

__int64 start, freq, end;
#define CHECK_TIME_START QueryPerformanceFrequency((LARGE_INTEGER*)&freq); QueryPerformanceCounter((LARGE_INTEGER*)&start)
#define CHECK_TIME_END(a) QueryPerformanceCounter((LARGE_INTEGER*)&end); a = (float)((float)(end - start) / (freq / 1000.0f))
float compute_time;


#define MATDIM 1024
double *pMatA, *pMatB, *pMatC;
void MultiplySquareMatrices_1(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize);
void MultiplySquareMatrices_2(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize);
void MultiplySquareMatrices_3(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize);
void MultiplySquareMatrices_4(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize);

void init_MatMat(void);

void main(void)
{
	init_MatMat();

	CHECK_TIME_START;
	MultiplySquareMatrices_1(pMatC, pMatA, pMatB, MATDIM);
	//Sleep(500);
	CHECK_TIME_END(compute_time);
	printf("MultiplySquareMatrices_1 : %f ms\n", compute_time);

	CHECK_TIME_START;
	MultiplySquareMatrices_2(pMatC, pMatA, pMatB, MATDIM);
	//Sleep(500);
	CHECK_TIME_END(compute_time);
	printf("MultiplySquareMatrices_2 = %f ms\n", compute_time);

	CHECK_TIME_START;
	MultiplySquareMatrices_3(pMatC, pMatA, pMatB, MATDIM);
	//Sleep(500);
	CHECK_TIME_END(compute_time);
	printf("MultiplySquareMatrices_3 = %f ms\n", compute_time);

	CHECK_TIME_START;
	MultiplySquareMatrices_4(pMatC, pMatA, pMatB, MATDIM);
	//Sleep(500);
	CHECK_TIME_END(compute_time);
	printf("MultiplySquareMatrices_4 = %f ms\n", compute_time);
}

void MultiplySquareMatrices_1(double *pDestMatrix, double *pLeftMatrix ,double *pRightMatrix, int MatSize)
{
	for (int i = 0; i < MatSize; i++) {
		for (int j = 0; j < MatSize; j++) {
			pDestMatrix[i * MatSize + j] = 0;
			for (int k = 0; k < MatSize; k++) {
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k] * pRightMatrix[k * MatSize + j];
			}
		}
	}
}

void MultiplySquareMatrices_2(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize)
{
	double* pTransMatrix = (double*)malloc(sizeof(double) * MatSize * MatSize);
	
	for (int i = 0; i < MatSize; i++) {
		for (int j = 0; j < MatSize; j++) {
			pTransMatrix[i * MatSize + j] = pRightMatrix[j * MatSize + i];
		}
	}

	for (int i = 0; i < MatSize; i++) {
		for (int j = 0; j < MatSize; j++) {
			pDestMatrix[i * MatSize + j] = 0;
			for (int k = 0; k < MatSize; k++) {
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k] * pTransMatrix[i * MatSize + k];
			}
		}
	}

	free(pTransMatrix);
}

void MultiplySquareMatrices_3(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize)
{
	double* pTransMatrix = (double*)malloc(sizeof(double) * MatSize * MatSize);

	for (int i = 0; i < MatSize; i++) {
		for (int j = 0; j < MatSize; j++) {
			pTransMatrix[i * MatSize + j] = pRightMatrix[j * MatSize + i];
		}
	}

	for (int i = 0; i < MatSize; i++) {
		for (int j = 0; j < MatSize; j++) {
			pDestMatrix[i * MatSize + j] = 0;
			for (int k = 0; k < MatSize; k+=8) {
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k] * pTransMatrix[i * MatSize + k];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 1] * pTransMatrix[i * MatSize + k + 1];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 2] * pTransMatrix[i * MatSize + k + 2];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 3] * pTransMatrix[i * MatSize + k + 3];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 4] * pTransMatrix[i * MatSize + k + 4];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 5] * pTransMatrix[i * MatSize + k + 5];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 6] * pTransMatrix[i * MatSize + k + 6];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 7] * pTransMatrix[i * MatSize + k + 7];
				
			}
		}
	}

	free(pTransMatrix);
}

void MultiplySquareMatrices_4(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize)
{
	double* pTransMatrix = (double*)malloc(sizeof(double) * MatSize * MatSize);

	for (int i = 0; i < MatSize; i++) {
		for (int j = 0; j < MatSize; j++) {
			pTransMatrix[i * MatSize + j] = pRightMatrix[j * MatSize + i];
		}
	}

	for (int i = 0; i < MatSize; i++) {
		for (int j = 0; j < MatSize; j++) {
			pDestMatrix[i * MatSize + j] = 0;
			for (int k = 0; k < MatSize; k += 16) {
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k] * pTransMatrix[i * MatSize + k];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 1] * pTransMatrix[i * MatSize + k + 1];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 2] * pTransMatrix[i * MatSize + k + 2];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 3] * pTransMatrix[i * MatSize + k + 3];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 4] * pTransMatrix[i * MatSize + k + 4];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 5] * pTransMatrix[i * MatSize + k + 5];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 6] * pTransMatrix[i * MatSize + k + 6];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 7] * pTransMatrix[i * MatSize + k + 7];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 8] * pTransMatrix[i * MatSize + k + 8];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 9] * pTransMatrix[i * MatSize + k + 9];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 10] * pTransMatrix[i * MatSize + k + 10];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 11] * pTransMatrix[i * MatSize + k + 11];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 8] * pTransMatrix[i * MatSize + k + 12];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 9] * pTransMatrix[i * MatSize + k + 13];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 10] * pTransMatrix[i * MatSize + k + 14];
				pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k + 11] * pTransMatrix[i * MatSize + k + 15];
			}
		}
	}

	free(pTransMatrix);
}

void init_MatMat(void)
{
	double *ptr;
	pMatA = (double *)malloc(sizeof(double)*MATDIM*MATDIM);
	pMatB = (double *)malloc(sizeof(double)*MATDIM*MATDIM);
	pMatC = (double *)malloc(sizeof(double)*MATDIM*MATDIM);
	srand((unsigned)time(NULL));
	ptr = pMatA;
	for (int i = 0; i < MATDIM*MATDIM; i++)
		*ptr++ = (double)rand() / ((double)RAND_MAX);
	ptr = pMatB;
	for (int i = 0; i < MATDIM*MATDIM; i++)
		*ptr++ = (double)rand() / ((double)RAND_MAX);
}
