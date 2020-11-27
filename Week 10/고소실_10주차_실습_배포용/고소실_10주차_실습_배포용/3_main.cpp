#include <stdio.h>
#include <random>
#include <time.h>
#include <math.h>
#include <time.h>
#include <Windows.h>

#define N 50
float Taylor_series(double x, int n);
double Taylor_series_ex(double x, int n);
__int64 start, freq, end;
#define CHECK_TIME_START QueryPerformanceFrequency((LARGE_INTEGER*)&freq); QueryPerformanceCounter((LARGE_INTEGER*)&start)
#define CHECK_TIME_END(a) QueryPerformanceCounter((LARGE_INTEGER*)&end); a = (float)((float)(end - start) / (freq / 1000.0f))
float compute_time;

void main(void)
{
	CHECK_TIME_START;
	printf("*** f<-8.3> = %.6e \n", Taylor_series(-8.3, N));
	CHECK_TIME_END(compute_time);
	printf("***Honer, float : %f ms\n\n", compute_time);

	CHECK_TIME_START;
	printf("*** f<-8.3> = %.6e \n", Taylor_series_ex(-8.3, N));
	CHECK_TIME_END(compute_time);
	printf("***Honer, double : %f ms\n\n", compute_time);

	printf("*** f<-8.3> = %.6e \n", pow(2.71828182846, -8.3));
}




float Taylor_series(double x, int n)
{

	/* Float */
	float ans = x / (n - 1);

	for (int i = n - 2; i > 0; i--) {
		ans = (ans + 1) * (x / i);
	}
	ans += 1;

	return ans;
}


double Taylor_series_ex(double x, int n)
{
	/* Double */
	double ans = x / (n - 1);

	for (int i = n - 2; i > 0; i--) {
		ans = (ans + 1) * (x / i);
	}
	ans += 1;

	return ans;
}
