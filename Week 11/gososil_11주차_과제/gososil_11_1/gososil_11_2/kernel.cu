#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <Windows.h>
#include <time.h>
#include <assert.h>

#define CUDA_CALL(x) { const cudaError_t a = (x); if(a != cudaSuccess) { printf("\nCuda Error: %s (err_num=%d) at line:%d\n", cudaGetErrorString(a), a, __LINE__); cudaDeviceReset(); assert(0);}}
typedef float TIMER_T;

#define USE_CPU_TIMER 1
#define USE_GPU_TIMER 1

#if USE_CPU_TIMER == 1
__int64 start, freq, end;
#define CHECK_TIME_START { QueryPerformanceFrequency((LARGE_INTEGER*)&freq); QueryPerformanceCounter((LARGE_INTEGER*)&start); }
#define CHECK_TIME_END(a) { QueryPerformanceCounter((LARGE_INTEGER*)&end); a = (float)((float)(end - start) / (freq / 1000.0f)); }
#else
#define CHECK_TIME_START
#define CHECK_TIME_END(a)
#endif

#if USE_GPU_TIMER == 1
cudaEvent_t cuda_timer_start, cuda_timer_stop;
#define CUDA_STREAM_0 (0)

void create_device_timer()
{
    CUDA_CALL(cudaEventCreate(&cuda_timer_start));
    CUDA_CALL(cudaEventCreate(&cuda_timer_stop));
}

void destroy_device_timer()
{
    CUDA_CALL(cudaEventDestroy(cuda_timer_start));
    CUDA_CALL(cudaEventDestroy(cuda_timer_stop));
}

inline void start_device_timer()
{
    cudaEventRecord(cuda_timer_start, CUDA_STREAM_0);
}

inline TIMER_T stop_device_timer()
{
    TIMER_T ms;
    cudaEventRecord(cuda_timer_stop, CUDA_STREAM_0);
    cudaEventSynchronize(cuda_timer_stop);

    cudaEventElapsedTime(&ms, cuda_timer_start, cuda_timer_stop);
    return ms;
}

#define CHECK_TIME_INIT_GPU() { create_device_timer(); }
#define CHECK_TIME_START_GPU() { start_device_timer(); }
#define CHECK_TIME_END_GPU(a) { a = stop_device_timer(); }
#define CHECK_TIME_DEST_GPU() { destroy_device_timer(); }
#else
#define CHECK_TIME_INIT_GPU()
#define CHECK_TIME_START_GPU()
#define CHECK_TIME_END_GPU(a)
#define CHECK_TIME_DEST_GPU()
#endif

#define BLOCK_SIZE 128
#define N_EQUATIONS (1 << 26)
#define EPSILON 1

TIMER_T compute_time = 0;
TIMER_T device_time = 0;

int N;
float* A;
float* B;
float* C;
float* X0_CPU; float* X1_CPU;
float* FX0_CPU; float* FX1_CPU;
float* X0_GPU; float* X1_GPU;
float* FX0_GPU; float* FX1_GPU;

cudaError_t find_roots_GPU(float* A, float* B, float* C,
    float* X0, float* X1, float* FX0, float* FX1, int n);

__global__ void Find_Root_Kernel(float* A, float* B, float* C,
    float* X0, float* X1, float* FX0, float* FX1, int n)
{
    /*Thread Idx*/
    //int row = blockDim.y * blockIdx.y + threadIdx.y;
    //int col = blockDim.x * blockIdx.x + threadIdx.x;
    //int tid = gridDim.x * blockDim.x * row + col;
    int tid = blockIdx.x* blockDim.x + threadIdx.x;

    float a, b, c, d, x0, x1, tmp;
    a = A[tid]; b = B[tid]; c = C[tid];
    d = sqrtf(b * b - 4.0f * a * c);
    tmp = 1.0f / (2.0f * a);
    X0[tid] = x0 = (-b - d) * tmp;
    X1[tid] = x1 = (-b + d) * tmp;
    FX0[tid] = (a * x0 + b) * x0 + c;
    FX1[tid] = (a * x1 + b) * x1 + c;
}


void find_roots_CPU(float* A, float* B, float* C,
    float* X0, float* X1, float* FX0, float* FX1, int n) {

    int i;
    float a, b, c, d, x0, x1, tmp;
    for (i = 0; i < n; i++) {
        a = A[i]; b = B[i]; c = C[i];
        d = sqrtf(b * b - 4.0f * a * c);
        tmp = 1.0f / (2.0f * a);
        X0[i] = x0 = (-b - d) * tmp;
        X1[i] = x1 = (-b + d) * tmp;
        FX0[i] = (a * x0 + b) * x0 + c;
        FX1[i] = (a * x1 + b) * x1 + c;
    }  
}

void init_bin_file(int size) {
    srand((unsigned)time(NULL));
    printf("***Binary File init Start!!\n");
    FILE *fpa = fopen("a.bin", "wb");
    FILE* fpb = fopen("b.bin", "wb");
    FILE* fpc = fopen("c.bin", "wb");

    fwrite(&size, sizeof(int), 1, fpa);
    fwrite(&size, sizeof(int), 1, fpb);
    fwrite(&size, sizeof(int), 1, fpc);

    float random_a; float random_b; float random_c;
    for (int i = 0; i < size; i++) {
        random_a = (float)(((double)rand() / RAND_MAX) * 40.0f - 20.0f);
        random_b = (float)(((double)rand() / RAND_MAX) * 40.0f - 20.0f);
        random_c = (float)(((double)rand() / RAND_MAX) *40.0f - 20.0f);

        if (random_b * random_b < random_a * random_c * 4.0f) {
            i--;
            continue;
        }

        fwrite(&random_a, sizeof(float), 1, fpa);
        fwrite(&random_b, sizeof(float), 1, fpb);
        fwrite(&random_c, sizeof(float), 1, fpc);
    }

    fclose(fpa);
    fclose(fpb);
    fclose(fpc);
    printf("***Binary File init End!!\n\n");
}

void read_bin_file()
{
    printf("***Binary File Read Start!!\n");
    FILE* fpa = fopen("a.bin", "rb");
    FILE* fpb = fopen("b.bin", "rb");
    FILE* fpc = fopen("c.bin", "rb");
    /*Todo*/
    fread(&N, sizeof(int), 1, fpa);
    fread(&N, sizeof(int), 1, fpb);
    fread(&N, sizeof(int), 1, fpc);

    A = new float[N];
    B = new float[N];
    C = new float[N];
    X0_CPU = new float[N]; X1_CPU = new float[N];
    FX0_CPU = new float[N]; FX1_CPU = new float[N];
    X0_GPU = new float[N]; X1_GPU = new float[N];
    FX0_GPU = new float[N]; FX1_GPU = new float[N];

    fread(A, sizeof(float), N, fpa);
    fread(B, sizeof(float), N, fpb);
    fread(C, sizeof(float), N, fpc);

    fclose(fpa);
    fclose(fpb);
    fclose(fpc);
    printf("***Binary File Read End!!\n\n");
}

int main()
{
    init_bin_file(N_EQUATIONS);
    read_bin_file();

    printf("The problem size is %d.\n", N);
    int i;

    //CPU 
    printf("***find_roots_CPU Start!!\n");
    CHECK_TIME_START;
    find_roots_CPU(A, B, C, X0_CPU, X1_CPU, FX0_CPU, FX1_CPU, N);
    CHECK_TIME_END(compute_time);
    printf("***find_roots_CPU End!!\n\n");

    printf("CPU time = %.6f\n", compute_time);
    //GPU

    printf("***find_roots_GPU Start!!\n");
    find_roots_GPU(A, B, C, X0_GPU, X1_GPU, FX0_GPU, FX1_GPU, N);
    printf("***find_roots_GPU End!!\n\n");

    printf("GPU time = %.6f\n", device_time);

    for (i = 0; i < N; i++)
    {
        if (fabs(FX0_CPU[i]) > EPSILON || fabs(FX1_CPU[i]) > EPSILON || fabs(FX0_GPU[i]) > EPSILON || fabs(FX1_GPU[i]) > EPSILON)
        {
            printf("%d\n", i);
            printf("%f != %f || %f != %f\n", FX0_CPU[i], FX1_CPU[i], FX0_GPU[i], FX1_GPU[i]);
            break;
        }
    }

    if (i == N)
        printf("***Kernel execution Success!!\n\n");

    // Write the output array into the output file.
    FILE* fpx0 = fopen("X0.bin", "wb");
    FILE* fpx1 = fopen("X1.bin", "wb");
    FILE* fpfx0 = fopen("FX0.bin", "wb");
    FILE* fpfx1 = fopen("FX1.bin", "wb");

    if (!fpx0 || !fpx1 || !fpfx0 || !fpfx1) {
        fprintf(stderr, "Error: cannot open the output file...\n");
        exit(-1);
    }

    fwrite(&N, sizeof(int), 1, fpx0);
    fwrite(&N, sizeof(int), 1, fpx1);
    fwrite(&N, sizeof(int), 1, fpfx0);
    fwrite(&N, sizeof(int), 1, fpfx1);

    for (i = 0; i < N; i++) {
        fwrite(&X0_GPU[i], sizeof(float), 1, fpx0);
        fwrite(&X1_GPU[i], sizeof(float), 1, fpx1);
        fwrite(&FX0_GPU[i], sizeof(float), 1, fpfx0);
        fwrite(&FX1_GPU[i], sizeof(float), 1, fpfx1);
    }

    fclose(fpx0);
    fclose(fpx1);
    fclose(fpfx0);
    fclose(fpfx1);
    printf("end!!\n\n");

    return 0;
}


cudaError_t find_roots_GPU(float* A, float* B, float* C,
    float* X0, float* X1, float* FX0, float* FX1, int n)
{
    float* copy_a;
    float* copy_b;
    float* copy_c;
    float* copy_x0;
    float* copy_x1;
    float* copy_fx0;
    float* copy_fx1;

    cudaError_t cudaStatus;
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    size_t size = n * sizeof(float);
    CUDA_CALL(cudaMalloc(&copy_a, size));
    CUDA_CALL(cudaMemcpy(copy_a, A, size, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&copy_b, size));
    CUDA_CALL(cudaMemcpy(copy_b, B, size, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMalloc(&copy_c, size));
    CUDA_CALL(cudaMemcpy(copy_c, C, size, cudaMemcpyHostToDevice));

    CUDA_CALL(cudaMalloc(&copy_x0, size));
    CUDA_CALL(cudaMalloc(&copy_x1, size));
    CUDA_CALL(cudaMalloc(&copy_fx0, size));
    CUDA_CALL(cudaMalloc(&copy_fx1, size));

    dim3 dimBlock(BLOCK_SIZE);
    dim3 dimGrid(N / dimBlock.x);

    CHECK_TIME_INIT_GPU();
    CHECK_TIME_START_GPU();

    Find_Root_Kernel << < dimGrid, dimBlock >> > (copy_a, copy_b, copy_c, copy_x0, copy_x1, copy_fx0, copy_fx1, N);

    CHECK_TIME_END_GPU(device_time)
    CHECK_TIME_DEST_GPU();

    CUDA_CALL(cudaGetLastError());
    CUDA_CALL(cudaDeviceSynchronize());

    CUDA_CALL(cudaMemcpy(X0, copy_x0, size, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(X1, copy_x1, size, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(FX0, copy_fx0, size, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(FX1, copy_fx1, size, cudaMemcpyDeviceToHost));

        Error:
    cudaFree(copy_a);
    cudaFree(copy_b);
    cudaFree(copy_c);
    cudaFree(copy_x0);
    cudaFree(copy_x1);
    cudaFree(copy_fx0);
    cudaFree(copy_fx1);

    return cudaStatus;
}