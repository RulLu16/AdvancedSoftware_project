#pragma once

#include "cuda_code.cuh"

#include <stdio.h>
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <Windows.h>
#include <time.h>
#include <assert.h>

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

TIMER_T compute_time = 0;
TIMER_T device_time = 0;

#define Window 2
#define BLOCK_SIZE 4

__constant__ float constant_gaussian_kernel[ 25 ];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	
//	Gaussian ���͸��� �ϴ� Ŀ��
//	shared memory�� ������� �ʴ´�
//	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void Gaussian_kernel_no_shared(IN unsigned char *d_bitmaps, OUT unsigned char *d_Gaussian, long width, long height) {
	/*Todo*/
	int row = blockDim.y * blockIdx.y + threadIdx.y;
	int col = blockDim.x * blockIdx.x + threadIdx.x;
	int w = floor((double)5 / 2);
	double mean = 0;

	for (int l = -w; l <= w; l++) {
		for (int k = -w; k <= w; k++) {

			if (row + l < 0 || row + l >= height || col + k < 0 || col + k >= width) continue;
			else
				mean += d_bitmaps[(row + l) * width + col +k] * constant_gaussian_kernel[(l + w) * 5 + k + w];		
		}
	}
	d_Gaussian[row*width + col] = mean;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	
//	Gaussian ���͸��� �ϴ� Ŀ��
//	shared memory�� ����Ѵ�.
//	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern __shared__ unsigned char sharedBuffer[];
__global__ void Gaussian_kernel_shared(INOUT unsigned char *d_bitmaps, OUT unsigned char *d_Gaussian, long width, long height) {
	/*Todo*/
	int row = blockDim.y * blockIdx.y + threadIdx.y;
	int col = blockDim.x * blockIdx.x + threadIdx.x;

	const unsigned block_id = blockIdx.y * gridDim.x + blockIdx.x;
	const unsigned thread_id = threadIdx.y * blockDim.x + threadIdx.x;
	const unsigned id = block_id * BLOCK_SIZE + thread_id;

	int w = 2;// floor((double)5 / 2);
	double mean = 0;
	
	if (thread_id == 0) { // ����
		for (int i = 1; i <= w; i++) {
			if (id - i < 0) {
				for (int j = 0; j < width; j++) {
					sharedBuffer[(w - i) * width + j] = 0;
				}
			}
			else if (id - i >= 0 && id - i <= height) {
				for (int j = 0; j < width; j++)
					sharedBuffer[(w - i) * width + j] = d_bitmaps[(id - i) * width + j];
			}
		}

	}
	/*else if (thread_id % BLOCK_SIZE == 0) { // ���� ��

	}
	else if (thread_id % BLOCK_SIZE == BLOCK_SIZE - 1) { // ������ ��

	}*/
	if (thread_id == BLOCK_SIZE - 1) { // �Ʒ� ��
		for (int i = 1; i <= w; i++) {

			if (id +i >= height) {
				for (int j = 0; j < width; j++)
					sharedBuffer[(BLOCK_SIZE + w + i - 1) * width + j] = 0;
			}
			else if (id + i < height) {
				for (int j = 0; j < width; j++)
					sharedBuffer[(BLOCK_SIZE + w + i - 1) * width + j] = d_bitmaps[(id + i) * width + j];
			}
		}
	}

	for (int i = 0; i < width; i++)
		//sharedBuffer[(thread_id + w) * width + i] = 0; d_bitmaps[id * width + i];
	//sharedBuffer[(thread_id + w) * width] = d_bitmaps[id * width];

	__syncthreads();

	for (int i = 0; i < width; i++) {
		mean = 0;

		for (int l = -w; l <= w; l++) {
			for (int k = -w; k <= w; k++) {

				if (i + l >= width || i + l < 0) continue;
				//else
					//mean += sharedBuffer[(thread_id + w + k) * width + i + l] * constant_gaussian_kernel[(k + w) * 5 + l + w]; // error
				//printf("%lf\n", constant_gaussian_kernel[(l + w) * 5 + k + w]);
			}
		}
		//d_Gaussian[id * width + i] = mean;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	
//	Constant variable �� gaussian kernel�� �����ϴ� �Լ�
//	�Ŀ� gaussian filtering ���� ����Ѵ�.
//	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Set_Gaussain_Kernel(){
	float _1 = 1.0f / 256.0f;
	float _4 = _1 * 4;
	float _6 = _1 * 6;
	float _16 = _1 * 16;
	float _24 = _1 * 24;
	float _36 = _1 * 36;

	float *p_gaussian_kernel = new float[25];

	p_gaussian_kernel[0] = p_gaussian_kernel[4] = p_gaussian_kernel[20] = p_gaussian_kernel[24] = _1;
	p_gaussian_kernel[1] = p_gaussian_kernel[3] = p_gaussian_kernel[5] = p_gaussian_kernel[9]= _4;
	p_gaussian_kernel[15] = p_gaussian_kernel[19] = p_gaussian_kernel[21] = p_gaussian_kernel[23] = _4;
	p_gaussian_kernel[2] = p_gaussian_kernel[10] = p_gaussian_kernel[14] = p_gaussian_kernel[22] = _6;
	p_gaussian_kernel[6] = p_gaussian_kernel[8] = p_gaussian_kernel[16] = p_gaussian_kernel[18] = _16;
	p_gaussian_kernel[7] = p_gaussian_kernel[11] =p_gaussian_kernel[13] = p_gaussian_kernel[17] = _24;
	p_gaussian_kernel[12] = _36;

	cudaMemcpyToSymbol( constant_gaussian_kernel, p_gaussian_kernel, sizeof( float ) * 25 );

	delete[] p_gaussian_kernel;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	
//	Ŀ���� �����ϱ� �� �ʿ��� �ڷ�� �غ� �� Ŀ���� ������ ����̽��� ����
//	Shared_flag �Է� �� NO_SHARED �� SHARED �� �� ���� ��ũ�θ� ������
//	flag���� �´� Ŀ���� ����
//	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float Do_Gaussian_on_GPU(IN unsigned char *p_bitmaps, OUT unsigned char *p_Gaussian, long width, long height, int Shared_flag)
{
	/*Todo*/
	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	CUDA_CALL(cudaStatus);

	Set_Gaussain_Kernel();
	unsigned char* d_bitmaps, * d_Gaussian;
	size_t mem_size;

	mem_size = width * height * sizeof(unsigned char);
	CUDA_CALL(cudaMalloc(&d_bitmaps, mem_size));
	CUDA_CALL(cudaMemcpy(d_bitmaps, p_bitmaps, mem_size, cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc(&d_Gaussian, mem_size));

	dim3 blockDIm(BLOCK_SIZE, BLOCK_SIZE);
	dim3 gridDim(width / blockDIm.x, height / blockDIm.y);
	CHECK_TIME_INIT_GPU();
	CHECK_TIME_START_GPU();

	switch (Shared_flag)
	{
	case NO_SHARED:
		Gaussian_kernel_no_shared << <gridDim, blockDIm >> > (d_bitmaps, d_Gaussian, width, height);
		break;
	case SHARED:
		Gaussian_kernel_shared << <gridDim, blockDIm, sizeof(unsigned char)* (BLOCK_SIZE * BLOCK_SIZE + (BLOCK_SIZE + 4) * width) >> > (d_bitmaps, d_Gaussian, width, height);
		break;
	}

	CUDA_CALL(cudaStatus = cudaDeviceSynchronize());
	CHECK_TIME_END_GPU(device_time);
	CHECK_TIME_DEST_GPU();

	CUDA_CALL(cudaMemcpy(p_Gaussian, d_Gaussian, mem_size, cudaMemcpyDeviceToHost));

	cudaFree(d_bitmaps);
	cudaFree(d_Gaussian);

	return device_time;
}