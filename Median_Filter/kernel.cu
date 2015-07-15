#ifndef __CUDACC__  
#define __CUDACC__
#endif

#define N 128
#define BLOCK_DIM 4
#define T 32

#define BLOCK_W 32
#define BLOCK_H 32
#define R 9
#define D (R*2+1)
#define D2 D*D
#define TILE_W (BLOCK_W - (2*R))
#define TILE_H (BLOCK_H - (2*R))

//#define TILE_W 12
//#define TILE_H 12

//#define BLOCK_W (TILE_W + (2*R))
//#define BLOCK_H (TILE_H + (2*R))


#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>

#include "Grid.h"
#include "ReadWrite.h"
#include "Filter.h"


//fix for size 1 larger than block

cudaError_t addWithCuda(int *outputImage, const int *inputImage, unsigned int inputSize, unsigned int outputSize);

__global__ void addKernel(int *inputImage, int *outputImage, int width, int height)//int *inputImage, int *outputImage, unsigned int width, unsigned int height)
{
	//int idc, val, min, max, inf, equal, sup, mxinf, minsup, estim;

	//int ib = threadIdx.y;
	//int jb = threadIdx.x;
	//int idx_h = __mul24(ib+R, blockDim.x+2) + jb + 1;
	//int offset = __mul24(blockDim.x, R);

	//int j = __mul24(blockIdx.x, blockDim.x) + jb;
	//int i = __mul24(blockIdx.y, blockDim.y) + ib;
	//if (i < N && j < N)
	//{
	//	/*extern*/ __shared__ int buff[40 * 40];
	//	buff[idx_h] = inputImage[i * N + j];
	//	if (ib < R)
	//	{
	//		buff[idx_h - offset] = 0;//padding
	//	}
	//	else if (ib >= (blockDim.y-R))
	//	{
	//		buff[idx_h + offset] = 0;//padding
	//	}
	//	__syncthreads();
	//	printf("dim[%d][%d] | %-10d : %-9d, %-9d, %-9d, %-9d, %-9d, %-9d, %-9d, %-9d, %-9d \n", i, j, inputImage[(i) * N + j],
	//		buff[idx_h-offset-1], buff[idx_h-offset], buff[idx_h-offset+1],
	//		buff[idx_h-1], buff[idx_h], buff[idx_h+1], 
	//		buff[idx_h+offset-1], buff[idx_h+offset], buff[idx_h+offset]);
	//	min = max = buff[ib * blockDim.x + jb];

	//	/*for (idc = 0; idc < 2 * R + 1; ++idc)
	//	{
	//	val = buff[__mul24(ib + idc, blockDim.x) + jb];
	//	if (val < min) min = val;
	//	if (val > max) max = val;
	//	}

	//	while (1)
	//	{
	//	estim = (min + max) / 2;
	//	inf = sup = equal = 0;
	//	mxinf = min;
	//	minsup = max;
	//	for (idc = 0; idc < 2 * R + 1; ++idc)
	//	{
	//	val = buff[__mul24(ib + idc, blockDim.x) + jb];
	//	if (val < estim)
	//	{
	//	inf++;
	//	if (val > mxinf) mxinf = val;
	//	}
	//	else if (val > estim)
	//	{
	//	sup++;
	//	if (val < minsup) minsup = val;
	//	}
	//	else equal++;
	//	}
	//	if ((inf <= (R + 1)) && (sup <= (R + 1))) break;
	//	else if (inf > sup) max = mxinf;
	//	else min = minsup;
	//	}
	//	if (inf >= R + 1) val = mxinf;
	//	else if (inf + equal >= R + 1) val = estim;
	//	else val = minsup;*/

	//	outputImage[__mul24(j, N) + i] = 0;//val;
	//}

//	__shared__ int smem[16 * 16];
//
//	int x = blockIdx.x*TILE_W + threadIdx.x - R;
//
//	int y = blockIdx.y*TILE_H + threadIdx.y - R;
//
//	//clamp to edge of image
////	if (blockIdx.x*blockDim.x + threadIdx.x < N + R && blockIdx.y*blockDim.y + threadIdx.y < N + R)
//	//{
//		x = max(0, x);
//
//		x = min(x, width - 1);
//
//		y = max(y, 0);
//
//		y = min(y, height - 1);
//
//
//
//		unsigned int index = y*width + x;
//
//		unsigned int bindex = threadIdx.y * blockDim.y + threadIdx.x;
//
//
//
//		//each thread copies its pixel of the block to shared memory
//
//		smem[bindex] = inputImage[index];
//
//		__syncthreads();
//
//		if ((threadIdx.x >= R) && (threadIdx.x < (BLOCK_W - R)) &&
//			(threadIdx.y >= R) && (threadIdx.y < (BLOCK_H - R)))
//		{
//			
//			min = max = smem[ib * blockDim.x + jb];
//
//				for (idc = 0; idc < 2 * R + 1; ++idc)
//				{
//					val = smem[__mul24(ib + idc, blockDim.x) + jb];
//				if (val < min) min = val;
//				if (val > max) max = val;
//				}
//
//				while (1)
//				{
//				estim = (min + max) / 2;
//				inf = sup = equal = 0;
//				mxinf = min;
//				minsup = max;
//				for (idc = 0; idc < 2 * R + 1; ++idc)
//				{
//					val = smem[__mul24(ib + idc, blockDim.x) + jb];
//				if (val < estim)
//				{
//				inf++;
//				if (val > mxinf) mxinf = val;
//				}
//				else if (val > estim)
//				{
//				sup++;
//				if (val < minsup) minsup = val;
//				}
//				else equal++;
//				}
//				if ((inf <= (R + 1)) && (sup <= (R + 1))) break;
//				else if (inf > sup) max = mxinf;
//				else min = minsup;
//				}
//				if (inf >= R + 1) val = mxinf;
//				else if (inf + equal >= R + 1) val = estim;
//				else val = minsup;
//
//			/*int      val, i, less, greater, equal,  min, max, guess, maxltguess, mingtguess;
//
//			min = max = smem[bindex];
//			for (i =-1; i<=1; i++) {
//			val = smem[__mul24(threadIdx.y + i, blockDim.x) + threadIdx.x];
//			if (val<min) min = val;
//			if (val> max) max = val;
//			}
//
//			while (1) {
//				guess = (min + max) / 2;
//				less = 0; greater = 0; equal = 0;
//				maxltguess = min;
//				mingtguess = max;
//				for (i = -1; i<=1; i++) {
//					val = smem[__mul24(threadIdx.y + i, blockDim.x) + threadIdx.x];
//					if (val<guess) {
//						less++;
//						if (val>maxltguess) maxltguess = val;
//					}
//					else if (val>guess) {
//						greater++;
//					if (val<mingtguess) mingtguess = val;
//					}
//					else equal++;
//				}
//				if (less <= (R + 1) && greater <= (R + 1)) break;
//				else if (less>greater) max = maxltguess;
//				else min = mingtguess;
//			}
//			if (less >= (R + 1)) val = maxltguess;
//			else if (less + equal >= (R + 1)) val = guess;
//			else val = mingtguess;*/
//			/*printf("dim[%d][%d] | %-10d : %-10d, %-10d, %-10d, %-10d, %-10d, %-10d, %-10d, %-10d, %-10d \n", y, x, inputImage[index],
//				smem[bindex - blockDim.x - 1], smem[bindex - blockDim.x], smem[bindex - blockDim.x + 1],
//				smem[bindex - 1], smem[bindex], smem[bindex + 1],
//				smem[bindex + blockDim.x - 1], smem[bindex + blockDim.x], smem[bindex + blockDim.x + 1]);*/
//			//outputImage[index] =  smem[bindex];
//			outputImage[index] = val;
//
//		}
	//}

	 //extern __shared__ int data[];
	__shared__ int data[(R * 2 + 1) * (T + 2 * R)];
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int d_col = col + R;
	int d_row = blockIdx.y;

	if (row < N && col < N)
	{

		for (int r = -R; r <= R; ++r)
		{
			//col += r;
			int i_row = row + r;
			//col = max(0, col);

			//col = min(col, N - 1);

			i_row = max(i_row, 0);

			i_row = min(i_row, N - 1);
			//data[1028 + col] = 0; //inputImage[row*N + col];
			data[(R + r) * (T + 2 * R) + d_col] = inputImage[i_row*N + col];
			if (col == 0)
			{
				for (int i = 0; i <= R;++i)
					data[(R + r) * (T + 2 * R) + col+i] = inputImage[i_row*N + col];
			}
			if (col == N - 1 || col == blockDim.x)
			{
				for (int i = 0; i <= R; ++i)
					data[(R + r) * (T + 2 * R) + d_col + i] = inputImage[i_row*N + col];
			}
		}
		__syncthreads();
		/*printf("dim : %d | [%d][%d] | %-10d : %-9d, %-9d, %-9d, %-9d, %-9d, %-9d, %-9d, %-9d, %-9d \n", blockDim.y, d_row, col, inputImage[(row) * N + col], 
			data[(0) * 514 + d_col - 1], data[(0) * 514 + d_col], data[(0) * 514 + d_col + 1],
			data[(1) * 514 + d_col - 1], data[(1) * 514 + d_col], data[(1) * 514 + d_col + 1],
			data[(2) * 514 + d_col - 1], data[(2) * 514 + d_col], data[(2) * 514 + d_col + 1]);*/
		//printf("dim : %d | [%d][%d] | %-10d : %-9d, %-9d, %-9d\n", blockDim.y, d_row, N, inputImage[(row) * N + col], data[(0) * 514 + N], data[(1)* 514 + N], data[(2) * 514 + N]);
		//int k = 0;
		//int *window = new int[D*D];
		//for (int dy = -R; dy <= R; ++dy)
		//{
		//	for (int dx = -R; dx <= R; ++dx)
		//	{
		//		window[k++] = data[(R + dy) * (T + 2 * R) + d_col + dx];
		//	}
		//}
		//for (int j = 0; j < (D*D+1)/2; ++j)
		//{
		//	//   find position of minimum element
		//	int min = j;
		//	for (int l = j + 1; l < (D*D); ++l)
		//		if (window[l] < window[min])
		//			min = l;
		//	const int temp = window[j];
		//	window[j] = window[min];
		//	window[min] = temp;
		//}
		/*int      val, i, less, greater, equal,  min, max, guess, maxltguess, mingtguess;
		
		min = max = data[__mul24(threadIdx.y + R, blockDim.x + 2 * R) + threadIdx.x + R];
		for (i =0; i<D; i++) {
			val = data[__mul24(threadIdx.y + i, blockDim.x) + threadIdx.x];
			if (val<min) min = val;
			if (val> max) max = val;
		}
		while (1) {
			guess = (min + max) / 2;
			less = 0; greater = 0; equal = 0;
			maxltguess = min;
			mingtguess = max;
			for (i = 0; i<D; i++) {
				val = data[__mul24(threadIdx.y + i, blockDim.x) + threadIdx.x];
				if (val<guess) {
					less++;
					if (val>maxltguess) maxltguess = val;
				}
				else if (val>guess) {
					greater++;
				if (val<mingtguess) mingtguess = val;
				}
				else equal++;
			}
			if (less <= (R + 1) && greater <= (R + 1)) break;
			else if (less>greater) max = maxltguess;
			else min = mingtguess;
		}
		if (less >= (R + 1)) val = maxltguess;
		else if (less + equal >= (R + 1)) val = guess;
		else val = mingtguess;*/

		//outputImage[(row)* N + col] = val;
		//outputImage[(row)* N + col] = data[(R)* (T + 2 * R) + d_col/*(1) * 514 + d_col*/];
		outputImage[(row)* N + col] = 1;// window[(D*D) / 2];
	}
}

__global__ void median_Kernel_3x3(int *inputImage, int *outputImage, int width, int height)//int *inputImage, int *outputImage, unsigned int width, unsigned int height)
{
__shared__ int smem[(BLOCK_DIM + 2 * R)][(BLOCK_DIM + 2 * R)];
int x_s = threadIdx.x + R;
int y_s = threadIdx.y + R;

int x_a = blockIdx.x*blockDim.x + threadIdx.x + R;
int y_a = blockIdx.y*blockDim.y + threadIdx.y + R;

int x = blockIdx.x*blockDim.x + threadIdx.x;
int y = blockIdx.y*blockDim.y + threadIdx.y;
if (x_a < N + 2 * R && y_a < N + 2 * R)
{
	smem[y_s][x_s] = inputImage[y_a*(N + 2 * R) + x_a];
	if (threadIdx.x == 0 || threadIdx.x == BLOCK_DIM - 1)
	{
		for (int i = 0; i < R; ++i)
		{
			smem[threadIdx.y][threadIdx.x + i] = inputImage[y*(N + 2 * R) + x + i];
		}
	}
	if (threadIdx.y == 0 || threadIdx.y == BLOCK_DIM - 1)
	{
		for (int i = 0; i < R; ++i)
		{
			smem[threadIdx.y + i][threadIdx.x] = inputImage[(y + i)*(N + 2 * R) + x];
		}
	}

	__syncthreads();
	if (x < N && y < N)
	{
		/*printf("dim[%d][%d] | %-10d : \n%-10d, %-10d, %-10d, \n%-10d, %-10d, %-10d, \n%-10d, %-10d, %-10d \n", y, x, inputImage[y_a*(N + 2 * R) + x_a],
			smem[y_a - 1][x_a - 1], smem[y_a - 1][x_a], smem[y_a - 1][x_a + 1],
			smem[y_a][x_a - 1], smem[y_a][x_a], smem[y_a][x_a + 1],
			smem[y_a + 1][x_a - 1], smem[y_a+1][x_a], smem[y_a + 1][x_a + 1]);*/
		outputImage[y*width + x] = smem[y_s][x_s];
	}
}
//how did you change your seat
}

__global__ void median_Kernel_5x5(int *inputImage, int *outputImage, int width, int height)//int *inputImage, int *outputImage, unsigned int width, unsigned int height)
{
	__shared__ int smem[(BLOCK_DIM + 2 * R)][(BLOCK_DIM + 2 * R)];
	int x_s = threadIdx.x + R;
	int y_s = threadIdx.y + R;

	int x_a = blockIdx.x*blockDim.x + threadIdx.x + R;
	int y_a = blockIdx.y*blockDim.y + threadIdx.y + R;

	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y*blockDim.y + threadIdx.y;
	if (x < N && y < N)
	{
		smem[y_s][x_s] = inputImage[y*(N)+x];
		if (threadIdx.x == 0 && threadIdx.y == 0 && x != N - 1 && y != N - 1)
		{
			for (int j = 0; j < R; ++j)
				for (int i = 0; i < R; ++i)
					smem[threadIdx.y + j][threadIdx.x + i] = inputImage[y*(N)+x];
		}
		else if (threadIdx.x == BLOCK_DIM - 1 && threadIdx.y == BLOCK_DIM - 1)
		{
			for (int j = 0; j < R; ++j)
				for (int i = 0; i < R; ++i)
					smem[y_s + j + 1][x_s + i + 1] = inputImage[y*(N)+x];

		}
		else if (threadIdx.x == 0 && threadIdx.y == BLOCK_DIM - 1)
		{
			for (int j = 0; j < R; ++j)
				for (int i = 0; i < R; ++i)
					smem[y_s + j + 1][threadIdx.x + i] = inputImage[y*(N)+x];

		}
		else if (threadIdx.x == BLOCK_DIM - 1 && threadIdx.y == 0)
		{
			for (int j = 0; j < R; ++j)
				for (int i = 0; i < R; ++i)
					smem[threadIdx.y + j][x_s + i + 1] = inputImage[y*(N)+x];

		}

		if (threadIdx.x == 0 || x != N - 1)
		{
			for (int i = 0; i < R; ++i)
			{
				smem[y_s][threadIdx.x + i] = inputImage[y*(N)+x];
			}
		}
		else if (threadIdx.x == BLOCK_DIM - 1)
		{
			for (int i = 1; i <= R; ++i)
			{
				smem[y_s][x_s + i] = inputImage[y*(N)+x];
			}
		}
		if (threadIdx.y == 0)
		{
			for (int i = 0; i < R; ++i)
			{
				smem[threadIdx.y + i][x_s] = inputImage[y*(N)+x];
			}
		}
		else if (threadIdx.y == BLOCK_DIM - 1)
		{
			for (int i = 1; i <= R; ++i)
			{
				smem[y_s + i][x_s] = inputImage[y*(N)+x];
			}
		}

		__syncthreads();

		const int dim = D;
		int window[dim*dim];
		for (int j = -R; j <= R; ++j)
			for (int i = -R; i <= R; ++i)
				window[(R - j)*dim + R - i] = smem[y_s + j][x_s + i];
		for (int j = 0; j < (dim*dim + 1) / 2; ++j)
		{
			//   Find position of minimum element
			int min = j;
			for (int l = j + 1; l < dim*dim; ++l)
				if (window[l] < window[min])
					min = l;
			//   Put found minimum element in its place
			const int temp = window[j];
			window[j] = window[min];
			window[min] = temp;
		}
		/*printf("dim[%d][%d] | %-10d : \n\t\t\t%-10d, %-10d, %-10d, \n\t\t\t%-10d, %-10d, %-10d, \n\t\t\t%-10d, %-10d, %-10d \n", y, x, inputImage[y*(N) + x],
		smem[y_s - 1][x_s - 1], smem[y_s - 1][x_s], smem[y_s - 1][x_s + 1],
		smem[y_s][x_s - 1], smem[y_s][x_s], smem[y_s][x_s + 1],
		smem[y_s + 1][x_s - 1], smem[y_s + 1][x_s], smem[y_s + 1][x_s + 1]);*/

		printf("dim[%d][%d] | %-10d : \n%-10d, %-10d, %-10d, %-10d, %-10d, \n%-10d, %-10d, %-10d, %-10d, %-10d, \n%-10d, %-10d, %-10d, %-10d, %-10d, \n%-10d, %-10d, %-10d, %-10d, %-10d, \n%-10d, %-10d, %-10d, %-10d, %-10d,\n", y, x, inputImage[y*(N)+x],
			smem[y_s - 2][x_s - 2], smem[y_s - 2][x_s - 1], smem[y_s - 2][x_s], smem[y_s - 2][x_s + 1], smem[y_s - 2][x_s + 2],
			smem[y_s - 1][x_s - 2], smem[y_s - 1][x_s - 1], smem[y_s - 1][x_s], smem[y_s - 1][x_s + 1], smem[y_s - 1][x_s + 2],

			smem[y_s][x_s - 2], smem[y_s][x_s - 1], smem[y_s][x_s], smem[y_s][x_s + 1], smem[y_s][x_s + 2],

			smem[y_s + 1][x_s - 2], smem[y_s + 1][x_s - 1], smem[y_s + 1][x_s], smem[y_s + 1][x_s + 1], smem[y_s + 1][x_s + 2],
			smem[y_s + 2][x_s - 2], smem[y_s + 2][x_s - 1], smem[y_s + 2][x_s], smem[y_s + 2][x_s + 1], smem[y_s + 2][x_s + 2]
			);

		outputImage[y*width + x] = window[dim*dim / 2];
	}
	//how did you change your seat
}

__global__ void median_Kernel_21x21(int *inputImage, int *outputImage, int width, int height)//int *inputImage, int *outputImage, unsigned int width, unsigned int height)
{

	__shared__ int smem[BLOCK_W * BLOCK_H];

	int x = blockIdx.x*TILE_W + threadIdx.x - R;

	int y = blockIdx.y*TILE_H + threadIdx.y - R;

	int a = blockIdx.x*blockDim.x + threadIdx.x;
	int b = blockIdx.y*blockDim.y + threadIdx.y;

	//clamp to edge of image
	if (x < N + 2 && y < N + 2)
	{
		x = max(0, x);

		x = min(x, width - 1);

		y = max(y, 0);

		y = min(y, height - 1);



		unsigned int index = y*width + x;

		unsigned int bindex = threadIdx.y * blockDim.y + threadIdx.x;



		//each thread copies its pixel of the block to shared memory

		smem[bindex] = inputImage[index];

		__syncthreads();
		if ((threadIdx.x >= R) && (threadIdx.x < BLOCK_W - R) && (threadIdx.y >= R) && (threadIdx.y < BLOCK_W - R))

		{

			if (blockIdx.x*TILE_W + threadIdx.x < N + 2 && blockIdx.y*TILE_H + threadIdx.y < N + 2)
			{
				/*printf("[%d][%d] dim[%d][%d] | %-10d : \n%-10d, %-10d, %-10d, %-10d, %-10d, \n%-10d, %-10d, %-10d, %-10d, %-10d, \n%-10d, %-10d, %-10d, %-10d, %-10d, \n%-10d, %-10d, %-10d, %-10d, %-10d, \n%-10d, %-10d, %-10d, %-10d, %-10d,\n", b, a, y, x, inputImage[y*(N)+x],
					smem[bindex - 2 - BLOCK_W - BLOCK_W], smem[bindex - 1 - BLOCK_W - BLOCK_W], smem[bindex - BLOCK_W - BLOCK_W], smem[bindex + 1 - BLOCK_W - BLOCK_W], smem[bindex + 2 - BLOCK_W - BLOCK_W],
					smem[bindex - 2 - BLOCK_W], smem[bindex - 1 - BLOCK_W], smem[bindex - BLOCK_W], smem[bindex + 1 - BLOCK_W], smem[bindex + 2 - BLOCK_W],

					smem[bindex - 2], smem[bindex - 1], smem[bindex], smem[bindex + 1], smem[bindex + 2],

					smem[bindex - 2 + BLOCK_W], smem[bindex - 1 + BLOCK_W], smem[bindex + BLOCK_W], smem[bindex + 1 + BLOCK_W], smem[bindex + 2 + BLOCK_W],
					smem[bindex - 2 + BLOCK_W + BLOCK_W], smem[bindex - 1 + BLOCK_W + BLOCK_W], smem[bindex + BLOCK_W + BLOCK_W], smem[bindex + 1 + BLOCK_W + BLOCK_W], smem[bindex + 2 + BLOCK_W + BLOCK_W]
					);*/

				const int dim = D;
				int window[dim*dim];
				for (int j = -R; j <= R; ++j)
					for (int i = -R; i <= R; ++i)
						window[(R - j)*dim + R - i] = smem[bindex + j*BLOCK_W + i];
				for (int j = 0; j < (dim*dim + 1) / 2; ++j)
				{
					//   Find position of minimum element
					int min = j;
					for (int l = j + 1; l < dim*dim; ++l)
						if (window[l] < window[min])
							min = l;
					//   Put found minimum element in its place
					const int temp = window[j];
					window[j] = window[min];
					window[min] = temp;
				}
				outputImage[y*width + x] = window[dim*dim / 2];
				//outputImage[index] = smem[bindex];
			}

			//outputImage[index] = smem[bindex];
		}
	}
}

int main()
{
	char* name = "Points_[1.0e+08]_Noise_[030]_Normal.bin";
	double t1, t2;
	Grid grid_p = Grid(N, N);
	printf("CONFIGURATION\n------------------\n");
	printf("Bin dimensions : %d - %d\n", N, N);
	printf("Filter dimensions : %d - %d\n", D, D);
	printf("OpenMP threads: %d\n", omp_get_max_threads());
	printf("Binning file\n");
	ReadWrite::LoadData_omp(grid_p, name);
	printf("original grid : %d\n", grid_p.Count());
	//grid_p.Print();
	printf("padding grid\n");
	t1 = omp_get_wtime();
	Grid padded = grid_p.Pad(R);
	t2 = omp_get_wtime();
	printf("Time for padding: %12.3f sec.\n", t2 - t1);
	Grid input = grid_p;
	//padded.Print();
	//grid_p.Print();
	t1 = omp_get_wtime();
	//Filter::m_Filter_half(padded, grid_p, D);
	t2 = omp_get_wtime();
	//printf("Time for half Filtering: %12.3f sec.\n", t2 - t1);
	 
    // Add vectors in parallel.
    //cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
	Grid cuda = Grid(N, N);
	printf("filtered serial grid : %d\n", grid_p.Count());
	//grid_p.Print();
	//padded.Print();
	//grid_p.Print();
	 
	//cuda.set(0);
	cudaError_t cudaStatus = addWithCuda(cuda.grid, input.grid, input.cols*input.rows, cuda.rows*cuda.cols);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }
	printf("filtered cuda grid : %d\n", cuda.Count());
	//cuda.Print();

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *outputImage, const int *inputImage, unsigned int inputSize, unsigned int outputSize)
{
	const int size = 32;
	int a[size*size], b[size*size], c[size*size];
	for (int j = 0; j < size*size; ++j)
	{
			a[j] = 1;
			b[j] = 1;
			c[j] = -1;
			if (j%size == 0)
				c[j] = 0;
			
	}

	/*for (int j = 0; j < size; ++j)
	{
		printf("\n");
		for (int i = 0; i < size; ++i)
			printf("%-3d ", c[j * size + i]);
	}
	printf("\n");*/
    int *dev_i = 0;
    int *dev_o = 0;
	int *dev_a = 0;
	int *dev_b = 0;
	int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_i, inputSize * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	
	cudaStatus = cudaMalloc((void**)&dev_o, outputSize * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	
	/*cudaStatus = cudaMalloc((void**)&dev_c, 64 * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

	cudaStatus = cudaMalloc((void**)&dev_a, 64 * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

	cudaStatus = cudaMalloc((void**)&dev_b, 64 * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}*/

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
    // Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_i, inputImage, inputSize * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	
	/*cudaStatus = cudaMemcpy(dev_a, a, 64 * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

	cudaStatus = cudaMemcpy(dev_b, b, 64 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}*/

	dim3 grid((N + TILE_W - 1) / TILE_W, (N + TILE_H - 1) / TILE_W);
	dim3 threads(BLOCK_H, BLOCK_W);
	//dim3 threads(T, 1);
	//dim3 threads(BLOCK_DIM, BLOCK_DIM);
	//dim3 grid((N + threads.x - 1) / threads.x, (N + threads.y - 1) / threads.y);
	//dim3 grid((N + threads.x - 1) / threads.x, (N + threads.y - 1) / threads.y);
	//dim3 grid(1, 1);
	cudaEventRecord(start);
	//double t1, t2;
	//t1 = omp_get_wtime();
    // Launch a kernel on the GPU with one thread for each element.
	//if (R == 1)
	//	median_Kernel_3x3 << <grid, threads/*, 128*128*sizeof(int)*/ >> >(dev_i, dev_o, N, N);
	//if (R == 2)
	//	median_Kernel_5x5 << <grid, threads/*, 128*128*sizeof(int)*/ >> >(dev_i, dev_o, N, N);
	//
	//if (R == 10)
		median_Kernel_21x21 << <grid, threads/*, 128*128*sizeof(int)*/ >> >(dev_i, dev_o, N, N);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	printf("cudatime = %12.12f ms\n", milliseconds);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	//t2 = omp_get_wtime();
	//printf("Time for cuda Filtering: %12.3f sec.\n", t2 - t1);
    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }
    // Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(outputImage, dev_o, outputSize * sizeof(int), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	
	/*cudaStatus = cudaMemcpy(c, dev_c, 64 * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }*/

	
Error:
    cudaFree(dev_i);
    cudaFree(dev_o);
	cudaFree(dev_a);
	cudaFree(dev_b);
	cudaFree(dev_c);
	/*for (int j = 0; j < size; ++j)
	{
		printf("\n");
		for (int i = 0; i < size; ++i)
			printf("%-3d ", c[j * size + i]);
	}*/
    return cudaStatus;
}
