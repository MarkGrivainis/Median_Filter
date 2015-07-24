#ifndef __CUDACC__  
#define __CUDACC__
#endif


#define BLOCK_DIM 4
#define T 16

//for the kernels that uses fixed sizes
#define N 4096
#define M 4096
#define BLOCK_W 32
#define BLOCK_H 32
#define R 10
#define D (R*2+1)
#define D2 D*D
#define TILE_W (BLOCK_W - (2*R))
#define TILE_H (BLOCK_H - (2*R))

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <string>

#include "Grid.h"
#include "ReadWrite.h"
#include "Filter.h"
//texture<int, cudaTextureType1D, cudaReadModeElementType> tex;
cudaError_t median_filter(Grid &output, const Grid &input, const int filter_radius, int padding_type);

//__global__ void median_Kernel_border_padding(int *inputImage, int *outputImage, int width, int height)//int *inputImage, int *outputImage, unsigned int width, unsigned int height)
//{
//	__shared__ int smem[(BLOCK_DIM + 2 * R)][(BLOCK_DIM + 2 * R)];
//	int x_s = threadIdx.x + R;
//	int y_s = threadIdx.y + R;
//
//	int x = blockIdx.x*blockDim.x + threadIdx.x;
//	int y = blockIdx.y*blockDim.y + threadIdx.y;
//	if (x < N && y < N)
//	{
//		smem[y_s][x_s] = inputImage[y*(N)+x];
//
//		if (threadIdx.x == 0)
//		{
//			if (threadIdx.y == 0)
//			{
//				for (int j = -1; j >= -R; --j)
//					for (int i = -1; i >= -R; --i)
//					{
//						int offset_x = max(0, x + i);
//						int offset_y = max(0, y + j);
//						smem[y_s + j][x_s + i] = inputImage[offset_y*(N)+offset_x];
//					}
//
//			}
//			if (threadIdx.y == BLOCK_DIM - 1 || y == N - 1)
//			{
//				for (int j = 1; j <= R; ++j)
//					for (int i = -1; i >= -R; --i)
//					{
//						int offset_x = max(0, x + i);
//						int offset_y = min(y + j, height - 1);
//						smem[y_s + j][x_s + i] = inputImage[offset_y*(N)+offset_x];
//					}
//			}
//			for (int i = -1; i >= -R; --i)
//			{
//				int offset_x = max(0, x + i);
//				smem[y_s][x_s + i] = inputImage[y*(N)+offset_x];
//			}
//		}
//		if (threadIdx.x == BLOCK_DIM - 1 || x == N - 1)
//		{
//			if (threadIdx.y == 0)
//			{
//				for (int j = -1; j >= -R; --j)
//					for (int i = 1; i <= R; ++i)
//					{
//						int offset_x = min(x + i, width - 1);
//						int offset_y = max(0, y + j);
//						smem[y_s + j][x_s + i] = inputImage[offset_y*(N)+offset_x];
//					}
//			}
//			if (threadIdx.y == BLOCK_DIM - 1 || y == N - 1)
//			{
//				for (int j = 1; j <= R; ++j)
//					for (int i = 1; i <= R; ++i)
//					{
//						int offset_x = min(x + i, width - 1);
//						int offset_y = min(y + j, height - 1);
//						smem[y_s + j][x_s + i] = inputImage[offset_y*(N)+offset_x];
//					}
//			}
//			for (int i = 1; i <= R; ++i)
//			{
//				int offset_x = min(x + i, width - 1);
//				smem[y_s][x_s + i] = inputImage[y*(N)+offset_x];
//			}
//		}
//		if (threadIdx.y == 0)
//		{
//			for (int i = -1; i >= -R; --i)
//			{
//				int offset_y = max(0, y + i);
//				smem[y_s + i][x_s] = inputImage[offset_y*(N)+x];
//			}
//		}
//		if (threadIdx.y == BLOCK_DIM - 1 || y == N - 1)
//		{
//			for (int i = 1; i <= R; ++i)
//			{
//				int offset_y = min(y + i, height - 1);
//				smem[y_s + i][x_s] = inputImage[offset_y*(N)+x];
//			}
//		}
//
//		__syncthreads();
//		const int dim = D;
//		const int k = dim*dim / 2;
//		int a = 0;
//		unsigned int window[dim*dim];
//
//
//		for (int h_offset = -R; h_offset <= R; h_offset++)
//			for (int w_offset = -R; w_offset <= R; w_offset++)
//				window[a++] = smem[y_s + h_offset][x_s + w_offset];
//
//		int low, high;
//		int median;
//		int middle, ll, hh;
//
//		low = 0; high = dim*dim - 1; median = (low + high) / 2;
//		for (;;)
//		{
//			if (high <= low)
//				break;
//			if (high == low + 1)
//			{
//				if (window[low] > window[high])
//				{
//					int temp = window[low];
//					window[low] = window[high];
//					window[high] = temp;
//				}
//				break;
//			}
//
//			middle = (low + high) / 2;
//			if (window[middle] > window[high])
//			{
//				int temp = window[middle];
//				window[middle] = window[high];
//				window[high] = temp;
//			}
//			if (window[low] > window[high])
//			{
//				int temp = window[low];
//				window[low] = window[high];
//				window[high] = temp;
//			}
//			if (window[middle] > window[low])
//			{
//				int temp = window[low];
//				window[low] = window[middle];
//				window[middle] = temp;
//			}
//			int temp = window[low + 1];
//			window[low + 1] = window[middle];
//			window[middle] = temp;
//
//			ll = low + 1;
//			hh = high;
//			for (;;)
//			{
//				do{
//					ll++;
//				} while (window[low] > window[ll]);
//				do
//				{
//					hh--;
//				} while (window[hh] > window[low]);
//				if (hh < ll)
//					break;
//				int temp = window[ll];
//				window[ll] = window[hh];
//				window[hh] = temp;
//			}
//
//			temp = window[low];
//			window[low] = window[hh];
//			window[hh] = temp;
//
//			if (hh <= median)
//				low = ll;
//			if (hh >= median)
//				high = hh - 1;
//		}
//		outputImage[y*width + x] = window[median];
//
//	}
//}

__global__ void median_Kernel_21(int* inputImage, int* outputImage, const int width, const int height, const int radius)
{
	const int dim = radius * 2 + 1;
	volatile __shared__ int smem[(1 + 2 * 10)][(T + 2 * 10)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y;
	int x_p = x + radius;
	int x_o = threadIdx.x + radius;
	int y_o = radius;
	if (x >= width)
		return;

	if (threadIdx.x == 0)
		for (int j = 0; j < (1 + 2 * radius); ++j)
			for (int i = 0; i < (T + 2 * radius); ++i)
				smem[j][i] = -1;

	__syncthreads();

	for (int j = -radius; j <= radius; ++j)
		for (int i = -radius; i <= radius; ++i)
			if (smem[y_o + j][x_o + i] == -1)
				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];
	
	//find the offset of the current window
	//int a =  (threadIdx.x)*dim*dim;
	int a = 0;
	int window[21*21];
	//populate the window with the values around the current index
	for (int h_offset = -radius; h_offset <= radius; h_offset++)
		for (int w_offset = -radius; w_offset <= radius; w_offset++)
			window[a++] = smem[y_o + h_offset][ x_o + w_offset];
	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	int low, high;
	int median;
	int middle, ll, hh;
	low = 0; high = dim*dim - 1; median = (low + high) / 2;
	for (;;)
	{
		if (high <= low)
			break;
		if (high == low + 1)
		{
			if (window[low] > window[high])
			{
				int temp = window[low];
				window[low] = window[high];
				window[high] = temp;
			}
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		if (window[middle] > window[high])
		{
			int temp = window[middle];
			window[middle] = window[high];
			window[high] = temp;
		}
		if (window[low] > window[high])
		{
			int temp = window[low];
			window[low] = window[high];
			window[high] = temp;
		}
		if (window[middle] > window[low])
		{
			int temp = window[low];
			window[low] = window[middle];
			window[middle] = temp;
		}
		//Swap low item (now in position middle) into position (low+1)
		int temp = window[low + 1];
		window[low + 1] = window[middle];
		window[middle] = temp;

		//Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for (;;)
		{
			do{
				ll++;
			} while (window[low] > window[ll]);
			do
			{
				hh--;
			} while (window[hh] > window[low]);
			if (hh < ll)
				break;
			int temp = window[ll];
			window[ll] = window[hh];
			window[hh] = temp;
		}

		//Swap middle item (in position low) back into correct position
		temp = window[low];
		window[low] = window[hh];
		window[hh] = temp;
		//Re-set active partition
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}
	
	outputImage[y*width + x] = window[median];
}

__global__ void median_Kernel_19(int* inputImage, int* outputImage, const int width, const int height, const int radius)
{
	const int dim = radius * 2 + 1;
	__shared__ int smem[(1 + 2 * 9)][(T + 2 * 9)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y;
	int x_p = x + radius;
	int x_o = threadIdx.x + radius;
	int y_o = radius;
	if (x >= width)
		return;

	if (threadIdx.x == 0)
		for (int j = 0; j < (1 + 2 * radius); ++j)
			for (int i = 0; i < (T + 2 * radius); ++i)
				smem[j][i] = 0;

	__syncthreads();

	for (int j = -radius; j <= radius; ++j)
		for (int i = -radius; i <= radius; ++i)
			if (smem[y_o + j][x_o + i] == 0)
				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];

	__syncthreads();

	//find the offset of the current window
	//int a =  (threadIdx.x)*dim*dim;
	int a = 0;
	int window[D*D];
	//populate the window with the values around the current index
	for (int h_offset = -radius; h_offset <= radius; h_offset++)
		for (int w_offset = -radius; w_offset <= radius; w_offset++)
			window[a++] = smem[y_o + h_offset][x_o + w_offset];

	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	int low, high;
	int median;
	int middle, ll, hh;
	low = 0; high = dim*dim - 1; median = (low + high) / 2;
	for (;;)
	{
		if (high <= low)
			break;
		if (high == low + 1)
		{
			if (window[low] > window[high])
			{
				int temp = window[low];
				window[low] = window[high];
				window[high] = temp;
			}
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		if (window[middle] > window[high])
		{
			int temp = window[middle];
			window[middle] = window[high];
			window[high] = temp;
		}
		if (window[low] > window[high])
		{
			int temp = window[low];
			window[low] = window[high];
			window[high] = temp;
		}
		if (window[middle] > window[low])
		{
			int temp = window[low];
			window[low] = window[middle];
			window[middle] = temp;
		}
		//Swap low item (now in position middle) into position (low+1)
		int temp = window[low + 1];
		window[low + 1] = window[middle];
		window[middle] = temp;

		//Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for (;;)
		{
			do{
				ll++;
			} while (window[low] > window[ll]);
			do
			{
				hh--;
			} while (window[hh] > window[low]);
			if (hh < ll)
				break;
			int temp = window[ll];
			window[ll] = window[hh];
			window[hh] = temp;
		}

		//Swap middle item (in position low) back into correct position
		temp = window[low];
		window[low] = window[hh];
		window[hh] = temp;
		//Re-set active partition
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}

	outputImage[y*width + x] = window[median];
}

__global__ void median_Kernel_17(int* inputImage, int* outputImage, const int width, const int height, const int radius)
{
	const int dim = radius * 2 + 1;
	__shared__ int smem[(1 + 2 * 8)][(T + 2 * 8)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y;
	int x_p = x + radius;
	int x_o = threadIdx.x + radius;
	int y_o = radius;
	if (x >= width)
		return;

	if (threadIdx.x == 0)
		for (int j = 0; j < (1 + 2 * radius); ++j)
			for (int i = 0; i < (T + 2 * radius); ++i)
				smem[j][i] = 0;

	__syncthreads();

	for (int j = -radius; j <= radius; ++j)
		for (int i = -radius; i <= radius; ++i)
			if (smem[y_o + j][x_o + i] == 0)
				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];

	__syncthreads();

	//find the offset of the current window
	//int a =  (threadIdx.x)*dim*dim;
	int a = 0;
	int window[19*19];
	//populate the window with the values around the current index
	for (int h_offset = -radius; h_offset <= radius; h_offset++)
		for (int w_offset = -radius; w_offset <= radius; w_offset++)
			window[a++] = smem[y_o + h_offset][x_o + w_offset];

	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	int low, high;
	int median;
	int middle, ll, hh;
	low = 0; high = dim*dim - 1; median = (low + high) / 2;
	for (;;)
	{
		if (high <= low)
			break;
		if (high == low + 1)
		{
			if (window[low] > window[high])
			{
				int temp = window[low];
				window[low] = window[high];
				window[high] = temp;
			}
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		if (window[middle] > window[high])
		{
			int temp = window[middle];
			window[middle] = window[high];
			window[high] = temp;
		}
		if (window[low] > window[high])
		{
			int temp = window[low];
			window[low] = window[high];
			window[high] = temp;
		}
		if (window[middle] > window[low])
		{
			int temp = window[low];
			window[low] = window[middle];
			window[middle] = temp;
		}
		//Swap low item (now in position middle) into position (low+1)
		int temp = window[low + 1];
		window[low + 1] = window[middle];
		window[middle] = temp;

		//Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for (;;)
		{
			do{
				ll++;
			} while (window[low] > window[ll]);
			do
			{
				hh--;
			} while (window[hh] > window[low]);
			if (hh < ll)
				break;
			int temp = window[ll];
			window[ll] = window[hh];
			window[hh] = temp;
		}

		//Swap middle item (in position low) back into correct position
		temp = window[low];
		window[low] = window[hh];
		window[hh] = temp;
		//Re-set active partition
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}

	outputImage[y*width + x] = window[median];
}

__global__ void median_Kernel_15(int* inputImage, int* outputImage, const int width, const int height, const int radius)
{
	const int dim = radius * 2 + 1;
	__shared__ int smem[(1 + 2 * 7)][(T + 2 * 7)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y;
	int x_p = x + radius;
	int x_o = threadIdx.x + radius;
	int y_o = radius;
	if (x >= width)
		return;

	if (threadIdx.x == 0)
		for (int j = 0; j < (1 + 2 * radius); ++j)
			for (int i = 0; i < (T + 2 * radius); ++i)
				smem[j][i] = 0;

	__syncthreads();

	for (int j = -radius; j <= radius; ++j)
		for (int i = -radius; i <= radius; ++i)
			if (smem[y_o + j][x_o + i] == 0)
				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];

	__syncthreads();

	//find the offset of the current window
	//int a =  (threadIdx.x)*dim*dim;
	int a = 0;
	int window[15*15];
	//populate the window with the values around the current index
	for (int h_offset = -radius; h_offset <= radius; h_offset++)
		for (int w_offset = -radius; w_offset <= radius; w_offset++)
			window[a++] = smem[y_o + h_offset][x_o + w_offset];

	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	int low, high;
	int median;
	int middle, ll, hh;
	low = 0; high = dim*dim - 1; median = (low + high) / 2;
	for (;;)
	{
		if (high <= low)
			break;
		if (high == low + 1)
		{
			if (window[low] > window[high])
			{
				int temp = window[low];
				window[low] = window[high];
				window[high] = temp;
			}
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		if (window[middle] > window[high])
		{
			int temp = window[middle];
			window[middle] = window[high];
			window[high] = temp;
		}
		if (window[low] > window[high])
		{
			int temp = window[low];
			window[low] = window[high];
			window[high] = temp;
		}
		if (window[middle] > window[low])
		{
			int temp = window[low];
			window[low] = window[middle];
			window[middle] = temp;
		}
		//Swap low item (now in position middle) into position (low+1)
		int temp = window[low + 1];
		window[low + 1] = window[middle];
		window[middle] = temp;

		//Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for (;;)
		{
			do{
				ll++;
			} while (window[low] > window[ll]);
			do
			{
				hh--;
			} while (window[hh] > window[low]);
			if (hh < ll)
				break;
			int temp = window[ll];
			window[ll] = window[hh];
			window[hh] = temp;
		}

		//Swap middle item (in position low) back into correct position
		temp = window[low];
		window[low] = window[hh];
		window[hh] = temp;
		//Re-set active partition
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}

	outputImage[y*width + x] = window[median];
}

__global__ void median_Kernel_13(int* inputImage, int* outputImage, const int width, const int height, const int radius)
{
	const int dim = radius * 2 + 1;
	__shared__ int smem[(1 + 2 * 6)][(T + 2 * 6)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y;
	int x_p = x + radius;
	int x_o = threadIdx.x + radius;
	int y_o = radius;
	if (x >= width)
		return;

	if (threadIdx.x == 0)
		for (int j = 0; j < (1 + 2 * radius); ++j)
			for (int i = 0; i < (T + 2 * radius); ++i)
				smem[j][i] = 0;

	__syncthreads();

	for (int j = -radius; j <= radius; ++j)
		for (int i = -radius; i <= radius; ++i)
			if (smem[y_o + j][x_o + i] == 0)
				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];

	__syncthreads();

	//find the offset of the current window
	//int a =  (threadIdx.x)*dim*dim;
	int a = 0;
	int window[13*13];
	//populate the window with the values around the current index
	for (int h_offset = -radius; h_offset <= radius; h_offset++)
		for (int w_offset = -radius; w_offset <= radius; w_offset++)
			window[a++] = smem[y_o + h_offset][x_o + w_offset];

	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	int low, high;
	int median;
	int middle, ll, hh;
	low = 0; high = dim*dim - 1; median = (low + high) / 2;
	for (;;)
	{
		if (high <= low)
			break;
		if (high == low + 1)
		{
			if (window[low] > window[high])
			{
				int temp = window[low];
				window[low] = window[high];
				window[high] = temp;
			}
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		if (window[middle] > window[high])
		{
			int temp = window[middle];
			window[middle] = window[high];
			window[high] = temp;
		}
		if (window[low] > window[high])
		{
			int temp = window[low];
			window[low] = window[high];
			window[high] = temp;
		}
		if (window[middle] > window[low])
		{
			int temp = window[low];
			window[low] = window[middle];
			window[middle] = temp;
		}
		//Swap low item (now in position middle) into position (low+1)
		int temp = window[low + 1];
		window[low + 1] = window[middle];
		window[middle] = temp;

		//Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for (;;)
		{
			do{
				ll++;
			} while (window[low] > window[ll]);
			do
			{
				hh--;
			} while (window[hh] > window[low]);
			if (hh < ll)
				break;
			int temp = window[ll];
			window[ll] = window[hh];
			window[hh] = temp;
		}

		//Swap middle item (in position low) back into correct position
		temp = window[low];
		window[low] = window[hh];
		window[hh] = temp;
		//Re-set active partition
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}

	outputImage[y*width + x] = window[median];
}

__global__ void median_Kernel_11(int* inputImage, int* outputImage, const int width, const int height, const int radius)
{
	const int dim = radius * 2 + 1;
	__shared__ int smem[(1 + 2 * 5)][(T + 2 * 5)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y;
	int x_p = x + radius;
	int x_o = threadIdx.x + radius;
	int y_o = radius;
	if (x >= width)
		return;

	if (threadIdx.x == 0)
		for (int j = 0; j < (1 + 2 * radius); ++j)
			for (int i = 0; i < (T + 2 * radius); ++i)
				smem[j][i] = 0;

	__syncthreads();

	for (int j = -radius; j <= radius; ++j)
		for (int i = -radius; i <= radius; ++i)
			if (smem[y_o + j][x_o + i] == 0)
				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];

	__syncthreads();

	//find the offset of the current window
	//int a =  (threadIdx.x)*dim*dim;
	int a = 0;
	int window[11*11];
	//populate the window with the values around the current index
	for (int h_offset = -radius; h_offset <= radius; h_offset++)
		for (int w_offset = -radius; w_offset <= radius; w_offset++)
			window[a++] = smem[y_o + h_offset][x_o + w_offset];

	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	int low, high;
	int median;
	int middle, ll, hh;
	low = 0; high = dim*dim - 1; median = (low + high) / 2;
	for (;;)
	{
		if (high <= low)
			break;
		if (high == low + 1)
		{
			if (window[low] > window[high])
			{
				int temp = window[low];
				window[low] = window[high];
				window[high] = temp;
			}
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		if (window[middle] > window[high])
		{
			int temp = window[middle];
			window[middle] = window[high];
			window[high] = temp;
		}
		if (window[low] > window[high])
		{
			int temp = window[low];
			window[low] = window[high];
			window[high] = temp;
		}
		if (window[middle] > window[low])
		{
			int temp = window[low];
			window[low] = window[middle];
			window[middle] = temp;
		}
		//Swap low item (now in position middle) into position (low+1)
		int temp = window[low + 1];
		window[low + 1] = window[middle];
		window[middle] = temp;

		//Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for (;;)
		{
			do{
				ll++;
			} while (window[low] > window[ll]);
			do
			{
				hh--;
			} while (window[hh] > window[low]);
			if (hh < ll)
				break;
			int temp = window[ll];
			window[ll] = window[hh];
			window[hh] = temp;
		}

		//Swap middle item (in position low) back into correct position
		temp = window[low];
		window[low] = window[hh];
		window[hh] = temp;
		//Re-set active partition
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}

	outputImage[y*width + x] = window[median];
}

__global__ void median_Kernel_9(int* inputImage, int* outputImage, const int width, const int height, const int radius)
{
	const int dim = radius * 2 + 1;
	__shared__ int smem[(1 + 2 * 4)][(T + 2 * 4)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y;
	int x_p = x + radius;
	int x_o = threadIdx.x + radius;
	int y_o = radius;
	if (x >= width)
		return;

	if (threadIdx.x == 0)
		for (int j = 0; j < (1 + 2 * radius); ++j)
			for (int i = 0; i < (T + 2 * radius); ++i)
				smem[j][i] = 0;

	__syncthreads();

	for (int j = -radius; j <= radius; ++j)
		for (int i = -radius; i <= radius; ++i)
			if (smem[y_o + j][x_o + i] == 0)
				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];

	__syncthreads();

	//find the offset of the current window
	//int a =  (threadIdx.x)*dim*dim;
	int a = 0;
	int window[9*9];
	//populate the window with the values around the current index
	for (int h_offset = -radius; h_offset <= radius; h_offset++)
		for (int w_offset = -radius; w_offset <= radius; w_offset++)
			window[a++] = smem[y_o + h_offset][x_o + w_offset];

	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	int low, high;
	int median;
	int middle, ll, hh;
	low = 0; high = dim*dim - 1; median = (low + high) / 2;
	for (;;)
	{
		if (high <= low)
			break;
		if (high == low + 1)
		{
			if (window[low] > window[high])
			{
				int temp = window[low];
				window[low] = window[high];
				window[high] = temp;
			}
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		if (window[middle] > window[high])
		{
			int temp = window[middle];
			window[middle] = window[high];
			window[high] = temp;
		}
		if (window[low] > window[high])
		{
			int temp = window[low];
			window[low] = window[high];
			window[high] = temp;
		}
		if (window[middle] > window[low])
		{
			int temp = window[low];
			window[low] = window[middle];
			window[middle] = temp;
		}
		//Swap low item (now in position middle) into position (low+1)
		int temp = window[low + 1];
		window[low + 1] = window[middle];
		window[middle] = temp;

		//Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for (;;)
		{
			do{
				ll++;
			} while (window[low] > window[ll]);
			do
			{
				hh--;
			} while (window[hh] > window[low]);
			if (hh < ll)
				break;
			int temp = window[ll];
			window[ll] = window[hh];
			window[hh] = temp;
		}

		//Swap middle item (in position low) back into correct position
		temp = window[low];
		window[low] = window[hh];
		window[hh] = temp;
		//Re-set active partition
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}

	outputImage[y*width + x] = window[median];
}

__global__ void median_Kernel_7(int* inputImage, int* outputImage, const int width, const int height, const int radius)
{
	const int dim = radius * 2 + 1;
	__shared__ int smem[(1 + 2 * 3)][(T + 2 * 3)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y;
	int x_p = x + radius;
	int x_o = threadIdx.x + radius;
	int y_o = radius;
	if (x >= width)
		return;

	if (threadIdx.x == 0)
		for (int j = 0; j < (1 + 2 * radius); ++j)
			for (int i = 0; i < (T + 2 * radius); ++i)
				smem[j][i] = 0;

	__syncthreads();

	for (int j = -radius; j <= radius; ++j)
		for (int i = -radius; i <= radius; ++i)
			if (smem[y_o + j][x_o + i] == 0)
				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];

	__syncthreads();

	//find the offset of the current window
	//int a =  (threadIdx.x)*dim*dim;
	int a = 0;
	int window[7*7];
	//populate the window with the values around the current index
	for (int h_offset = -radius; h_offset <= radius; h_offset++)
		for (int w_offset = -radius; w_offset <= radius; w_offset++)
			window[a++] = smem[y_o + h_offset][x_o + w_offset];

	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	int low, high;
	int median;
	int middle, ll, hh;
	low = 0; high = dim*dim - 1; median = (low + high) / 2;
	for (;;)
	{
		if (high <= low)
			break;
		if (high == low + 1)
		{
			if (window[low] > window[high])
			{
				int temp = window[low];
				window[low] = window[high];
				window[high] = temp;
			}
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		if (window[middle] > window[high])
		{
			int temp = window[middle];
			window[middle] = window[high];
			window[high] = temp;
		}
		if (window[low] > window[high])
		{
			int temp = window[low];
			window[low] = window[high];
			window[high] = temp;
		}
		if (window[middle] > window[low])
		{
			int temp = window[low];
			window[low] = window[middle];
			window[middle] = temp;
		}
		//Swap low item (now in position middle) into position (low+1)
		int temp = window[low + 1];
		window[low + 1] = window[middle];
		window[middle] = temp;

		//Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for (;;)
		{
			do{
				ll++;
			} while (window[low] > window[ll]);
			do
			{
				hh--;
			} while (window[hh] > window[low]);
			if (hh < ll)
				break;
			int temp = window[ll];
			window[ll] = window[hh];
			window[hh] = temp;
		}

		//Swap middle item (in position low) back into correct position
		temp = window[low];
		window[low] = window[hh];
		window[hh] = temp;
		//Re-set active partition
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}

	outputImage[y*width + x] = window[median];
}

__global__ void median_Kernel_5(int* inputImage, int* outputImage, const int width, const int height, const int radius)
{
	const int dim = radius * 2 + 1;
	__shared__ int smem[(1 + 2 * 2)][(T + 2 * 2)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y;
	int x_p = x + radius;
	int x_o = threadIdx.x + radius;
	int y_o = radius;
	if (x >= width)
		return;

	if (threadIdx.x == 0)
		for (int j = 0; j < (1 + 2 * radius); ++j)
			for (int i = 0; i < (T + 2 * radius); ++i)
				smem[j][i] = 0;

	__syncthreads();

	for (int j = -radius; j <= radius; ++j)
		for (int i = -radius; i <= radius; ++i)
			if (smem[y_o + j][x_o + i] == 0)
				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];

	__syncthreads();

	//find the offset of the current window
	//int a =  (threadIdx.x)*dim*dim;
	int a = 0;
	int window[5*5];
	//populate the window with the values around the current index
	for (int h_offset = -radius; h_offset <= radius; h_offset++)
		for (int w_offset = -radius; w_offset <= radius; w_offset++)
			window[a++] = smem[y_o + h_offset][x_o + w_offset];

	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	int low, high;
	int median;
	int middle, ll, hh;
	low = 0; high = dim*dim - 1; median = (low + high) / 2;
	for (;;)
	{
		if (high <= low)
			break;
		if (high == low + 1)
		{
			if (window[low] > window[high])
			{
				int temp = window[low];
				window[low] = window[high];
				window[high] = temp;
			}
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		if (window[middle] > window[high])
		{
			int temp = window[middle];
			window[middle] = window[high];
			window[high] = temp;
		}
		if (window[low] > window[high])
		{
			int temp = window[low];
			window[low] = window[high];
			window[high] = temp;
		}
		if (window[middle] > window[low])
		{
			int temp = window[low];
			window[low] = window[middle];
			window[middle] = temp;
		}
		//Swap low item (now in position middle) into position (low+1)
		int temp = window[low + 1];
		window[low + 1] = window[middle];
		window[middle] = temp;

		//Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for (;;)
		{
			do{
				ll++;
			} while (window[low] > window[ll]);
			do
			{
				hh--;
			} while (window[hh] > window[low]);
			if (hh < ll)
				break;
			int temp = window[ll];
			window[ll] = window[hh];
			window[hh] = temp;
		}

		//Swap middle item (in position low) back into correct position
		temp = window[low];
		window[low] = window[hh];
		window[hh] = temp;
		//Re-set active partition
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}

	outputImage[y*width + x] = window[median];
}

__global__ void median_Kernel_3(int* inputImage, int* outputImage, const int width, const int height, const int radius)
{
	const int dim = radius * 2 + 1;
	__shared__ int smem[(1 + 2 * 1)][(T + 2 * 1)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y;
	int x_p = x + radius;
	int x_o = threadIdx.x + radius;
	int y_o = radius;
	if (x >= width)
		return;

	if (threadIdx.x == 0)
		for (int j = 0; j < (1 + 2 * radius); ++j)
			for (int i = 0; i < (T + 2 * radius); ++i)
				smem[j][i] = 0;

	__syncthreads();


	//smem[0][x_o - 1] = inputImage[(0)*(width + 2 * radius) + x_p - 1]; smem[0][x_o] = inputImage[(0)*(width + 2 * radius) + x_p]; smem[0][x_o + 1] = inputImage[(0)*(width + 2 * radius) + x_p + 1];
	//smem[1][x_o - 1] = inputImage[(1)*(width + 2 * radius) + x_p - 1]; smem[1][x_o] = inputImage[(1)*(width + 2 * radius) + x_p]; smem[1][x_o + 1] = inputImage[(1)*(width + 2 * radius) + x_p + 1];
	//smem[2][x_o - 1] = inputImage[(2)*(width + 2 * radius) + x_p - 1]; smem[2][x_o] = inputImage[(2)*(width + 2 * radius) + x_p]; smem[2][x_o + 1] = inputImage[(2)*(width + 2 * radius) + x_p + 1];
	for (int j = -radius; j <= radius; ++j)
		for (int i = -radius; i <= radius; ++i)
			if (smem[y_o + j][x_o + i] == 0)
				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];

	__syncthreads();

	//find the offset of the current window
	//int a =  (threadIdx.x)*dim*dim;
	int a = 0;
	int window[3*3];
	//populate the window with the values around the current index
//	window[0] = smem[y_o-1][x_o-1]; window[1] = smem[y_o-1][x_o]; window[2] = smem[y_o-1][x_o+1];
//	window[3] = smem[y_o][x_o - 1]; window[4] = smem[y_o][x_o]; window[5] = smem[y_o][x_o + 1];
//	window[6] = smem[y_o + 1][x_o - 1]; window[7] = smem[y_o + 1][x_o]; window[8] = smem[y_o + 1][x_o + 1];
	for (int h_offset = -radius; h_offset <= radius; h_offset++)
		for (int w_offset = -radius; w_offset <= radius; w_offset++)
			window[a++] = smem[y_o + h_offset][x_o + w_offset];

	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	int low, high;
	int median;
	int middle, ll, hh;
	low = 0; high = dim*dim - 1; median = (low + high) / 2;
	for (;;)
	{
		if (high <= low)
			break;
		if (high == low + 1)
		{
			if (window[low] > window[high])
			{
				int temp = window[low];
				window[low] = window[high];
				window[high] = temp;
			}
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		if (window[middle] > window[high])
		{
			int temp = window[middle];
			window[middle] = window[high];
			window[high] = temp;
		}
		if (window[low] > window[high])
		{
			int temp = window[low];
			window[low] = window[high];
			window[high] = temp;
		}
		if (window[middle] > window[low])
		{
			int temp = window[low];
			window[low] = window[middle];
			window[middle] = temp;
		}
		//Swap low item (now in position middle) into position (low+1)
		int temp = window[low + 1];
		window[low + 1] = window[middle];
		window[middle] = temp;

		//Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for (;;)
		{
			do{
				ll++;
			} while (window[low] > window[ll]);
			do
			{
				hh--;
			} while (window[hh] > window[low]);
			if (hh < ll)
				break;
			int temp = window[ll];
			window[ll] = window[hh];
			window[hh] = temp;
		}

		//Swap middle item (in position low) back into correct position
		temp = window[low];
		window[low] = window[hh];
		window[hh] = temp;
		//Re-set active partition
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}

	outputImage[y*width + x] = window[median];
}

__global__ void median_Kernel_border_prepadded(int *inputImage, int *outputImage, const int width, const int height, const int radius)//int *inputImage, int *outputImage, unsigned int width, unsigned int height)
{
	const int dim = radius * 2 + 1;
	//extern __shared__ int shared[];
	volatile __shared__ int smem[(BLOCK_DIM + 2 * R)*(BLOCK_DIM + 2 * R)];
	//int *smem = &shared[0];
	//int *window = &shared[(blockDim.x + 2 * radius)*(blockDim.y + 2 * radius)];
	int x_s = threadIdx.x + radius;
	int y_s = threadIdx.y + radius;

	int x_a = blockIdx.x*blockDim.x + threadIdx.x + radius;
	int y_a = blockIdx.y*blockDim.y + threadIdx.y + radius;

	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y*blockDim.y + threadIdx.y;
	if (x >= width || y >= height)
		return;
	//	if (smem[y_s*(BLOCK_DIM + radius * 2) + x_s] == 0)
			smem[y_s*(BLOCK_DIM + radius * 2) + x_s] = inputImage[y_a*(width + 2 * radius) + x_a];

		//pad the left column
		if (threadIdx.x == 0)
		{
			//if it is the first row pad the top left corner
			if (threadIdx.y == 0)
			{
				for (int j = 0; j < radius; ++j)
					for (int i = 0; i < radius; ++i)
					//	if (smem[(threadIdx.y + j)*(BLOCK_DIM + radius * 2) + threadIdx.x + i] == 0)
							smem[(threadIdx.y + j)*(BLOCK_DIM + radius * 2) + threadIdx.x + i] = inputImage[(y + j)*(width + 2 * radius) + x + i];
			}
			//if it the last row pad the bottow left corner
			if (threadIdx.y == BLOCK_DIM - 1 || y == width - 1)
			{
				for (int j = 0; j < radius; ++j)
					for (int i = 0; i < radius; ++i)
					//	if (smem[(y_s + 1 + j)*(BLOCK_DIM + radius * 2) + threadIdx.x + i] == 0)
							smem[(y_s + 1 + j)*(BLOCK_DIM + radius * 2) + threadIdx.x + i] = inputImage[(y_a + 1 + j)*(width + 2 * radius) + x + i];
			}
			//pad to the left of the column
			for (int i = 0; i < radius; ++i)
			{
				//if (smem[y_s*(BLOCK_DIM + radius * 2) + threadIdx.x + i] == 0)
					smem[y_s*(BLOCK_DIM + radius * 2) + threadIdx.x + i] = inputImage[y_a*(width + 2 * radius) + x + i];
			}
		}
		//pad the right column
		if (threadIdx.x == BLOCK_DIM - 1 || x == width - 1)
		{
			//if it is the first row pad the top right corner
			if (threadIdx.y == 0)
			{
				for (int j = 0; j < radius; ++j)
					for (int i = 0; i < radius; ++i)
						//if (smem[(threadIdx.y + j)*(BLOCK_DIM + radius * 2) + x_s + 1 + i] == 0)
							smem[(threadIdx.y + j)*(BLOCK_DIM + radius * 2) + x_s + 1 + i] = inputImage[(y + j)*(width + 2 * radius) + x_a + 1 + i];
			}
			//if it is the last row pad the bottom right corner
			if (threadIdx.y == BLOCK_DIM - 1 || y == width - 1)
			{
				for (int j = 0; j < radius; ++j)
					for (int i = 0; i < radius; ++i)
						//if (smem[(y_s + 1 + j)*(BLOCK_DIM + radius * 2) + x_s + 1 + i] == 0)
							smem[(y_s + 1 + j)*(BLOCK_DIM + radius * 2) + x_s + 1 + i] = inputImage[(y_a + 1 + j)*(width + 2 * radius) + x_a + 1 + i];
			}
			//otherwise pad to the right
			for (int i = 0; i < radius; ++i)
			{
				//if (smem[y_s*(BLOCK_DIM + radius * 2) + x_s + 1 + i] == 0)
					smem[y_s*(BLOCK_DIM + radius * 2) + x_s + 1 + i] = inputImage[y_a*(width + 2 * radius) + x_a + 1 + i];
			}
		}
		//pad the top row
		if (threadIdx.y == 0)
		{
			//pad above the current index
			for (int i = 0; i < radius; ++i)
			{
			//	if (smem[(threadIdx.y + i)*(BLOCK_DIM + radius * 2) + x_s] == 0)
					smem[(threadIdx.y + i)*(BLOCK_DIM + radius * 2) + x_s] = inputImage[(y + i)*(width + 2 * radius) + x_a];
			}
		}
		//pad the bottom row
		if (threadIdx.y == BLOCK_DIM - 1 || y == width - 1)
		{
			for (int i = 0; i < radius; ++i)
			{
				//pad the below the index
				//if (smem[(y_s + 1 + i)*(BLOCK_DIM + radius * 2) + x_s] == 0)
					smem[(y_s + 1 + i)*(BLOCK_DIM + radius * 2) + x_s] = inputImage[(y_a + 1 + i)*(width + 2 * radius) + x_a];
			}
		}

		__syncthreads();
		//find the offset of the current window
		//int a = (threadIdx.x + threadIdx.y*BLOCK_DIM)*dim*dim;
		int window[D*D];
		int a = 0;
		//populate the window with the values around the current index
		for (int h_offset = -radius; h_offset <= radius; h_offset++)
			for (int w_offset = -radius; w_offset <= radius; w_offset++)
				window[a++] = smem[(y_s + h_offset) * (BLOCK_DIM + 2 * radius) + x_s + w_offset];

		int low, high;
		int median;
		int middle, ll, hh;
		/*
		* This Quickselect routine is based on the algorithm described in
		* "Numerical recipes in C", Second Edition,
		* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
		* This code by Nicolas Devillard - 1998. Public domain.
		*/
		//low = (threadIdx.x + threadIdx.y*BLOCK_DIM)*dim*dim; high = (threadIdx.x + threadIdx.y*BLOCK_DIM)*dim*dim+(dim*dim - 1); median = (low + high) / 2;
		low = 0; high = D*D - 1; median = (low + high) / 2;
		//for (;;)
		//{
		//	if (high <= low)
		//		break;
		//	if (high == low + 1)
		//	{
		//		if (window[low] > window[high])
		//		{
		//			int temp = window[low];
		//			window[low] = window[high];
		//			window[high] = temp;
		//		}
		//		break;
		//	}
		//	// Find median of low, middle and high items; swap into position low
		//	middle = (low + high) / 2;
		//	if (window[middle] > window[high])
		//	{
		//		int temp = window[middle];
		//		window[middle] = window[high];
		//		window[high] = temp;
		//	}
		//	if (window[low] > window[high])
		//	{
		//		int temp = window[low];
		//		window[low] = window[high];
		//		window[high] = temp;
		//	}
		//	if (window[middle] > window[low])
		//	{
		//		int temp = window[low];
		//		window[low] = window[middle];
		//		window[middle] = temp;
		//	}
		//	//Swap low item (now in position middle) into position (low+1)
		//	int temp = window[low+1];
		//	window[low+1] = window[middle];
		//	window[middle] = temp;

		//	//Nibble from each end towards middle, swapping items when stuck
		//	ll = low + 1;
		//	hh = high;
		//	for (;;)
		//	{
		//		do{
		//			ll++;
		//		} while (window[low] > window[ll]);
		//		do
		//		{
		//			hh--;
		//		}
		//		while (window[hh] > window[low]);
		//		if (hh < ll)
		//			break;
		//		int temp = window[ll];
		//		window[ll] = window[hh];
		//		window[hh] = temp;
		//	}

		//	//Swap middle item (in position low) back into correct position
		//	temp = window[low];
		//	window[low] = window[hh];
		//	window[hh] = temp;
		//	//Re-set active partition
		//	if (hh <= median)
		//		low = ll;
		//	if (hh >= median)
		//		high = hh - 1;
		//}

		outputImage[y*width + x] = window[median];
}

int main(int argc, char* argv[])
{
	char* name;
	int bin_dim;
	int dim;
	if (argc < 3) {
		printf("The program requires you to specify the file, bin dimension and filter radius \n\n");
		return 1;
	}
	else
	{
		name = argv[1];
		try{
			bin_dim = std::stoi(argv[2]);
		}
		catch (std::exception const & e)
		{
			printf("--- ERROR : '%s' is not a valid integer\n", argv[1]);
			return 1;
		}
		try{
			dim = std::stoi(argv[3]);
		}
		catch (std::exception const & e)
		{
			printf("--- ERROR :'%s' is not a valid integer\n", argv[2]);
			return 1;
		}
	}

	int radius = (dim - 1) / 2;
	double t1, t2;
	Grid grid_p = Grid(bin_dim, bin_dim);
	Grid cuda = grid_p;

	printf("\n================== CONFIGURATION ==================\n\n");
	printf("--- Bin dimensions \t\t: %d x %d\n", bin_dim, bin_dim);
	printf("--- Filter dimensions \t\t: %d x %d\n", dim, dim);
	printf("--- OpenMP threads \t\t: %d\n", omp_get_max_threads());

	printf("\n=================== BINNING ===================\n\n");
	printf("--- Serial Binning ---\n\n");
	int success = ReadWrite::LoadData_s(grid_p, name);
	printf("--- OMP Binning ---\n\n");
	success = ReadWrite::LoadData_omp(grid_p, name);
	if (success == -1)
		return 2;
	Grid input = grid_p;
	printf("\n================= Serial Filtering =================\n\n");
	
	t1 = omp_get_wtime();
	Grid padded = grid_p.Pad(radius);
	t2 = omp_get_wtime();
	printf("--- Time for padding \t\t: %-12.5f sec.\n", t2 - t1);

	//three slower sorting algorithms
	/*t1 = omp_get_wtime();
	Filter::m_Filter_half(padded, grid_p, dim);
	t2 = omp_get_wtime();
	printf("--- Time for half bubble sort \t: %-12.5f sec.\n", t2 - t1);
	printf("--- Count after filter \t\t: %d\n", grid_p.Count());

	t1 = omp_get_wtime();
	Filter::m_Filter_quickselect(padded, grid_p, dim);
	t2 = omp_get_wtime();
	printf("--- Time for Quick Select \t: %-12.5f sec.\n", t2 - t1);
	printf("--- Count after filter \t\t: %d\n", grid_p.Count());*/
	
	t1 = omp_get_wtime();
	//Filter::m_Filter_quickselect2(padded, grid_p, dim);
	t2 = omp_get_wtime();
	printf("--- Time for Quick Select \t: %-12.5f sec.\n", t2 - t1);
	printf("--- Count after filter \t\t: %d\n", grid_p.Count());

	//padded.Print();
	printf("\n================== Cuda Filtering ==================\n\n");
	cudaError_t cudaStatus;

	/*-------------------
	*	The two commented out kernels use constant values for the grid
	*	and filter size
	---------------------*/

	/*cudaStatus = median_filter(cuda, input, R, 1);
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "addWithCuda failed!");
	return 1;
	}
	printf("--- Count after filter \t\t: %d\n", cuda.Count());
	
	cudaStatus = median_filter(cuda, input, R, 2);
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "addWithCuda failed!");
	return 1;
	}
	printf("--- Count after filter \t\t: %d\n", cuda.Count());*/
	cudaStatus = median_filter(cuda, padded, radius, 3);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addWithCuda failed!");
		return 1;
	}
	printf("--- Count after filter \t\t: %d\n", cuda.Count());

	printf("\n===================== Writing csv =====================\n\n");
	printf("--- Writing unfiltered grid\n");
	//ReadWrite::WriteData(input, "unfiltered.csv");
	//printf("--- Writing serial grid\n");
	//ReadWrite::WriteData(grid_p, "serial.csv");
	printf("--- Writing cuda grid\n");
	//ReadWrite::WriteData(cuda, "cuda.csv");

	printf("\n======================= Grids ======================\n\n");
	if (bin_dim < 20)
	{
		printf("--- Initial Grid ---\n");
		input.Print();
		printf("---- Serial Grid ---\n");
		grid_p.Print();
		printf("----- Cuda Grid ----\n");
		cuda.Print();
	}
	else
	{
		printf("--- Initial Subset ---\n");
		input.PrintRange(0, 20, 0, 20);
		printf("---- Serial Subset ---\n");
		grid_p.PrintRange(0, 20, 0, 20);
		printf("----- Cuda Subset ----\n");
		cuda.PrintRange(0, 20, 0, 20);
	}


    //cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}

// Helper function for using CUDA to apply a median filter
cudaError_t median_filter(Grid &output, const Grid &input, const int filter_radius, int padding_type)
{

    int *dev_i = 0;
    int *dev_o = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for input array.
	cudaStatus = cudaMalloc((void**)&dev_i, input.rows*input.cols * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	
	cudaStatus = cudaMalloc((void**)&dev_o, output.rows*output.cols * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
    // Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_i, input.grid, input.rows*input.cols * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	//cudaStatus = cudaBindTexture(NULL, tex, dev_i, input.rows*input.cols*sizeof(int));
	
	cudaThreadSetLimit(cudaLimitMallocHeapSize, (BLOCK_DIM + 2 * filter_radius)*(BLOCK_DIM + 2 * filter_radius) * (filter_radius * 2 + 1)*(filter_radius * 2 + 1) * sizeof(int));

	dim3 grid;
	dim3 threads;
	dim3 gr, th;
	if (padding_type == 1)
	{
//		grid = dim3((N + TILE_W - 1) / TILE_W, (N + TILE_H - 1) / TILE_W);
//		threads = dim3(BLOCK_H, BLOCK_W);
	}
	else
	{
		threads = dim3(BLOCK_DIM, BLOCK_DIM);
		grid = dim3((output.cols + threads.x - 1) / threads.x, (output.rows + threads.y - 1) / threads.y);
		th = dim3(T, 1);
		gr = dim3((output.cols + th.x - 1) / th.x, output.rows);
	}
	dim3 grid_warmup(1, 1);

	if (padding_type == 1)
	{
		//		median_Kernel_threads_padding <<<grid_warmup, threads>>>(dev_i, dev_o, output.cols, output.rows);
	}
	else if (padding_type == 2)
	{
		//	median_Kernel_border_padding <<<grid_warmup, threads>>>(dev_i, dev_o, output.cols, output.rows);
	}
	else
	{ 
		//median_Kernel_border_prepadded << <grid_warmup, threads, (threads.x +2*filter_radius)*(threads.y+2*filter_radius)*sizeof(int) >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		//median_Kernel_texture << <grid_warmup, th, ((T + 2 * filter_radius)*(1 + 2 * filter_radius) + (T*(filter_radius * 2 + 1)*(filter_radius * 2 + 1)))*sizeof(int) >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
	}

	cudaEventRecord(start);

    // launch the different kernels
	if (padding_type == 1)
	{
		//median_Kernel_threads_padding << <grid_warmup, threads >> >(dev_i, dev_o, output.cols, output.rows);
	}
	else if (padding_type == 2)
	{
	//	median_Kernel_border_padding << <grid, threads >> >(dev_i, dev_o, output.cols, output.rows);
	}
	else
	{
		//median_Kernel_border_prepadded << <grid, threads/*, ((threads.x + 2 * filter_radius)*(threads.y + 2 * filter_radius) + threads.x*threads.y*(filter_radius * 2 + 1)*(filter_radius * 2 + 1))*sizeof(int)*/ >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		if (filter_radius == 10)
		{
			median_Kernel_21 << <grid_warmup, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);

			median_Kernel_21 << <gr, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		}
		else if (filter_radius == 9)
			median_Kernel_19 << <gr, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		else if (filter_radius == 8)
			median_Kernel_17 << <gr, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		else if (filter_radius == 7)
			median_Kernel_15 << <gr, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		else if (filter_radius == 6)
			median_Kernel_13 << <gr, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		else if (filter_radius == 5)
			median_Kernel_11 << <gr, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		else if (filter_radius == 4)
			median_Kernel_9 << <gr, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		else if (filter_radius == 3)
			median_Kernel_7 << <gr, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		else if (filter_radius == 2)
			median_Kernel_5 << <gr, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
		else if (filter_radius == 1)
			median_Kernel_3 << <gr, th >> >(dev_i, dev_o, output.cols, output.rows, filter_radius);
	
	}

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	printf("--- Cuda Quick Select \t: %-12.5f sec\n", milliseconds/1000);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

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
	cudaStatus = cudaMemcpy(output.grid, dev_o, output.rows*output.cols * sizeof(int), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	
Error:
//1	cudaUnbindTexture(tex);
    cudaFree(dev_i);
    cudaFree(dev_o);

    return cudaStatus;
}

//the two kernels below use constant values
//__global__ void median_Kernel_border_padding(int *inputImage, int *outputImage, int width, int height)//int *inputImage, int *outputImage, unsigned int width, unsigned int height)
//{
//	__shared__ int smem[(BLOCK_DIM + 2 * R)][(BLOCK_DIM + 2 * R)];
//	int x_s = threadIdx.x + R;
//	int y_s = threadIdx.y + R;
//
//	int x = blockIdx.x*blockDim.x + threadIdx.x;
//	int y = blockIdx.y*blockDim.y + threadIdx.y;
//	if (x < N && y < N)
//	{
//		smem[y_s][x_s] = inputImage[y*(N)+x];
//
//		if (threadIdx.x == 0)
//		{
//			if (threadIdx.y == 0)
//			{
//				for (int j = -1; j >= -R; --j)
//					for (int i = -1; i >= -R; --i)
//					{
//						int offset_x = max(0, x + i);
//						int offset_y = max(0, y + j);
//						smem[y_s + j][x_s + i] = inputImage[offset_y*(N)+offset_x];
//					}
//
//			}
//			if (threadIdx.y == BLOCK_DIM - 1 || y == N - 1)
//			{
//				for (int j = 1; j <= R; ++j)
//					for (int i = -1; i >= -R; --i)
//					{
//						int offset_x = max(0, x + i);
//						int offset_y = min(y + j, height - 1);
//						smem[y_s + j][x_s + i] = inputImage[offset_y*(N)+offset_x];
//					}
//			}
//			for (int i = -1; i >= -R; --i)
//			{
//				int offset_x = max(0, x + i);
//				smem[y_s][x_s + i] = inputImage[y*(N)+offset_x];
//			}
//		}
//		if (threadIdx.x == BLOCK_DIM - 1 || x == N - 1)
//		{
//			if (threadIdx.y == 0)
//			{
//				for (int j = -1; j >= -R; --j)
//					for (int i = 1; i <= R; ++i)
//					{
//						int offset_x = min(x + i, width - 1);
//						int offset_y = max(0, y + j);
//						smem[y_s + j][x_s + i] = inputImage[offset_y*(N)+offset_x];
//					}
//			}
//			if (threadIdx.y == BLOCK_DIM - 1 || y == N - 1)
//			{
//				for (int j = 1; j <= R; ++j)
//					for (int i = 1; i <= R; ++i)
//					{
//						int offset_x = min(x + i, width - 1);
//						int offset_y = min(y + j, height - 1);
//						smem[y_s + j][x_s + i] = inputImage[offset_y*(N)+offset_x];
//					}
//			}
//			for (int i = 1; i <= R; ++i)
//			{
//				int offset_x = min(x + i, width - 1);
//				smem[y_s][x_s + i] = inputImage[y*(N)+offset_x];
//			}
//		}
//		if (threadIdx.y == 0)
//		{
//			for (int i = -1; i >= -R; --i)
//			{
//				int offset_y = max(0, y + i);
//				smem[y_s + i][x_s] = inputImage[offset_y*(N)+x];
//			}
//		}
//		if (threadIdx.y == BLOCK_DIM - 1 || y == N - 1)
//		{
//			for (int i = 1; i <= R; ++i)
//			{
//				int offset_y = min(y + i, height - 1);
//				smem[y_s + i][x_s] = inputImage[offset_y*(N)+x];
//			}
//		}
//
//		__syncthreads();
//		const int dim = D;
//		const int k = dim*dim / 2;
//		int a = 0;
//		unsigned int window[dim*dim];
//
//
//		for (int h_offset = -R; h_offset <= R; h_offset++)
//			for (int w_offset = -R; w_offset <= R; w_offset++)
//				window[a++] = smem[y_s + h_offset][x_s + w_offset];
//
//		int low, high;
//		int median;
//		int middle, ll, hh;
//
//		low = 0; high = dim*dim - 1; median = (low + high) / 2;
//		for (;;)
//		{
//			if (high <= low)
//				break;
//			if (high == low + 1)
//			{
//				if (window[low] > window[high])
//				{
//					int temp = window[low];
//					window[low] = window[high];
//					window[high] = temp;
//				}
//				break;
//			}
//
//			middle = (low + high) / 2;
//			if (window[middle] > window[high])
//			{
//				int temp = window[middle];
//				window[middle] = window[high];
//				window[high] = temp;
//			}
//			if (window[low] > window[high])
//			{
//				int temp = window[low];
//				window[low] = window[high];
//				window[high] = temp;
//			}
//			if (window[middle] > window[low])
//			{
//				int temp = window[low];
//				window[low] = window[middle];
//				window[middle] = temp;
//			}
//			int temp = window[low + 1];
//			window[low + 1] = window[middle];
//			window[middle] = temp;
//
//			ll = low + 1;
//			hh = high;
//			for (;;)
//			{
//				do{
//					ll++;
//				} while (window[low] > window[ll]);
//				do
//				{
//					hh--;
//				} while (window[hh] > window[low]);
//				if (hh < ll)
//					break;
//				int temp = window[ll];
//				window[ll] = window[hh];
//				window[hh] = temp;
//			}
//
//			temp = window[low];
//			window[low] = window[hh];
//			window[hh] = temp;
//
//			if (hh <= median)
//				low = ll;
//			if (hh >= median)
//				high = hh - 1;
//		}
//		outputImage[y*width + x] = window[median];
//
//	}
//}
//
//__global__ void median_Kernel_threads_padding(int *inputImage, int *outputImage, int width, int height)//int *inputImage, int *outputImage, unsigned int width, unsigned int height)
//{
//
//	__shared__ int smem[BLOCK_W * BLOCK_H];
//
//	int x = blockIdx.x*TILE_W + threadIdx.x - R;
//
//	int y = blockIdx.y*TILE_H + threadIdx.y - R;
//
//	int a = x;
//	int b = y;
//
//	int c = blockIdx.x*blockDim.x + threadIdx.x;
//	int d = blockIdx.y*blockDim.y + threadIdx.y;
//
//	{
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
//		smem[bindex] = inputImage[index];
//
//		__syncthreads();
//		if ((threadIdx.x >= R) && (threadIdx.x < BLOCK_W - R) && (threadIdx.y >= R) && (threadIdx.y < BLOCK_W - R) && (a < N && b < N))
//
//		{
//
//			const int dim = D;
//			const int k = dim*dim / 2;
//			int a = 0;
//			unsigned int window[dim*dim];
//
//
//			for (int h_offset = -R; h_offset <= R; h_offset++)
//				for (int w_offset = -R; w_offset <= R; w_offset++)
//					window[a++] = smem[bindex + h_offset*BLOCK_W + w_offset];
//			int low, high;
//			int median;
//			int middle, ll, hh;
//
//			low = 0; high = dim*dim - 1; median = (low + high) / 2;
//			for (;;)
//			{
//				if (high <= low)
//					break;
//				if (high == low + 1)
//				{
//					if (window[low] > window[high])
//					{
//						int temp = window[low];
//						window[low] = window[high];
//						window[high] = temp;
//					}
//					break;
//				}
//
//				middle = (low + high) / 2;
//				if (window[middle] > window[high])
//				{
//					int temp = window[middle];
//					window[middle] = window[high];
//					window[high] = temp;
//				}
//				if (window[low] > window[high])
//				{
//					int temp = window[low];
//					window[low] = window[high];
//					window[high] = temp;
//				}
//				if (window[middle] > window[low])
//				{
//					int temp = window[low];
//					window[low] = window[middle];
//					window[middle] = temp;
//				}
//				int temp = window[low + 1];
//				window[low + 1] = window[middle];
//				window[middle] = temp;
//
//				ll = low + 1;
//				hh = high;
//				for (;;)
//				{
//					do{
//						ll++;
//					} while (window[low] > window[ll]);
//					do
//					{
//						hh--;
//					} while (window[hh] > window[low]);
//					if (hh < ll)
//						break;
//					int temp = window[ll];
//					window[ll] = window[hh];
//					window[hh] = temp;
//				}
//
//				temp = window[low];
//				window[low] = window[hh];
//				window[hh] = temp;
//
//				if (hh <= median)
//					low = ll;
//				if (hh >= median)
//					high = hh - 1;
//			}
//
//			outputImage[index] = window[median];
//
//		}
//
//	}
//}

//these were early attemps, windows timeout cause a lot of self doubt

//__global__ void addKernel(int *inputImage, int *outputImage, int width, int height)//int *inputImage, int *outputImage, unsigned int width, unsigned int height)
//{
//	/*int from = 0, to = dim*dim - 1;
//	while (from < to)
//	{
//	int r = from, w = to;
//	int mid = window[(r + w) / 2];
//
//	while (r < w)
//	{
//	if (window[r] >= mid)
//	{
//	int temp = window[w];
//	window[w] = window[r];
//	window[r] = temp;
//	w--;
//	}
//	else
//	{
//	r++;
//	}
//	}
//	if (window[r] > mid)
//	{
//	r--;
//	}
//	if (k <= r)
//	{
//	to = r;
//	}
//	else
//	{
//	from = r + 1;
//	}
//	}*/
//
//	//int idc, val, min, max, inf, equal, sup, mxinf, minsup, estim;
//
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
//
//	//	__shared__ int smem[16 * 16];
//	//
//	//	int x = blockIdx.x*TILE_W + threadIdx.x - R;
//	//
//	//	int y = blockIdx.y*TILE_H + threadIdx.y - R;
//	//
//	//	//clamp to edge of image
//	////	if (blockIdx.x*blockDim.x + threadIdx.x < N + R && blockIdx.y*blockDim.y + threadIdx.y < N + R)
//	//	//{
//	//		x = max(0, x);
//	//
//	//		x = min(x, width - 1);
//	//
//	//		y = max(y, 0);
//	//
//	//		y = min(y, height - 1);
//	//
//	//
//	//
//	//		unsigned int index = y*width + x;
//	//
//	//		unsigned int bindex = threadIdx.y * blockDim.y + threadIdx.x;
//	//
//	//
//	//
//	//		//each thread copies its pixel of the block to shared memory
//	//
//	//		smem[bindex] = inputImage[index];
//	//
//	//		__syncthreads();
//	//
//	//		if ((threadIdx.x >= R) && (threadIdx.x < (BLOCK_W - R)) &&
//	//			(threadIdx.y >= R) && (threadIdx.y < (BLOCK_H - R)))
//	//		{
//	//			
//	//			min = max = smem[ib * blockDim.x + jb];
//	//
//	//				for (idc = 0; idc < 2 * R + 1; ++idc)
//	//				{
//	//					val = smem[__mul24(ib + idc, blockDim.x) + jb];
//	//				if (val < min) min = val;
//	//				if (val > max) max = val;
//	//				}
//	//
//	//				while (1)
//	//				{
//	//				estim = (min + max) / 2;
//	//				inf = sup = equal = 0;
//	//				mxinf = min;
//	//				minsup = max;
//	//				for (idc = 0; idc < 2 * R + 1; ++idc)
//	//				{
//	//					val = smem[__mul24(ib + idc, blockDim.x) + jb];
//	//				if (val < estim)
//	//				{
//	//				inf++;
//	//				if (val > mxinf) mxinf = val;
//	//				}
//	//				else if (val > estim)
//	//				{
//	//				sup++;
//	//				if (val < minsup) minsup = val;
//	//				}
//	//				else equal++;
//	//				}
//	//				if ((inf <= (R + 1)) && (sup <= (R + 1))) break;
//	//				else if (inf > sup) max = mxinf;
//	//				else min = minsup;
//	//				}
//	//				if (inf >= R + 1) val = mxinf;
//	//				else if (inf + equal >= R + 1) val = estim;
//	//				else val = minsup;
//	//
//	//			/*int      val, i, less, greater, equal,  min, max, guess, maxltguess, mingtguess;
//	//
//	//			min = max = smem[bindex];
//	//			for (i =-1; i<=1; i++) {
//	//			val = smem[__mul24(threadIdx.y + i, blockDim.x) + threadIdx.x];
//	//			if (val<min) min = val;
//	//			if (val> max) max = val;
//	//			}
//	//
//	//			while (1) {
//	//				guess = (min + max) / 2;
//	//				less = 0; greater = 0; equal = 0;
//	//				maxltguess = min;
//	//				mingtguess = max;
//	//				for (i = -1; i<=1; i++) {
//	//					val = smem[__mul24(threadIdx.y + i, blockDim.x) + threadIdx.x];
//	//					if (val<guess) {
//	//						less++;
//	//						if (val>maxltguess) maxltguess = val;
//	//					}
//	//					else if (val>guess) {
//	//						greater++;
//	//					if (val<mingtguess) mingtguess = val;
//	//					}
//	//					else equal++;
//	//				}
//	//				if (less <= (R + 1) && greater <= (R + 1)) break;
//	//				else if (less>greater) max = maxltguess;
//	//				else min = mingtguess;
//	//			}
//	//			if (less >= (R + 1)) val = maxltguess;
//	//			else if (less + equal >= (R + 1)) val = guess;
//	//			else val = mingtguess;*/
//	//			/*printf("dim[%d][%d] | %-10d : %-10d, %-10d, %-10d, %-10d, %-10d, %-10d, %-10d, %-10d, %-10d \n", y, x, inputImage[index],
//	//				smem[bindex - blockDim.x - 1], smem[bindex - blockDim.x], smem[bindex - blockDim.x + 1],
//	//				smem[bindex - 1], smem[bindex], smem[bindex + 1],
//	//				smem[bindex + blockDim.x - 1], smem[bindex + blockDim.x], smem[bindex + blockDim.x + 1]);*/
//	//			//outputImage[index] =  smem[bindex];
//	//			outputImage[index] = val;
//	//
//	//		}
//	//}
//
//	//extern __shared__ int data[];
//	//__shared__ int data[(R * 2 + 1) * (T + 2 * R)];
//	//int col = blockIdx.x * blockDim.x + threadIdx.x;
//	//int row = blockIdx.y * blockDim.y + threadIdx.y;
//	//int d_col = col + R;
//	//int d_row = blockIdx.y;
//
//	//if (row < N && col < N)
//	//{
//
//	//	for (int r = -R; r <= R; ++r)
//	//	{
//	//		//col += r;
//	//		int i_row = row + r;
//	//		//col = max(0, col);
//
//	//		//col = min(col, N - 1);
//
//	//		i_row = max(i_row, 0);
//
//	//		i_row = min(i_row, N - 1);
//	//		//data[1028 + col] = 0; //inputImage[row*N + col];
//	//		data[(R + r) * (T + 2 * R) + d_col] = inputImage[i_row*N + col];
//	//		if (col == 0)
//	//		{
//	//			for (int i = 0; i <= R;++i)
//	//				data[(R + r) * (T + 2 * R) + col+i] = inputImage[i_row*N + col];
//	//		}
//	//		if (col == N - 1 || col == blockDim.x)
//	//		{
//	//			for (int i = 0; i <= R; ++i)
//	//				data[(R + r) * (T + 2 * R) + d_col + i] = inputImage[i_row*N + col];
//	//		}
//	//	}
//	//	__syncthreads();
//	/*printf("dim : %d | [%d][%d] | %-10d : %-9d, %-9d, %-9d, %-9d, %-9d, %-9d, %-9d, %-9d, %-9d \n", blockDim.y, d_row, col, inputImage[(row) * N + col],
//	data[(0) * 514 + d_col - 1], data[(0) * 514 + d_col], data[(0) * 514 + d_col + 1],
//	data[(1) * 514 + d_col - 1], data[(1) * 514 + d_col], data[(1) * 514 + d_col + 1],
//	data[(2) * 514 + d_col - 1], data[(2) * 514 + d_col], data[(2) * 514 + d_col + 1]);*/
//	//printf("dim : %d | [%d][%d] | %-10d : %-9d, %-9d, %-9d\n", blockDim.y, d_row, N, inputImage[(row) * N + col], data[(0) * 514 + N], data[(1)* 514 + N], data[(2) * 514 + N]);
//	//int k = 0;
//	//int *window = new int[D*D];
//	//for (int dy = -R; dy <= R; ++dy)
//	//{
//	//	for (int dx = -R; dx <= R; ++dx)
//	//	{
//	//		window[k++] = data[(R + dy) * (T + 2 * R) + d_col + dx];
//	//	}
//	//}
//	//for (int j = 0; j < (D*D+1)/2; ++j)
//	//{
//	//	//   find position of minimum element
//	//	int min = j;
//	//	for (int l = j + 1; l < (D*D); ++l)
//	//		if (window[l] < window[min])
//	//			min = l;
//	//	const int temp = window[j];
//	//	window[j] = window[min];
//	//	window[min] = temp;
//	//}
//	/*int      val, i, less, greater, equal,  min, max, guess, maxltguess, mingtguess;
//
//	min = max = data[__mul24(threadIdx.y + R, blockDim.x + 2 * R) + threadIdx.x + R];
//	for (i =0; i<D; i++) {
//	val = data[__mul24(threadIdx.y + i, blockDim.x) + threadIdx.x];
//	if (val<min) min = val;
//	if (val> max) max = val;
//	}
//	while (1) {
//	guess = (min + max) / 2;
//	less = 0; greater = 0; equal = 0;
//	maxltguess = min;
//	mingtguess = max;
//	for (i = 0; i<D; i++) {
//	val = data[__mul24(threadIdx.y + i, blockDim.x) + threadIdx.x];
//	if (val<guess) {
//	less++;
//	if (val>maxltguess) maxltguess = val;
//	}
//	else if (val>guess) {
//	greater++;
//	if (val<mingtguess) mingtguess = val;
//	}
//	else equal++;
//	}
//	if (less <= (R + 1) && greater <= (R + 1)) break;
//	else if (less>greater) max = maxltguess;
//	else min = mingtguess;
//	}
//	if (less >= (R + 1)) val = maxltguess;
//	else if (less + equal >= (R + 1)) val = guess;
//	else val = mingtguess;*/
//
//	//outputImage[(row)* N + col] = val;
//	//outputImage[(row)* N + col] = data[(R)* (T + 2 * R) + d_col/*(1) * 514 + d_col*/];
//	//outputImage[(row)* N + col] = 1;// window[(D*D) / 2];
//	//}
//}

//__global__ void median_Kernel_1d(int* inputImage, int* outputImage, const int width, const int height, const int radius)
//{
//	const int dim = radius * 2 + 1;
//	extern __shared__ int shared[];
//	int *smem = &shared[0];
//	int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
//	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
//	int x = blockIdx.x*blockDim.x + threadIdx.x;
//	int y = blockIdx.y;
//	int x_o = threadIdx.x + radius;
//	int y_o = radius;
//	if (x >= width)
//		return;
//	//for (int j = -radius; j <= radius; ++j)
//	//for (int i = -radius; i <= radius; ++i)
//	//if (smem[(y_o + j)*(blockDim.x + 2 * radius) + x_o + i] == 0)
//	//smem[(y_o + j)*(blockDim.x + 2 * radius) + x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_o + i];
//	for (int j = -radius; j <= radius; ++j)
//		if (smem[(y_o + j)*(blockDim.x + 2 * radius) + x_o] == 0)
//			smem[(y_o + j)*(blockDim.x + 2 * radius) + x_o] = inputImage[(y + radius + j)*(width + 2 * radius) + x_o];
//	if (threadIdx.x == 0)
//		for (int j = -radius; j <= radius; ++j)
//			for (int i = -radius; i < 0; ++i)
//				if (smem[(y_o + j)*(blockDim.x + 2 * radius) + x_o + i] == 0)
//					smem[(y_o + j)*(blockDim.x + 2 * radius) + x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_o + i];
//	if (threadIdx.x == blockDim.x - 1 || threadIdx.x == width - 1)
//		for (int j = -radius; j <= radius; ++j)
//			for (int i = 1; i <= radius; ++i)
//				if (smem[(y_o + j)*(blockDim.x + 2 * radius) + x_o + i] == 0)
//					smem[(y_o + j)*(blockDim.x + 2 * radius) + x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_o + i];
//
//	__syncthreads();
//	//find the offset of the current window
//	int a = (threadIdx.x)*dim*dim;
//	//populate the window with the values around the current index
//	for (int h_offset = -radius; h_offset <= radius; h_offset++)
//		for (int w_offset = -radius; w_offset <= radius; w_offset++)
//			window[a++] = smem[(y_o + h_offset) * (blockDim.x + 2 * radius) + x_o + w_offset];
//
//	int low, high;
//	int median;
//	int middle, ll, hh;
//	/*
//	* This Quickselect routine is based on the algorithm described in
//	* "Numerical recipes in C", Second Edition,
//	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
//	* This code by Nicolas Devillard - 1998. Public domain.
//	*/
//	low = (threadIdx.x)*dim*dim; high = (threadIdx.x)*dim*dim + (dim*dim - 1); median = (low + high) / 2;
//	for (;;)
//	{
//		if (high <= low)
//			break;
//		if (high == low + 1)
//		{
//			if (window[low] > window[high])
//			{
//				int temp = window[low];
//				window[low] = window[high];
//				window[high] = temp;
//			}
//			break;
//		}
//		// Find median of low, middle and high items; swap into position low
//		middle = (low + high) / 2;
//		if (window[middle] > window[high])
//		{
//			int temp = window[middle];
//			window[middle] = window[high];
//			window[high] = temp;
//		}
//		if (window[low] > window[high])
//		{
//			int temp = window[low];
//			window[low] = window[high];
//			window[high] = temp;
//		}
//		if (window[middle] > window[low])
//		{
//			int temp = window[low];
//			window[low] = window[middle];
//			window[middle] = temp;
//		}
//		//Swap low item (now in position middle) into position (low+1)
//		int temp = window[low + 1];
//		window[low + 1] = window[middle];
//		window[middle] = temp;
//
//		//Nibble from each end towards middle, swapping items when stuck
//		ll = low + 1;
//		hh = high;
//		for (;;)
//		{
//			do{
//				ll++;
//			} while (window[low] > window[ll]);
//			do
//			{
//				hh--;
//			} while (window[hh] > window[low]);
//			if (hh < ll)
//				break;
//			int temp = window[ll];
//			window[ll] = window[hh];
//			window[hh] = temp;
//		}
//
//		//Swap middle item (in position low) back into correct position
//		temp = window[low];
//		window[low] = window[hh];
//		window[hh] = temp;
//		//Re-set active partition
//		if (hh <= median)
//			low = ll;
//		if (hh >= median)
//			high = hh - 1;
//	}
//
//	outputImage[y*width + x] = window[median];
//}


//__global__ void median_Kernel_torben(int* inputImage, int* outputImage, const int width, const int height, const int radius)
//{
//	const int dim = radius * 2 + 1;
//	__shared__ int smem[(1 + 2 * R)][(T + 2 * R)];
//	//int *smem = &shared[0];
//	//int *window = &shared[(blockDim.x+2*radius)*(1+2*radius)];
//	//int *window = &shared[(blockDim.x + 2 * radius)*(1 + 2 * radius)];
//	int x = blockIdx.x*blockDim.x + threadIdx.x;
//	int y = blockIdx.y;
//	int x_p = x + radius;
//	int x_o = threadIdx.x + radius;
//	int y_o = radius;
//	if (x >= width)
//		return;
//
//	if (threadIdx.x == 0)
//		for (int j = 0; j < (1 + 2 * R); ++j)
//			for (int i = 0; i < (T + 2 * R); ++i)
//				smem[j][i] = 0;
//
//	__syncthreads();
//
//	for (int j = -radius; j <= radius; ++j)
//		for (int i = -radius; i <= radius; ++i)
//			if (smem[y_o + j][x_o + i] == 0)
//				smem[y_o + j][x_o + i] = inputImage[(y + radius + j)*(width + 2 * radius) + x_p + i];
//
//	__syncthreads();
//	/*
//
//	* The following code is public domain.
//
//	* Algorithm by Torben Mogensen, implementation by N. Devillard.
//
//	* This code in public domain.
//
//	*/
//
//
//	int i, less, greater, equal;
//	int min, max, guess, maxltguess, mingtguess;
//	min = max = smem[0][threadIdx.x];
//	for (i = 1; i < dim*dim; i++) {
//		int row = i / dim;
//		int col = i%dim;
//		if (smem[row][threadIdx.x + col] < min) min = smem[row][threadIdx.x + col];
//		if (smem[row][threadIdx.x + col] > max) max = smem[row][threadIdx.x + col];
//	}
//	while (1) {
//		guess = (min + max) / 2;
//		less = 0; greater = 0; equal = 0;
//		maxltguess = min;
//		mingtguess = max;
//		for (i = 0; i<dim*dim; i++)
//		{
//			int row = i / dim;
//			int col = i%dim;
//			if (smem[row][threadIdx.x + col]<guess)
//			{
//				less++;
//				if (smem[row][threadIdx.x + col]>maxltguess) maxltguess = smem[row][threadIdx.x + col];
//			}
//			else if (smem[row][threadIdx.x + col]>guess)
//			{
//				greater++;
//
//				if (smem[row][threadIdx.x + col] <mingtguess) mingtguess = smem[row][threadIdx.x + col];
//
//			}
//			else equal++;
//		}
//		if (less <= (dim*dim + 1) / 2 && greater <= (dim*dim + 1) / 2) break;
//		else if (less>greater) max = maxltguess;
//		else min = mingtguess;
//	}
//	int val = 0;
//	if (less >= (dim*dim + 1) / 2) val = maxltguess;
//	else if (less + equal >= (dim*dim + 1) / 2) val = guess;
//	else val = mingtguess;
//
//	outputImage[y*width + x] = val;
//
//	
//}