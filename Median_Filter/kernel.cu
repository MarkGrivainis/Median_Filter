#ifndef __CUDACC__  
#define __CUDACC__
#endif

#define N 6
#define BLOCK_DIM 8

#define TILE_W 6
#define TILE_H 6
#define R 1
#define D (R*2+1)
#define BLOCK_W (TILE_W + (2*R))
#define BLOCK_H (TILE_H + (2*R))


#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>

#include "Grid.h"
#include "ReadWrite.h"
#include "Filter.h"


cudaError_t addWithCuda(int *outputImage, const int *inputImage, unsigned int paddedsize, unsigned int gridSize);

__global__ void addKernel(int *inputImage, int *outputImage, unsigned int width, unsigned int height)
{
	__shared__ int smem[BLOCK_W*BLOCK_H];

	int x = blockIdx.x*TILE_W + threadIdx.x;
	int y = blockIdx.y*TILE_H + threadIdx.y;
	//		if (x < 5 && y < 5)
	//		{
	//# if __CUDA_ARCH__>=200
	//	printf("%d \n", inputImage[x ]);
	//#endif
//}		
	//printf("x : %d, y : %d \n", x, y);
		x = max(0, x);
		x = min(x, width - 1);
		y = max(y, 0);
		y = min(y, height - 1);
		//printf("x : %d, y : %d \n", x, y);

		unsigned int index = y * width + x;
		unsigned int bindex = threadIdx.y*blockDim.y + threadIdx.x;
		//printf("current [%d][%d] : %d\n", threadIdx.y, threadIdx.x, inputImage[index]);
		smem[bindex] = inputImage[index];
		__syncthreads();

		if ((threadIdx.x >= R) && (threadIdx.x < (BLOCK_W - R)) &&
			(threadIdx.y >= R) && (threadIdx.y < (BLOCK_H - R)))
		{
			printf("current [%d][%d] : %d\n", threadIdx.y, threadIdx.x, smem[bindex]);
			float sum = 0;
			int k = 0;
			int *window = new int[9];
			for (int dy = -R; dy <= R; ++dy)
			{
				for (int dx = -R; dx <= R; ++dx)
				{
					window[k++] = smem[bindex + (dy*blockDim.x) + dx];
				}
			}
			for (int j = 0; j < 5; ++j)
			{
				//   Find position of minimum element
				int min = j;
				for (int l = j + 1; l < 9; ++l)
					if (window[l] < window[min])
						min = l;
				//   Put found minimum element in its place
				const int temp = window[j];
				window[j] = window[min];
				window[min] = temp;
			}

			//outputImage[index] = window[4];
			outputImage[index] = smem[bindex];

		}
}




int main()
{
	int size = 3;
	char* name = "Points_[1.0e+08]_Noise_[030]_Normal.bin";
	double t1, t2;
	//Grid grid_s = Grid(512, 512);
	Grid grid_p = Grid(6, 6);
	//Grid grid_p_b = Grid(4096, 4096);
	printf("OpenMP threads: %d\n", omp_get_max_threads());
	//t1 = omp_get_wtime();
	printf("Binning file\n");
	ReadWrite::LoadData_omp(grid_p, name);
	//grid_p.Print();
	printf("padding grid\n");
	int radius = (int)(size - 1) / 2;
	Grid padded = grid_p.Pad(radius);
	padded.Print();
	//t1 = omp_get_wtime();
	//Filter::m_Filter_fullsort(padded, grid_p, radius);
	//t2 = omp_get_wtime();
	//printf("Time for full Filtering: %12.3f sec.\n", t2 - t1);
	//printf("filtered grid\n");
	//grid_p.Print();
	t1 = omp_get_wtime();
	//Filter::m_Filter_half(padded, grid_p, size);
	t2 = omp_get_wtime();
	printf("Time for half Filtering: %12.3f sec.\n", t2 - t1);
	//grid_p.Print();
	//printf("grid[0][0] = %d\n", grid_p.grid[0]);
	//t2 = omp_get_wtime();
	//printf("Time for omp binning: %12.3f sec, checksum=%d (must be 100000000).\n", t2 - t1, grid_p.Count());
	//t1 = omp_get_wtime();
	//ReadWrite::LoadData_omp(grid_s, name);
	//ReadWrite::LoadData_s(grid_s, name);
	//ReadWrite::LoadData_omp_buffer(grid_p_b, name);
	//t2 = omp_get_wtime();
	//printf("Time for serial binning: %12.3f sec, checksum=%d (must be 100000000).\n", t2 - t1, grid_s.Count());
	//grid_s.Print();
	//t1 = 0.0;
	//t2 = 0.0;	
	//printf("Writing unfiltered");
	//ReadWrite::WriteData(grid_p, "unfiltered.csv");
	
	
	
	/*grid_p.Print();
	printf("grid_count : %d\n", grid_p.Count());
	printf("filtering");
	t1 = omp_get_wtime();
	Filter::m_Filter_extended(grid_p, radius);
	t2 = omp_get_wtime();
	printf("Time for serial Filtering: %12.3f sec.\n", t2 - t1);
	printf("writing filtered\n");*/
	//ReadWrite::WriteData(grid_p, "filtered.csv");
	//grid_p.Print();
	//grid_p_b.Print();
    //const int arraySize = 5;
    /*int *a = new int[N*N];
	int *b = new int[N*N];
	int *c = new int[N*N];
	unsigned int arraySize = N*N;

	for (int x = 0; x < N*N; ++x)
	{
		a[x] = 1;
		b[x] = 0;
		c[x] = 0;
	}*/

    // Add vectors in parallel.
    //cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
	printf("filtered grid\n");
	grid_p.Print();
	//padded.Print();
	cudaError_t cudaStatus = addWithCuda(grid_p.grid, padded.grid, padded.cols*padded.rows, grid_p.rows*grid_p.cols);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }
	grid_p.Print();
	//printf("c[0][0] = %d\n", c[0]);

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
cudaError_t addWithCuda(int *outputImage, const int *inputImage, unsigned int paddedsize, unsigned int gridSize)
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

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_o, gridSize * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_i, paddedsize * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_i, inputImage, paddedsize * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

	dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
	dim3 dimGrid((N + dimBlock.x)/dimBlock.x, (N + dimBlock.y)/dimBlock.y);
	cudaEventRecord(start);
    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<dimGrid, dimBlock>>>(dev_i, dev_o, 8, 8);
	cudaEventRecord(stop);

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
    cudaStatus = cudaMemcpy(outputImage, dev_o, gridSize * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	printf("cudatime = %12.12f\n", milliseconds/1000);
Error:
    cudaFree(dev_i);
    cudaFree(dev_o);
    
    return cudaStatus;
}
