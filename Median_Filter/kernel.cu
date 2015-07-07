#define N 512
#define BLOCK_DIM 32

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
//#include "Filter.h"
#include "Grid.h"
#include "ReadWrite.h"
#include "Filter.h"


cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

__global__ void addKernel(int *c, const int *a, const int *b)
{
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int i = col + row * N;
	if (col < N && row < N) 
		c[i] = a[i] + b[i];
}



int main()
{
	char* name = "Points_[1.0e+08]_Noise_[030]_Normal.bin";
	double t1, t2;
	//Grid grid_s = Grid(4096, 4096);
	Grid grid_p = Grid(100, 100);
	//Grid grid_p_b = Grid(4096, 4096);
	printf("OpenMP threads: %d\n", omp_get_max_threads());
	//t1 = omp_get_wtime();
	printf("Binning file");
	ReadWrite::LoadData_omp(grid_p, name);
	//t2 = omp_get_wtime();
	//printf("Time for omp binning: %12.3f sec, checksum=%d (must be 100000000).\n", t2 - t1, grid_p.Count());
	//t1 = omp_get_wtime();
	//ReadWrite::LoadData_s(grid_s, name);
	//ReadWrite::LoadData_omp_buffer(grid_p_b, name);
	//t2 = omp_get_wtime();
	//printf("Time for serial binning: %12.3f sec, checksum=%d (must be 100000000).\n", t2 - t1, grid_s.Count());
	//grid_s.Print();
	//t1 = 0.0;
	//t2 = 0.0;	
	printf("Writing unfiltered");
	ReadWrite::WriteData(grid_p, "unfiltered.csv");
	
	
	
	//grid_p.Print();
	printf("filtering");
	t1 = omp_get_wtime();
	Filter::m_Filter_extended(grid_p, 3);
	t2 = omp_get_wtime();
	printf("Time for serial Filtering: %12.3f sec.\n", t2 - t1);
	printf("writing filtered\n");
	ReadWrite::WriteData(grid_p, "filtered.csv");
	//grid_p.Print();
	//grid_p_b.Print();
    //const int arraySize = 5;
    //int *a = new int[N*N];
	//int *b = new int[N*N];
	//int *c = new int[N*N];
	//int arraySize = N*N;

	//for (int x = 0; x < N*N; ++x)
	//{
	//	a[x] = 1;
	//	b[x] = 0;
	//	c[x] = 0;
	//}

    // Add vectors in parallel.
    //cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
    //if (cudaStatus != cudaSuccess) {
    //    fprintf(stderr, "addWithCuda failed!");
     //   return 1;
   // }
	

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    /*cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }*/

    return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
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
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

	dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
	dim3 dimGrid((int)ceil(N / (float)dimBlock.x), (int)ceil(N / (float)dimBlock.y));
    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<dimGrid, dimBlock>>>(dev_c, dev_a, dev_b);

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
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}
