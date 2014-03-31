#include "cuda_memory.h"
#include "stdio.h"
#include "cuda_kernels.h"
#include "model.h"

constExecutionParameters hConstParams;
varExecutionParameters hVarParams;
devMemory hDevMem;

//=============================================================================
void cuAllocateDataMemory()
{
	printf(" Cuda_module: cu_init_mem_()\n");
	
	if(hConstParams.Nx <= 0 || hConstParams.Ny <= 0)
		printf(" Cuda module: Warning: Nx (%d) or Ny (%d) invalid.\n", hConstParams.Nx, hConstParams.Ny);
	
	const size_t elemcount = hConstParams.Nx * hConstParams.Ny;
	const size_t sizeInts = elemcount * sizeof(fortInteger);
	const size_t sizeReals = elemcount * sizeof(fortReal8);
	
	// Allocate device memory
	checkCudaErrors( cudaMalloc((void **)&hDevMem.ocean_eta, 	sizeInts) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.ocean_u,   	sizeInts) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.ocean_v,   	sizeInts) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.SWM_Coef_eta, 	sizeReals * 9) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.SWM_Coef_u,   	sizeReals * 11) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.SWM_Coef_v,   	sizeReals * 11) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.impl_eta, 	sizeReals) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.impl_u,   	sizeReals) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.impl_v,   	sizeReals) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.G_eta, 		sizeReals * 2) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.G_u,   		sizeReals * 2) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.G_v,   		sizeReals * 2) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.SWM_eta, 		sizeReals * 2) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.SWM_u,   		sizeReals * 2) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.SWM_v,   		sizeReals * 2) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.diag_psi,		sizeReals) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.F_eta, 		sizeReals) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.F_x,   		sizeReals) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.F_y,   		sizeReals) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.H_u,   		sizeReals) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.H_v,   		sizeReals) );
	checkCudaErrors( cudaMalloc((void **)&hDevMem.cos_lat_v, 	sizeof(fortReal8) * hConstParams.Ny) );
}

void cuFreeDataMemory()
{
	cudaDeviceSynchronize(); // wait for other tasks to finish
	cudaFree(hDevMem.ocean_eta);
	cudaFree(hDevMem.ocean_u);
	cudaFree(hDevMem.ocean_v);
	cudaFree(hDevMem.SWM_Coef_eta);
	cudaFree(hDevMem.SWM_Coef_u);
	cudaFree(hDevMem.SWM_Coef_v);
	cudaFree(hDevMem.impl_eta);
	cudaFree(hDevMem.impl_u);
	cudaFree(hDevMem.impl_v);
	cudaFree(hDevMem.G_eta);
	cudaFree(hDevMem.G_u);
	cudaFree(hDevMem.G_v);
	cudaFree(hDevMem.SWM_eta);
	cudaFree(hDevMem.SWM_u);
	cudaFree(hDevMem.SWM_v);
	cudaFree(hDevMem.diag_psi);
	cudaFree(hDevMem.F_eta);
	cudaFree(hDevMem.F_x);
	cudaFree(hDevMem.F_y);
	cudaFree(hDevMem.H_u);
	cudaFree(hDevMem.H_v);
	cudaFree(hDevMem.cos_lat_v);
}


void __checkCudaErrors( cudaError err, const char *file, const int line )
{
	if( cudaSuccess != err) {
		fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",
			file, line, (int)err, cudaGetErrorString( err ) );
		exit(-1);
	}
}
