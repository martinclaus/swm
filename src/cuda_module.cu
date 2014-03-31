#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cuda_module.h"
#include "cuda_memory.h"
#include "cuda_kernels.h"
#include "model.h"

//=============================================================================
enum module_state {
		not_initialized, 	// didn't try
		initialized,		// tried and succeeded
		not_available		// tried and failed
} state;

//=============================================================================
dim3 grid, blocksize;
cudaStream_t stream[3];

//=============================================================================
// Testing device capabilities
void cu_init_()
{
	// find devices
	int devCount = 0;
	cudaGetDeviceCount(&devCount);

	if(devCount == 0) {
		state = not_available;
		printf(" Cuda module not available: No device found.\n");
		return;
	}

	// use first device
	const int devId = 0;
	cudaDeviceProp deviceProp;

	cudaSetDevice(devId);
	cudaGetDeviceProperties(&deviceProp, devId);

	if(deviceProp.major <= 1 || (deviceProp.major <= 1 && deviceProp.minor < 3)) {
		state = not_available;
		printf(" Cuda Error: compute capability of GPU < 1.3 : No support for double precision.\n");
		return;
	}
	
	printf(" Using GPU Device %d: %s with compute capability %d.%d (%d Multiprocessors)\n", 
		devId, deviceProp.name, deviceProp.major, deviceProp.minor, deviceProp.multiProcessorCount);

	// Create execution streams
#ifdef CUDA_USE_STREAMS
	for(int i = 0; i < 3; i++) {
		cudaStreamCreate(&stream[i]);
	}
#endif
		
	state = initialized;	
}

//=============================================================================
// The grid is fitted to the blocksize. If the grid size is not a multiple of the blocksize, it is padded.
// This allows for use of the optimal block size and also odd grid dimensions. The extra elements of the grid
// are not computed by the kernel.
void fitGridSize() {
	// Set kernel size
	blocksize = dim3(8, 4);
	
	// Padding
	int kernelNx = hConstParams.Nx;
	if(hConstParams.Nx % blocksize.x != 0)
	  kernelNx = (hConstParams.Nx / blocksize.x + 1) * blocksize.x;
	
	int kernelNy = hConstParams.Ny;
	if(kernelNy % blocksize.y != 0)
	  kernelNy = (hConstParams.Ny / blocksize.y + 1) * blocksize.y;
	
	grid = dim3(kernelNx / blocksize.x, kernelNy / blocksize.y);
	
	// be verbose
	if(kernelNx == hConstParams.Nx && kernelNy == hConstParams.Ny)
		printf(" Cuda_module: No Padding necessary.\n");
	else
		printf(" Cuda_module: Padded to size %dx%d.\n", kernelNx, kernelNy);
}

//=============================================================================
void cu_setconstants_(fortInteger *pNx, fortInteger *pNy, fortReal8 *pDt, fortReal8 *pAB_C1, fortReal8 *pAB_C2, fortReal8 *A, fortReal8 *dTheta, fortReal8 *dLambda)
{	
	if(state != initialized)
		return;
	
	// Set initial values
	hVarParams.N0 = 0;
	hVarParams.N0p1 = 1;
	hVarParams.NG0 = 1;
	hVarParams.NG0m1 = 0;

	// Set constants
	hConstParams.Nx = *pNx;
	hConstParams.Ny = *pNy;
	hConstParams.dt = *pDt;
	hConstParams.AB_C1 = *pAB_C1;
	hConstParams.AB_C2 = *pAB_C2;
	hConstParams.A = *A;
	hConstParams.dTheta = *dTheta;
	hConstParams.dLambda = *dLambda;

	fitGridSize();
}

//=============================================================================
void cu_setfields_(
		bool *ocean_eta, bool *ocean_u, bool *ocean_v,
		fortReal8 *SWM_Coef_eta, fortReal8 *SWM_Coef_u, fortReal8 *SWM_Coef_v,
		fortReal8 *impl_eta, fortReal8 *impl_u, fortReal8 *impl_v,
		fortReal8 *G_eta, fortReal8 *G_u, fortReal8 *G_v,
		fortReal8 *SWM_eta, fortReal8 *SWM_u, fortReal8 *SWM_v,
		fortReal8 *F_eta, fortReal8 *F_x, fortReal8 *F_y,
		fortReal8 *H_u, fortReal8 *H_v,
		fortReal8 *cos_lat_v
		)
{	
	if(state != initialized)
		return;
	
	if(hConstParams.Nx <= 0 || hConstParams.Ny <= 0) {
	  printf(" Cuda module: Warning: Nx (%d) or Ny (%d) invalid.\n", hConstParams.Nx, hConstParams.Ny);
	  return;
	}
	
 	cuAllocateDataMemory();
	
	const size_t elemcount = hConstParams.Nx * hConstParams.Ny;
	const size_t sizeInts = elemcount * sizeof(fortInteger);
	const size_t sizeReals = elemcount * sizeof(fortReal8);
		
	// Convert ocean field from boolean to integers
	fortInteger *tmp_ocean = (fortInteger *) malloc(sizeInts);
	
	for(int i = 0; i < elemcount; i++)
	  tmp_ocean[i] = ocean_eta[i];
	checkCudaErrors( cudaMemcpy(hDevMem.ocean_eta, tmp_ocean, sizeInts, cudaMemcpyHostToDevice) );
	
	for(int i = 0; i < elemcount; i++)
	  tmp_ocean[i] = ocean_u[i];
	checkCudaErrors( cudaMemcpy(hDevMem.ocean_u,   tmp_ocean, sizeInts, cudaMemcpyHostToDevice) );
	
	for(int i = 0; i < elemcount; i++)
	  tmp_ocean[i] = ocean_v[i];
	checkCudaErrors( cudaMemcpy(hDevMem.ocean_v,   tmp_ocean, sizeInts, cudaMemcpyHostToDevice) );
	
	free(tmp_ocean);
		
	// Copy the data to the device
	checkCudaErrors( cudaMemcpy(hDevMem.SWM_Coef_eta, SWM_Coef_eta, sizeReals * 9,  cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.SWM_Coef_u,   SWM_Coef_u,   sizeReals * 11, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.SWM_Coef_v,   SWM_Coef_v,   sizeReals * 11, cudaMemcpyHostToDevice) );
	
	checkCudaErrors( cudaMemcpy(hDevMem.impl_eta, impl_eta, sizeReals, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.impl_u,   impl_u,   sizeReals, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.impl_v,   impl_v,   sizeReals, cudaMemcpyHostToDevice) );

	checkCudaErrors( cudaMemcpy(hDevMem.G_eta, G_eta, sizeReals * 2, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.G_u,   G_u,   sizeReals * 2, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.G_v,   G_v,   sizeReals * 2, cudaMemcpyHostToDevice) );
	
	checkCudaErrors( cudaMemcpy(hDevMem.SWM_eta, SWM_eta, sizeReals * 2, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.SWM_u,   SWM_u,   sizeReals * 2, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.SWM_v,   SWM_v,   sizeReals * 2, cudaMemcpyHostToDevice) );
	
	checkCudaErrors( cudaMemcpy(hDevMem.F_eta, F_eta, sizeReals, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.F_x,   F_x,   sizeReals, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.F_y,   F_y,   sizeReals, cudaMemcpyHostToDevice) );
	
	checkCudaErrors( cudaMemcpy(hDevMem.H_u,   H_u,   sizeReals, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.H_v,   H_v,   sizeReals, cudaMemcpyHostToDevice) );
	
	checkCudaErrors( cudaMemcpy(hDevMem.cos_lat_v,   cos_lat_v,   hConstParams.Ny * sizeof(fortReal8), cudaMemcpyHostToDevice) );
}

//=============================================================================
void cu_setforcing_(fortReal8* F_x, fortReal8* F_y, fortReal8* F_eta)
{
	if(state != initialized)
		return;

	if(hDevMem.F_eta == 0 || hDevMem.F_x == 0 || hDevMem.F_y == 0)
		printf(" Cuda module: Warning: forcing fields not allocated in GPU memory.\n");
	
	const size_t fieldsize = hConstParams.Nx * hConstParams.Ny * sizeof(fortReal8);
#ifdef CUDA_USE_STREAMS
	checkCudaErrors( cudaMemcpyAsync(hDevMem.F_eta, F_eta, fieldsize, cudaMemcpyHostToDevice, stream[0]) );
	checkCudaErrors( cudaMemcpyAsync(hDevMem.F_x,   F_x,   fieldsize, cudaMemcpyHostToDevice, stream[1]) );
	checkCudaErrors( cudaMemcpyAsync(hDevMem.F_y,   F_y,   fieldsize, cudaMemcpyHostToDevice, stream[2]) );
#else
	checkCudaErrors( cudaMemcpy(hDevMem.F_eta, F_eta, fieldsize, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.F_x,   F_x,   fieldsize, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(hDevMem.F_y,   F_y,   fieldsize, cudaMemcpyHostToDevice) );	
#endif
}

//=============================================================================
void cu_advance_()
{	
	// flipflop
	hVarParams.N0 = (hVarParams.N0 + 1) % 2;
	hVarParams.N0p1 = (hVarParams.N0p1 + 1) % 2;
	hVarParams.NG0 = (hVarParams.NG0 + 1) % 2;
	hVarParams.NG0m1 = (hVarParams.NG0m1 + 1) % 2;
}

//=============================================================================
void cu_timestep_()
{
	if(state != initialized)
		return;
	
	static int itt = 0;
		
	if(state != initialized)
		return;
	
	if(itt++ < 2) {
		cudaDeviceSynchronize();
		cuExecuteEulerForwardKernel(grid, blocksize, stream);
	} else {
		cudaDeviceSynchronize();
		cuExecuteAdamBashForthKernel(grid, blocksize, stream);
	}
}

//=============================================================================
void cu_computestreamfunction_()
{
	if(state != initialized)
		return;

	cudaDeviceSynchronize(); // wait until calculation of u and v have finished
#ifdef  CALC_LIB_ELLIPTIC_SOLVER
	printf("Cuda Warning: CALC_LIB_ELLIPCTIC_SOLVER is defined but is is not implemented yet.\n");
#endif
	cuExecuteStreamFunctionKernel(grid, blocksize, stream);
}

//=============================================================================
void cu_copytohost_eta_(fortReal8 *h_eta)
{
	if(state != initialized)
		return;
	
	size_t sizeReals = hConstParams.Nx*hConstParams.Ny*sizeof(fortReal8);
#ifdef CUDA_USE_STREAMS
 	cudaStreamSynchronize(stream[0]); // wait for stream to complete calculation of eta
	checkCudaErrors( cudaMemcpyAsync(h_eta, hDevMem.SWM_eta, sizeReals * 2, cudaMemcpyDeviceToHost, stream[0]) );
#else
	checkCudaErrors( cudaMemcpy(h_eta, hDevMem.SWM_eta, sizeReals * 2, cudaMemcpyDeviceToHost) );
#endif
}

//=============================================================================
void cu_copytohost_u_(fortReal8 *h_u)
{
	if(state != initialized)
		return;
	
	size_t sizeReals = hConstParams.Nx*hConstParams.Ny*sizeof(fortReal8);
#ifdef CUDA_USE_STREAMS
	cudaStreamSynchronize(stream[1]); // wait for stream to complete calculation of u
	checkCudaErrors( cudaMemcpyAsync(h_u,   hDevMem.SWM_u,   sizeReals * 2, cudaMemcpyDeviceToHost, stream[1]) );
#else
	checkCudaErrors( cudaMemcpy(h_u, hDevMem.SWM_u, sizeReals * 2, cudaMemcpyDeviceToHost) );
#endif
}

//=============================================================================
void cu_copytohost_v_(fortReal8 *h_v)
{
	if(state != initialized)
		return;
	
	size_t sizeReals = hConstParams.Nx*hConstParams.Ny*sizeof(fortReal8);
#ifdef CUDA_USE_STREAMS
	cudaStreamSynchronize(stream[2]); // wait for stream to complete calculation of v
	checkCudaErrors( cudaMemcpyAsync(h_v,   hDevMem.SWM_v,   sizeReals * 2, cudaMemcpyDeviceToHost, stream[2]) );
#else
	checkCudaErrors( cudaMemcpy(h_v, hDevMem.SWM_v, sizeReals * 2, cudaMemcpyDeviceToHost) );
#endif
}

//=============================================================================
void cu_copytohost_psi_(fortReal8 *h_psi)
{
	if(state != initialized)
		return;
	
	size_t sizeReals = hConstParams.Nx*hConstParams.Ny*sizeof(fortReal8);
	cudaDeviceSynchronize(); // wait for all kernels to complete execution
	checkCudaErrors( cudaMemcpy(h_psi, hDevMem.diag_psi, sizeReals, cudaMemcpyDeviceToHost) );
}

//=============================================================================
void cu_finish_()
{	
	if(state == initialized) 
	{
		printf(" Freeing GPU resources\n");
#ifdef CUDA_USE_STREAMS
		for(int i = 0; i < 3; i++) {
			cudaStreamDestroy(stream[i]);
		}
#endif		
		cuFreeDataMemory();		
		// de allocate
		state = not_initialized;
	}
}

//=============================================================================
void cu_testsizes_(size_t *sizeInteger, size_t *sizeReal8)
{
	printf(" Cuda_module: cu_testsizes_()\n");
	// check types: fortran <-> C
	if(	*sizeInteger != sizeof(fortInteger)	 ||
		*sizeReal8   != sizeof(fortReal8)   )
	{
		printf(" Cuda module not available: Size of data types incompatible.\n");
		printf(" \t\tC\tF90\n\tInteger\t%lu\t%lu\n\tReal8\t%lu\t%lu\n", sizeof(fortInteger), *sizeInteger, sizeof(fortReal8), *sizeReal8);
		
		if(state == initialized)
			cu_finish_();
		
		state = not_available;
	}
	
#ifdef QUADRATIC_BOTTOM_FRICTION
	printf(" Warning: Cuda module uses linear bottom friction but quadratic friction is defined");
#endif	
}