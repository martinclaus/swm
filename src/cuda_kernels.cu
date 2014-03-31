#include "cuda_kernels.h"
#include "cuda_memory.h"
#include <cuda.h>
#include <stdio.h>
#include "model.h"

//================== Kernel Declaration: AdamBashforth, Euler, Psi ==============
__global__ void adamBashforth_eta(const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem);
__global__ void adamBashforth_u(  const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem);
__global__ void adamBashforth_v(  const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem);

__global__ void eulerForward_eta(const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem);
__global__ void eulerForward_u(  const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem);
__global__ void eulerForward_v(  const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem);

__global__ void computeStreamFunctionA(  const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem);
__global__ void computeStreamFunctionB(  const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem);
__global__ void computeStreamFunctionC(  const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem);


//=============================================================================
void cuExecuteAdamBashForthKernel(dim3 const grid, dim3 const blocksize, cudaStream_t const *stream) {
#ifdef CUDA_USE_STREAMS
	adamBashforth_eta<<< grid, blocksize, 0, stream[0] >>>(hConstParams, hVarParams, hDevMem);
	adamBashforth_u  <<< grid, blocksize, 0, stream[1] >>>(hConstParams, hVarParams, hDevMem);
	adamBashforth_v  <<< grid, blocksize, 0, stream[2] >>>(hConstParams, hVarParams, hDevMem);
#else
	adamBashforth_eta<<< grid, blocksize>>>(hConstParams, hVarParams, hDevMem);
	adamBashforth_u  <<< grid, blocksize>>>(hConstParams, hVarParams, hDevMem);
	adamBashforth_v  <<< grid, blocksize>>>(hConstParams, hVarParams, hDevMem);	
#endif
}

//=============================================================================
void cuExecuteEulerForwardKernel(dim3 const grid, dim3 const blocksize, cudaStream_t const *stream) {
#ifdef CUDA_USE_STREAMS
	eulerForward_eta<<<grid, blocksize, 0, stream[0] >>>(hConstParams, hVarParams, hDevMem);
	eulerForward_u<<<grid, blocksize, 0, stream[1] >>>(hConstParams, hVarParams, hDevMem);
	eulerForward_v<<<grid, blocksize, 0, stream[2] >>>(hConstParams, hVarParams, hDevMem);
#else
	eulerForward_eta<<<grid, blocksize>>>(hConstParams, hVarParams, hDevMem);
	eulerForward_u<<<grid, blocksize>>>(hConstParams, hVarParams, hDevMem);
	eulerForward_v<<<grid, blocksize>>>(hConstParams, hVarParams, hDevMem);
#endif	
}

//=============================================================================
void cuExecuteStreamFunctionKernel(dim3 const grid, dim3 const blocksize, cudaStream_t const *streams)
{
	cudaDeviceSynchronize();
	computeStreamFunctionA<<<grid, blocksize>>>(hConstParams, hVarParams, hDevMem);
	computeStreamFunctionB<<<grid, blocksize>>>(hConstParams, hVarParams, hDevMem);
	computeStreamFunctionC<<<grid, blocksize>>>(hConstParams, hVarParams, hDevMem);
}

//=============================================================================
// Makros to access arrays in fortran style
// Array access is 'circular', meaning index violations non-existent
// Constant input:
#define ocean_eta(i,j) 		devMem.ocean_eta[((j) % Ny) * Nx + ((i) % Nx)]
#define ocean_u(i,j) 		devMem.ocean_u[((j) % Ny) * Nx + ((i) % Nx)]
#define ocean_v(i,j) 		devMem.ocean_v[((j) % Ny) * Nx + ((i) % Nx)]
#define impl_eta(i,j)		devMem.impl_eta[((j) % Ny) * Nx + ((i) % Nx)]
#define impl_u(i,j)		devMem.impl_u[((j) % Ny) * Nx + ((i) % Nx)]
#define impl_v(i,j)		devMem.impl_v[((j) % Ny) * Nx + ((i) % Nx)]
#define SWM_Coef_eta(i,j,k) 	devMem.SWM_Coef_eta[(k) * Nx * 9 + ((j) % Ny) * 9 + (i)]
#define SWM_Coef_u(i,j,k) 	devMem.SWM_Coef_u[(k) * Nx * 11 + ((j) % Ny) * 11 + (i)]
#define SWM_Coef_v(i,j,k) 	devMem.SWM_Coef_v[(k) * Nx * 11 + ((j) % Ny) * 11 + (i)]
#define cos_lat_v(j)		devMem.cos_lat_v[j % Ny]
#define H_v(i,j)		devMem.H_v[((j) % Ny) * Nx + ((i) % Nx)]
#define H_u(i,j)		devMem.H_u[((j) % Ny) * Nx + ((i) % Nx)]

// Externally changed input:
#define F_eta(i,j) 		devMem.F_eta[((j) % Ny) * Nx + ((i) % Nx)]
#define F_x(i,j) 		devMem.F_x[((j) % Ny) * Nx + ((i) % Nx)]
#define F_y(i,j) 		devMem.F_y[((j) % Ny) * Nx + ((i) % Nx)]

// Cuda Output:
#define G_eta(i,j,k) 		devMem.G_eta[(k) * Nx * Ny + ((j) % Ny) * Nx + ((i) % Nx)]
#define G_u(i,j,k) 		devMem.G_u[(k) * Nx * Ny + ((j) % Ny) * Nx + ((i) % Nx)]
#define G_v(i,j,k) 		devMem.G_v[(k) * Nx * Ny + ((j) % Ny) * Nx + ((i) % Nx)]
#define SWM_eta(i,j,k)		devMem.SWM_eta[(k) * Nx * Ny + ((j) % Ny) * Nx + ((i) % Nx)]
#define SWM_u(i,j,k) 		devMem.SWM_u[(k) * Nx * Ny + ((j) % Ny) * Nx + ((i) % Nx)]
#define SWM_v(i,j,k) 		devMem.SWM_v[(k) * Nx * Ny + ((j) % Ny) * Nx + ((i) % Nx)]
#define psi(i,j)		devMem.diag_psi[((j) % Ny) * Nx + ((i) % Nx)]

//=============================================================================
__global__ void adamBashforth_eta(const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem)
{
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int Nx = constParams.Nx;
	const unsigned int Ny = constParams.Ny;
	const unsigned int N0 = varParams.N0;
	const unsigned int N0p1 = varParams.N0p1;
	const unsigned int NG0 = varParams.NG0;
	const unsigned int NG0m1 = varParams.NG0m1;
	
	if( i < Nx && j < Ny && ocean_eta(i, j) == 1 )
	{
		G_eta(i, j, NG0) = 
				SWM_eta(i, j, N0) 	* SWM_Coef_eta(0,i,j) +
				SWM_eta(i + 1, 	j, N0)	* SWM_Coef_eta(1,i,j) +
				SWM_eta(i - 1, 	j, N0) 	* SWM_Coef_eta(2,i,j) +
				SWM_eta(i, j + 1, N0) 	* SWM_Coef_eta(3,i,j) +
				SWM_eta(i, j - 1, N0) 	* SWM_Coef_eta(4,i,j) +
				SWM_u(	i + 1, 	j, N0) 	* SWM_Coef_eta(5,i,j) +
				SWM_u(	i, j, N0) 	* SWM_Coef_eta(6,i,j) +
				SWM_v(	i, j + 1, N0) 	* SWM_Coef_eta(7,i,j) +
				SWM_v(	i, j, N0) 	* SWM_Coef_eta(8,i,j) +
				F_eta(	i, j);
		  
		SWM_eta(i, j, N0p1) =
				(SWM_eta(i, j, N0) + constParams.dt * (constParams.AB_C1 * G_eta(i, j, NG0) - constParams.AB_C2 * G_eta(i, j, NG0m1)))/impl_eta(i,j);
	}
}

//=============================================================================
__global__ void adamBashforth_u(const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem)
{
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int Nx = constParams.Nx;
	const unsigned int Ny = constParams.Ny;
	const unsigned int N0 = varParams.N0;
	const unsigned int N0p1 = varParams.N0p1;
	const unsigned int NG0 = varParams.NG0;
	const unsigned int NG0m1 = varParams.NG0m1;
    
	if( i < Nx && j < Ny && ocean_u(i, j) == 1 )
	{
		G_u(i, j, NG0) = 
				SWM_u(	i, 	j, N0) 		* SWM_Coef_u(0, i, j) +
				SWM_u(	i + 1, 	j, N0) 		* SWM_Coef_u(1, i, j) +
				SWM_u(	i - 1, 	j, N0) 		* SWM_Coef_u(2, i, j) +
				SWM_u(	i, 	j + 1, N0) 	* SWM_Coef_u(3, i, j) +
				SWM_u(	i, 	j - 1, N0) 	* SWM_Coef_u(4, i, j) +
				SWM_v(	i, 	j, N0) 		* SWM_Coef_u(5, i, j) +
				SWM_v(	i - 1,	j, N0) 		* SWM_Coef_u(6, i, j) +
				SWM_v(	i - 1, 	j + 1, N0) 	* SWM_Coef_u(7, i, j) +
				SWM_v(	i, 	j + 1, N0) 	* SWM_Coef_u(8, i, j) +
				SWM_eta(i, 	j, N0)		* SWM_Coef_u(9,i, j) +
				SWM_eta(i - 1, 	j, N0)		* SWM_Coef_u(10,i, j) +
				F_x(i, j);
		
		SWM_u(i, j, N0p1) =
				(SWM_u(i, j, N0) + constParams.dt * (constParams.AB_C1 * G_u(i, j, NG0) - constParams.AB_C2 * G_u(i, j, NG0m1)))/impl_u(i,j);
	}
}

//=============================================================================
__global__ void adamBashforth_v(const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem)
{
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int Nx = constParams.Nx;
	const unsigned int Ny = constParams.Ny;
	const unsigned int N0 = varParams.N0;
	const unsigned int N0p1 = varParams.N0p1;
	const unsigned int NG0 = varParams.NG0;
	const unsigned int NG0m1 = varParams.NG0m1;
  
	if( i < Nx && j < Ny && ocean_v(i, j) == 1 )
	{
		G_v(i, j, NG0) = 
				SWM_v(	i, 	j, N0) 		* SWM_Coef_v(0, i, j) +
				SWM_v(	i + 1, 	j, N0) 		* SWM_Coef_v(1, i, j) +
				SWM_v(	i - 1, 	j, N0) 		* SWM_Coef_v(2, i, j) +
				SWM_v(	i, 	j + 1, N0) 	* SWM_Coef_v(3, i, j) +
				SWM_v(	i, 	j - 1, N0) 	* SWM_Coef_v(4, i, j) +
				SWM_u(	i + 1, 	j - 1, N0) 	* SWM_Coef_v(5, i, j) +
				SWM_u(	i, 	j - 1, N0) 	* SWM_Coef_v(6, i, j) +
				SWM_u(	i, 	j, N0) 		* SWM_Coef_v(7, i, j) +
				SWM_u(	i + 1, 	j, N0) 		* SWM_Coef_v(8, i, j) +
				SWM_eta(i, 	j, N0)		* SWM_Coef_v(9,i, j) +
				SWM_eta(i, 	j - 1, N0)	* SWM_Coef_v(10,i, j) +
				F_y(i, j);
		
		SWM_v(i, j, N0p1) =
				(SWM_v(i, j, N0) + constParams.dt * (constParams.AB_C1 * G_v(i, j, NG0) - constParams.AB_C2 * G_v(i, j, NG0m1)))/impl_v(i,j);
	}
}

//=============================================================================
// Kernels for euler forward iteration
//=============================================================================
__global__ void eulerForward_eta(const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem)
{
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int Nx = constParams.Nx;
	const unsigned int Ny = constParams.Ny;
	const unsigned int N0 = varParams.N0;
	const unsigned int N0p1 = varParams.N0p1;
	const unsigned int NG0 = varParams.NG0;
    
	if( i < Nx && j < Ny && ocean_eta(i, j) == 1 )
	{
		G_eta(i, j, NG0) = 
				SWM_eta(i, j, N0) 	* SWM_Coef_eta(0,i,j) +
				SWM_eta(i + 1, j, N0) 	* SWM_Coef_eta(1,i,j) +
				SWM_eta(i - 1, j, N0) 	* SWM_Coef_eta(2,i,j) +
				SWM_eta(i, j + 1, N0)	* SWM_Coef_eta(3,i,j) +
				SWM_eta(i, j - 1, N0) 	* SWM_Coef_eta(4,i,j) +
				SWM_u(i + 1, j, N0) 	* SWM_Coef_eta(5,i,j) +
				SWM_u(i, j, N0)		* SWM_Coef_eta(6,i,j) +
				SWM_v(i, j + 1, N0) 	* SWM_Coef_eta(7,i,j) +
				SWM_v(i, j, N0) 	* SWM_Coef_eta(8,i,j) +
				F_eta(i, j);
		
		SWM_eta(i, j, N0p1) = (SWM_eta(i, j, N0) + constParams.dt * G_eta(i, j, NG0))/impl_eta(i,j);
	}
}

//=============================================================================
__global__ void eulerForward_u(const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem)
{
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int Nx = constParams.Nx;
	const unsigned int Ny = constParams.Ny;
	const unsigned int N0 = varParams.N0;
	const unsigned int N0p1 = varParams.N0p1;
	const unsigned int NG0 = varParams.NG0;
    
	if( i < Nx && j < Ny && ocean_u(i, j) == 1 )
	{
		G_u(i, j, NG0) = 
				SWM_u(i, j, N0) 	* SWM_Coef_u(0,i,j) +
				SWM_u(i + 1, j, N0) 	* SWM_Coef_u(1,i,j) +
				SWM_u(i - 1, j, N0) 	* SWM_Coef_u(2,i,j) +
				SWM_u(i, j + 1, N0)	* SWM_Coef_u(3,i,j) +
				SWM_u(i, j - 1, N0) 	* SWM_Coef_u(4,i,j) +
				SWM_v(i, j, N0) 	* SWM_Coef_u(5,i,j) +
				SWM_v(i - 1,j, N0) 	* SWM_Coef_u(6,i,j) +
				SWM_v(i - 1, j + 1, N0) * SWM_Coef_u(7,i,j) +
				SWM_v(i, j + 1, N0) 	* SWM_Coef_u(8,i,j) +
				SWM_eta(i, j, N0) 	* SWM_Coef_u(9,i,j) +
				SWM_eta(i - 1, j, N0) 	* SWM_Coef_u(10,i,j) +
				F_x(i, j);
		
		SWM_u(i, j, N0p1) = (SWM_u(i, j, N0) + constParams.dt * G_u(i, j, NG0))/impl_u(i,j);
	}
}    

//=============================================================================
__global__ void eulerForward_v(const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem)
{
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int Nx = constParams.Nx;
	const unsigned int Ny = constParams.Ny;
	const unsigned int N0 = varParams.N0;
	const unsigned int N0p1 = varParams.N0p1;
	const unsigned int NG0 = varParams.NG0;
    
	if( i < Nx && j < Ny && ocean_v(i, j) == 1 )
	{
		G_v(i, j, NG0) = 
				SWM_v(i, j, N0) 	* SWM_Coef_v(0,i,j) +
				SWM_v(i + 1, j, N0) 	* SWM_Coef_v(1,i,j) +
				SWM_v(i - 1, j, N0) 	* SWM_Coef_v(2,i,j) +
				SWM_v(i, j + 1, N0)	* SWM_Coef_v(3,i,j) +
				SWM_v(i, j - 1, N0) 	* SWM_Coef_v(4,i,j) +
				SWM_u(i + 1, j - 1, N0) * SWM_Coef_v(5,i,j) +
				SWM_u(i,j - 1, N0) 	* SWM_Coef_v(6,i,j) +
				SWM_u(i, j, N0) 	* SWM_Coef_v(7,i,j) +
				SWM_u(i + 1, j, N0) 	* SWM_Coef_v(8,i,j) +
				SWM_eta(i, j, N0) 	* SWM_Coef_v(9,i,j) +
				SWM_eta(i, j - 1, N0) 	* SWM_Coef_v(10,i,j) +
				F_y(i, j);
		
		SWM_v(i, j, N0p1) = (SWM_v(i, j, N0) + constParams.dt * G_v(i, j, NG0))/impl_v(i,j);
	}
}
//=============================================================================

//=============================================================================
__global__ void computeStreamFunctionA(  const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem)
{
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int Nx = constParams.Nx;
	const unsigned int Ny = constParams.Ny;
	
	if( j == 0 && i < Nx - 1)
	{
		psi(i,j) = 0;
		for(int x = i; x < Nx - 1; x++) {
			psi(i,j) += 
#ifdef BAROTROPIC
					H_v(x, j) *
#endif
					ocean_v(x, j) * SWM_v(x, j, varParams.N0p1);    
		}
		psi(i,j) *= (-1) * constParams.A * cos_lat_v(j) * constParams.dLambda;
	}
}

//=============================================================================
__global__ void computeStreamFunctionB(  const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem)
{
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int Nx = constParams.Nx;
	const unsigned int Ny = constParams.Ny;
	
	if( i == Nx-1 && j < Ny && j > 0 )
	{
		psi(i,j) = 0;
		for(int y = 0; y < j; y++) {
			psi(i,j) += 
#ifdef BAROTROPIC
					H_u(i, y) *
#endif
					ocean_u(i, y) * SWM_u(i, y, varParams.N0p1);    
		}
		psi(i,j) *= (-1) * constParams.A * constParams.dTheta;
	}
}

//=============================================================================
__global__ void computeStreamFunctionC(  const constExecutionParameters constParams, const varExecutionParameters varParams, const devMemory devMem)
{
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int Nx = constParams.Nx;
	const unsigned int Ny = constParams.Ny;
	
	if( i < Nx-1 && j < Ny && j > 0 )
	{
		double sum1 = 0;
		for(int x = i; x < Nx-1; x++) {
			sum1 +=
#ifdef BAROTROPIC
					H_v(x, j) *
#endif
					SWM_v(x, j, varParams.N0p1);
		}
		sum1 *= (-1) * constParams.A * cos_lat_v(j) * constParams.dLambda;
				
		double sum2 = 0;
		for(int y = 0; y < j-1; y++) 	{
			sum2 +=
#ifdef BAROTROPIC
					H_u(i, y) *
#endif
					SWM_u(i, y, varParams.N0p1);
		}
		sum2 *= constParams.A * constParams.dTheta;
		
		psi(i,j) = (sum1 + psi(Nx-1, j) - sum2 + psi(i,0)) / 2.;
	}
}