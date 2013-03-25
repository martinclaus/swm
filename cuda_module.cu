#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cuda_module.h"

// simple kernel function that adds two vectors
__global__ void vect_add(float *a, float *b, int N)
{
   int idx = threadIdx.x;
   if (idx<N) a[idx] = a[idx] + b[idx];
}


void cu_init_()
{
	if(state == not_initialized) {
		
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
		
        printf(" Using GPU Device %d: %s with compute capability %d.%d (%d Multiprocessors)\n", 
        		devId, deviceProp.name, deviceProp.major, deviceProp.minor, deviceProp.multiProcessorCount);
		state = initialized;
	}
}

void cu_testsizes_(size_t *sizeInteger, size_t *sizeReal8)
{
	// check types: fortran <-> C
	if(		*sizeInteger != sizeof(fortInteger)	||
			*sizeReal8   != sizeof(fortReal8)   )
	{
		printf(" Cuda module not available: Size of data types incompatible.\n");
		printf(" \t\tC\tF90\n\tInteger\t%lu\t%lu\n\tReal8\t%lu\t%lu\n", sizeof(fortInteger), *sizeInteger, sizeof(fortReal8), *sizeReal8);
		if(state == initialized)
			cu_finish_();
		
		state = not_available;
	}
}

void cu_setVars_(	fortReal8 *PI, fortReal8 *D2R, fortReal8 *A, fortReal8 *OMEGA,
					fortReal8 *G, fortReal8 *RHO0, fortReal8 *r, fortReal8 *k,
					fortReal8 *Ah, fortReal8 *missval, 
					fortReal8 *gamma_new, fortReal8 *gamma_new_sponge,
					fortReal8 *new_sponge_efolding,
					fortInteger *Nx, fortInteger *Ny, fortInteger *Nt,
					fortReal8 *dt, fortReal8 *dLambda, fortReal8 *dTheta)
{
/*  allocate(u(1:Nx, 1:Ny, 1:Ns))
  allocate(v(1:Nx, 1:Ny, 1:Ns))
  allocate(eta(1:Nx, 1:Ny, 1:Ns))
  allocate(H(1:Nx, 1:Ny))
  allocate(H_u(1:Nx, 1:Ny))
  allocate(H_v(1:Nx, 1:Ny))
  allocate(H_eta(1:Nx, 1:Ny))
  allocate(land_H(1:Nx, 1:Ny))
  allocate(land_u(1:Nx, 1:Ny))
  allocate(land_v(1:Nx, 1:Ny))
  allocate(land_eta(1:Nx, 1:Ny))
  allocate(ocean_H(1:Nx, 1:Ny))
  allocate(ocean_u(1:Nx, 1:Ny))
  allocate(ocean_v(1:Nx, 1:Ny))
  allocate(ocean_eta(1:Nx, 1:Ny))
  allocate(ip1(1:Nx))
  allocate(im1(1:Nx))
  allocate(jp1(1:Ny))
  allocate(jm1(1:Ny))
  allocate(lat_eta(1:Ny))
  allocate(lat_u(1:Ny))
  allocate(lat_v(1:Ny))
  allocate(lat_H(1:Ny))
  allocate(lon_eta(1:Nx))
  allocate(lon_u(1:Nx))
  allocate(lon_v(1:Nx))
  allocate(lon_H(1:Nx))
  allocate(cosTheta_v(1:Ny))
  allocate(cosTheta_u(1:Ny))
  allocate(tanTheta_v(1:Ny))
  allocate(tanTheta_u(1:Ny))		
	 */
}

/*
call CU_setDomain(im1, ip1, jm1, jp1, &
                  lat_u, lat_v, lat_eta, lat_h, &
                  lon_u, lon_v, lon_eta, lon_h, &
                  cosTheta_u, cosTheta_v, &
                  tanTheta_u, tanTheta_v, &
                  H_u, H_v, H_eta, H, &
                  land_u, land_v, land_eta, land_h, &
                  ocean_u, ocean_v, ocean_eta, ocean_h)

*/
void cu_setdomain_()
{
	if(state != initialized)
		return;
	
	// copy
	printf(" Cuda Memcpy took 123ms\n");
}

void cu_advance_()
{
	if(state != initialized)
		return;
	
}

void cu_timestep_()
{
	if(state != initialized)
		return;
	
}

void kernel_wrapper_(float *a, float *b, int *Np)
{
	if(state != initialized)
		return;
	
	float  *a_d, *b_d;  // declare GPU vector copies
   
   int blocks = 1;     // uses 1 block of
   int N = *Np;        // N threads on GPU

   // Allocate memory on GPU
   cudaMalloc( (void **)&a_d, sizeof(float) * N );
   cudaMalloc( (void **)&b_d, sizeof(float) * N );

   // copy vectors from CPU to GPU
   cudaMemcpy( a_d, a, sizeof(float) * N, cudaMemcpyHostToDevice );
   cudaMemcpy( b_d, b, sizeof(float) * N, cudaMemcpyHostToDevice );

   // call function on GPU
   vect_add<<< blocks, N >>>( a_d, b_d, N);

   // copy vectors back from GPU to CPU
   cudaMemcpy( a, a_d, sizeof(float) * N, cudaMemcpyDeviceToHost );
   cudaMemcpy( b, b_d, sizeof(float) * N, cudaMemcpyDeviceToHost );

   // free GPU memory
   cudaFree(a_d);
   cudaFree(b_d);
   return;
}


void cu_copytohost_u_(fortReal8 *h_u)
{
	printf( " called copyToHost_u() with %d = %f\n", h_u, h_u[0]);
}

void cu_copytohost_v_(fortReal8 *h_v)
{
	printf( " called copyToHost_v()with %d = %f\n", h_v, h_v[0]);
}

void cu_copytohost_eta_(fortReal8 *h_eta)
{
	printf( " called copyToHost_eta() with %d = %f\n", h_eta, h_eta[0]);
}

void cu_finish_()
{
	if(state == initialized) {
		printf(" Freeing GPU resources\n");
		// de allocate
		state = not_initialized;
	}
}
