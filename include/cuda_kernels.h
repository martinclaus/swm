#ifndef CUDA_KERNELS_H_
#define CUDA_KERNELS_H_

/* 
 * Execute adam bashforth scheme on the GPU for every point.
 * The grid and blocksize can be larger than the actual grid size for performance reasons.
 */
void cuExecuteAdamBashForthKernel(dim3 const grid, dim3 const blocksize, cudaStream_t const *streams);

/*
 * Execute the Euler Forward scheme on the GPU.
 * The grid and blocksize can be larger than the actual grid size for performance reasons.
 */
void cuExecuteEulerForwardKernel(dim3 const grid, dim3 const blocksize, cudaStream_t const *streams);

/*
 * Compute the values of psi on the GPU.
 * The grid and blocksize can be larger than the actual grid size for performance reasons.
 */
void cuExecuteStreamFunctionKernel(dim3 const grid, dim3 const blocksize, cudaStream_t const *streams);

#endif