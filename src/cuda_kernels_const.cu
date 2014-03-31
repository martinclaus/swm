//=============================================================================
// Device memory: Put execution parameters into special constant memory region for faster access
// Const Parameters
__constant__  fortReal8 dc_dt;
__constant__  fortReal8 dc_AB_C1;
__constant__  fortReal8 dc_AB_C2;
__constant__  fortInteger dc_Nx;
__constant__  fortInteger dc_Ny;

// Variable Parameters
__constant__  int dc_NG0;
__constant__  int dc_NG0m1;
__constant__  int dc_N0;
__constant__  int dc_N0p1;

// DevMem
__constant__  fortInteger *dc_ocean_eta;
__constant__  fortInteger *dc_ocean_u;
__constant__  fortInteger *dc_ocean_v;
__constant__  fortReal8 *dc_SWM_Coef_eta;
__constant__  fortReal8 *dc_SWM_Coef_u;
__constant__  fortReal8 *dc_SWM_Coef_v;
__constant__  fortReal8 *dc_impl_eta;
__constant__  fortReal8 *dc_impl_u;
__constant__  fortReal8 *dc_impl_v;
__constant__  fortReal8 *dc_G_eta;
__constant__  fortReal8 *dc_G_u;
__constant__  fortReal8 *dc_G_v;
__constant__  fortReal8 *dc_SWM_eta;
__constant__  fortReal8 *dc_SWM_u;
__constant__  fortReal8 *dc_SWM_v;
__constant__  fortReal8 *dc_F_eta;
__constant__  fortReal8 *dc_F_x;
__constant__  fortReal8 *dc_F_y;

//=============================================================================
void cuOptimizeConstParams()
{
	checkCudaErrors( cudaMemcpyToSymbol(dc_dt, &hConstParams.dt, sizeof(fortReal8)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_AB_C1, &hConstParams.AB_C1, sizeof(fortReal8)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_AB_C2, &hConstParams.AB_C2, sizeof(fortReal8)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_Nx, &hConstParams.Nx, sizeof(fortInteger)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_Ny, &hConstParams.Ny, sizeof(fortInteger)) );
}

//=============================================================================
void cuOptimizeVarParams()
{
	checkCudaErrors( cudaMemcpyToSymbol(dc_NG0, &hVarParams.NG0, sizeof(int)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_NG0m1, &hVarParams.NG0m1, sizeof(int)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_N0, &hVarParams.N0, sizeof(int)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_N0p1, &hVarParams.N0p1, sizeof(int)) );
}

//=============================================================================
void cuOptimizeDevMem()
{
	checkCudaErrors( cudaMemcpyToSymbol(dc_ocean_eta, 	&hDevMem.ocean_eta, 	sizeof(fortInteger*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_ocean_u, 	&hDevMem.ocean_u, 	sizeof(fortInteger*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_ocean_v, 	&hDevMem.ocean_v, 	sizeof(fortInteger*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_SWM_Coef_eta, 	&hDevMem.SWM_Coef_eta, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_SWM_Coef_u, 	&hDevMem.SWM_Coef_u, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_SWM_Coef_v, 	&hDevMem.SWM_Coef_v, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_impl_eta, 	&hDevMem.impl_eta, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_impl_u, 	&hDevMem.impl_u, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_impl_v, 	&hDevMem.impl_v, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_G_eta, 		&hDevMem.G_eta, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_G_u, 		&hDevMem.G_u, 		sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_G_v, 		&hDevMem.G_v, 		sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_SWM_eta, 	&hDevMem.SWM_eta, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_SWM_u, 		&hDevMem.SWM_u, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_SWM_v, 		&hDevMem.SWM_v, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_F_eta, 		&hDevMem.F_eta, 	sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_F_x, 		&hDevMem.F_x, 		sizeof(fortReal8*)) );
	checkCudaErrors( cudaMemcpyToSymbol(dc_F_y, 		&hDevMem.F_y, 		sizeof(fortReal8*)) );
}

//=============================================================================
// Makros to access arrays in fortran style
// Array access is 'circular', meaning index violations non-existent
// Constant input:
#define ocean_eta(i,j) 		dc_ocean_eta[((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define ocean_u(i,j) 		dc_ocean_u[((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define ocean_v(i,j) 		dc_ocean_v[((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define impl_eta(i,j)		dc_impl_eta[((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define impl_u(i,j)		dc_impl_u[((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define impl_v(i,j)		dc_impl_v[((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define SWM_Coef_eta(i,j,k) 	dc_SWM_Coef_eta[(k) * dc_Nx * 9 + ((j) % dc_Ny) * 9 + (i)]
#define SWM_Coef_u(i,j,k) 	dc_SWM_Coef_u[(k) * dc_Nx * 11 + ((j) % dc_Ny) * 11 + (i)]
#define SWM_Coef_v(i,j,k) 	dc_SWM_Coef_v[(k) * dc_Nx * 11 + ((j) % dc_Ny) * 11 + (i)]

// Externally changed input:
#define F_eta(i,j) 		dc_F_eta[((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define F_x(i,j) 		dc_F_x[((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define F_y(i,j) 		dc_F_y[((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]

// Cuda Output:
#define G_eta(i,j,k) 		dc_G_eta[(k) * dc_Nx * dc_Ny + ((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define G_u(i,j,k) 		dc_G_u[(k) * dc_Nx * dc_Ny + ((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define G_v(i,j,k) 		dc_G_v[(k) * dc_Nx * dc_Ny + ((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define SWM_eta(i,j,k)		dc_SWM_eta[(k) * dc_Nx * dc_Ny + ((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define SWM_u(i,j,k) 		dc_SWM_u[(k) * dc_Nx * dc_Ny + ((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]
#define SWM_v(i,j,k) 		dc_SWM_v[(k) * dc_Nx * dc_Ny + ((j) % dc_Ny) * dc_Nx + ((i) % dc_Nx)]

//=============================================================================
// Kernels for adamBashforth timestep()
//=============================================================================
__global__ void adamBashforth_eta(/*struct constExecutionParameters const * const constParams, struct varExecutionParameters const * const varParams, struct devMemory const * const devMem*/)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
    
	if( i < dc_Nx && j < dc_Ny && ocean_eta(i, j) == 1 )
	{
		G_eta(i, j, dc_NG0) = 
				SWM_eta(i, j, dc_N0) 		* SWM_Coef_eta(0,i,j) +
				SWM_eta(i + 1, 	j, dc_N0)	* SWM_Coef_eta(1,i,j) +
				SWM_eta(i - 1, 	j, dc_N0) 	* SWM_Coef_eta(2,i,j) +
				SWM_eta(i, j + 1, dc_N0) 	* SWM_Coef_eta(3,i,j) +
				SWM_eta(i, j - 1, dc_N0) 	* SWM_Coef_eta(4,i,j) +
				SWM_u(	i + 1, 	j, dc_N0) 	* SWM_Coef_eta(5,i,j) +
				SWM_u(	i, j, dc_N0) 		* SWM_Coef_eta(6,i,j) +
				SWM_v(	i, j + 1, dc_N0) 	* SWM_Coef_eta(7,i,j) +
				SWM_v(	i, j, dc_N0) 		* SWM_Coef_eta(8,i,j) +
				F_eta(	i, j);
		  
		SWM_eta(i, j, dc_N0p1) =
				(SWM_eta(i, j, dc_N0) + dc_dt * (dc_AB_C1 * G_eta(i, j, dc_NG0) - dc_AB_C2 * G_eta(i, j, dc_NG0m1)))/impl_eta(i,j);
	}
}

//=============================================================================
__global__ void adamBashforth_u(/*struct constExecutionParameters const * const constParams, struct varExecutionParameters const * const varParams, struct devMemory const * const devMem*/)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
    
	if( i < dc_Nx && j < dc_Ny && ocean_u(i, j) == 1 )
	{
		G_u(i, j, dc_NG0) = 
				SWM_u(	i, 	j, dc_N0) 	* SWM_Coef_u(0, i, j) +
				SWM_u(	i + 1, 	j, dc_N0) 	* SWM_Coef_u(1, i, j) +
				SWM_u(	i - 1, 	j, dc_N0) 	* SWM_Coef_u(2, i, j) +
				SWM_u(	i, 	j + 1, dc_N0) 	* SWM_Coef_u(3, i, j) +
				SWM_u(	i, 	j - 1, dc_N0) 	* SWM_Coef_u(4, i, j) +
				SWM_v(	i, 	j, dc_N0) 	* SWM_Coef_u(5, i, j) +
				SWM_v(	i - 1,	j, dc_N0) 	* SWM_Coef_u(6, i, j) +
				SWM_v(	i - 1, 	j + 1, dc_N0) 	* SWM_Coef_u(7, i, j) +
				SWM_v(	i, 	j + 1, dc_N0) 	* SWM_Coef_u(8, i, j) +
				SWM_eta(i, 	j, dc_N0)	* SWM_Coef_u(9,i, j) +
				SWM_eta(i - 1, 	j, dc_N0)	* SWM_Coef_u(10,i, j) +
				F_x(i, j);
		
		SWM_u(i, j, dc_N0p1) =
				(SWM_u(i, j, dc_N0) + dc_dt * (dc_AB_C1 * G_u(i, j, dc_NG0) - dc_AB_C2 * G_u(i, j, dc_NG0m1)))/impl_u(i,j);
	}
}

//=============================================================================
__global__ void adamBashforth_v(/*struct constExecutionParameters const * const constParams, struct varExecutionParameters const * const varParams, struct devMemory const * const devMem*/)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
    
	if( i < dc_Nx && j < dc_Ny && ocean_v(i, j) == 1 )
	{
		G_v(i, j, dc_NG0) = 
				SWM_v(	i, 	j, dc_N0) 	* SWM_Coef_v(0, i, j) +
				SWM_v(	i + 1, 	j, dc_N0) 	* SWM_Coef_v(1, i, j) +
				SWM_v(	i - 1, 	j, dc_N0) 	* SWM_Coef_v(2, i, j) +
				SWM_v(	i, 	j + 1, dc_N0) 	* SWM_Coef_v(3, i, j) +
				SWM_v(	i, 	j - 1, dc_N0) 	* SWM_Coef_v(4, i, j) +
				SWM_u(	i + 1, 	j - 1, dc_N0) 	* SWM_Coef_v(5, i, j) +
				SWM_u(	i, 	j - 1, dc_N0) 	* SWM_Coef_v(6, i, j) +
				SWM_u(	i, 	j, dc_N0) 	* SWM_Coef_v(7, i, j) +
				SWM_u(	i + 1, 	j, dc_N0) 	* SWM_Coef_v(8, i, j) +
				SWM_eta(i, 	j, dc_N0)	* SWM_Coef_v(9,i, j) +
				SWM_eta(i, 	j - 1, dc_N0)	* SWM_Coef_v(10,i, j) +
				F_y(i, j);
		
		SWM_v(i, j, dc_N0p1) =
				(SWM_v(i, j, dc_N0) + dc_dt * (dc_AB_C1 * G_v(i, j, dc_NG0) - dc_AB_C2 * G_v(i, j, dc_NG0m1)))/impl_v(i,j);
	}
}

//=============================================================================
// Kernels for euler forward iteration
//=============================================================================
__global__ void eulerForward_eta(/*struct constExecutionParameters const * const constParams, struct varExecutionParameters const * const varParams, struct devMemory const * const devMem*/)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
    
	if( i < dc_Nx && j < dc_Ny && ocean_eta(i, j) == 1 )
	{
		G_eta(i, j, dc_NG0) = 
				SWM_eta(i, j, dc_N0) 	* SWM_Coef_eta(0,i,j) +
				SWM_eta(i + 1, j, dc_N0) * SWM_Coef_eta(1,i,j) +
				SWM_eta(i - 1, j, dc_N0) * SWM_Coef_eta(2,i,j) +
				SWM_eta(i, j + 1, dc_N0)* SWM_Coef_eta(3,i,j) +
				SWM_eta(i, j - 1, dc_N0) * SWM_Coef_eta(4,i,j) +
				SWM_u(i + 1, j, dc_N0) 	* SWM_Coef_eta(5,i,j) +
				SWM_u(i, j, dc_N0)	* SWM_Coef_eta(6,i,j) +
				SWM_v(i, j + 1, dc_N0) 	* SWM_Coef_eta(7,i,j) +
				SWM_v(i, j, dc_N0) 	* SWM_Coef_eta(8,i,j) +
				F_eta(i, j);
		
		SWM_eta(i, j, dc_N0p1) = (SWM_eta(i, j, dc_N0) + dc_dt * G_eta(i, j, dc_NG0))/impl_eta(i,j);
	}
}

//=============================================================================
__global__ void eulerForward_u(/*struct constExecutionParameters const * const constParams, struct varExecutionParameters const * const varParams, struct devMemory const * const devMem*/)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
    
	if( i < dc_Nx && j < dc_Ny && ocean_u(i, j) == 1 )
	{
		G_u(i, j, dc_NG0) = 
				SWM_u(i, j, dc_N0) 	* SWM_Coef_u(0,i,j) +
				SWM_u(i + 1, j, dc_N0) 	* SWM_Coef_u(1,i,j) +
				SWM_u(i - 1, j, dc_N0) 	* SWM_Coef_u(2,i,j) +
				SWM_u(i, j + 1, dc_N0)	* SWM_Coef_u(3,i,j) +
				SWM_u(i, j - 1, dc_N0) 	* SWM_Coef_u(4,i,j) +
				SWM_v(i, j, dc_N0) 	* SWM_Coef_u(5,i,j) +
				SWM_v(i - 1,j, dc_N0) 	* SWM_Coef_u(6,i,j) +
				SWM_v(i - 1, j + 1, dc_N0) * SWM_Coef_u(7,i,j) +
				SWM_v(i, j + 1, dc_N0) 	* SWM_Coef_u(8,i,j) +
				SWM_eta(i, j, dc_N0) 	* SWM_Coef_u(9,i,j) +
				SWM_eta(i - 1, j, dc_N0) * SWM_Coef_u(10,i,j) +
				F_x(i, j);
		
		SWM_u(i, j, dc_N0p1) = (SWM_u(i, j, dc_N0) + dc_dt * G_u(i, j, dc_NG0))/impl_u(i,j);
	}
}    

//=============================================================================
__global__ void eulerForward_v(/*struct constExecutionParameters const * const constParams, struct varExecutionParameters const * const varParams, struct devMemory const * const devMem*/)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
    
	if( i < dc_Nx && j < dc_Ny && ocean_v(i, j) == 1 )
	{
		G_v(i, j, dc_NG0) = 
				SWM_v(i, j, dc_N0) 	* SWM_Coef_v(0,i,j) +
				SWM_v(i + 1, j, dc_N0) 	* SWM_Coef_v(1,i,j) +
				SWM_v(i - 1, j, dc_N0) 	* SWM_Coef_v(2,i,j) +
				SWM_v(i, j + 1, dc_N0)	* SWM_Coef_v(3,i,j) +
				SWM_v(i, j - 1, dc_N0) 	* SWM_Coef_v(4,i,j) +
				SWM_u(i + 1, j - 1, dc_N0) * SWM_Coef_v(5,i,j) +
				SWM_u(i,j - 1, dc_N0) 	* SWM_Coef_v(6,i,j) +
				SWM_u(i, j, dc_N0) 	* SWM_Coef_v(7,i,j) +
				SWM_u(i + 1, j, dc_N0) 	* SWM_Coef_v(8,i,j) +
				SWM_eta(i, j, dc_N0) 	* SWM_Coef_v(9,i,j) +
				SWM_eta(i, j - 1, dc_N0) 	* SWM_Coef_v(10,i,j) +
				F_y(i, j);
		
		SWM_v(i, j, dc_N0p1) = (SWM_v(i, j, dc_N0) + dc_dt * G_v(i, j, dc_NG0))/impl_v(i,j);
	}
}
//=============================================================================
