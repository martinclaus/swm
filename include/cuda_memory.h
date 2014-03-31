#ifndef CUDA_MEMORY_H_
#define CUDA_MEMORY_H_

#include "cuda_fortran_types.h"

//=============================================================================
// These functions are responsible for allocation and de-allocation of memory on the GPU 
// The size of the data fields are defined inside cuAllocateDataMemory()*/
// The addresses of the allocated data fields are stored in a global structure devMemory
void cuAllocateDataMemory();

//=============================================================================
void cuFreeDataMemory();

//=============================================================================
// This data is set in cu_setconstants_ defined in cuda_module.h
struct constExecutionParameters {
  fortReal8 	dt;
  fortReal8 	AB_C1;
  fortReal8 	AB_C2;
  fortReal8 	A;
  fortReal8 	dTheta;
  fortReal8 	dLambda;
  fortInteger	Nx;
  fortInteger	Ny;
};
extern constExecutionParameters hConstParams;

//=============================================================================
// This data is set once in cu_setconstants_ and then altered in calls to 
// cu_advance defined in cuda_module.h
struct varExecutionParameters {
  int NG0;
  int NG0m1;
  int N0;
  int N0p1;
};
extern varExecutionParameters hVarParams;

//=============================================================================
// Contains pointers to all GPU data fields both constant and variable data
struct devMemory {
  fortInteger *ocean_eta;
  fortInteger *ocean_u;
  fortInteger *ocean_v;
  fortReal8 *SWM_Coef_eta;
  fortReal8 *SWM_Coef_u;
  fortReal8 *SWM_Coef_v;
  fortReal8 *impl_eta;
  fortReal8 *impl_u;
  fortReal8 *impl_v;
  fortReal8 *G_eta;
  fortReal8 *G_u;
  fortReal8 *G_v;
  fortReal8 *SWM_eta;
  fortReal8 *SWM_u;
  fortReal8 *SWM_v;
  fortReal8 *diag_psi;
  fortReal8 *F_eta;
  fortReal8 *F_x;
  fortReal8 *F_y;
  fortReal8 *H_u;
  fortReal8 *H_v;
  fortReal8 *cos_lat_v;
};
extern devMemory hDevMem;

//=============================================================================
// This will output the proper CUDA error strings in the event that a CUDA host call returns an error
#define checkCudaErrors(err)           __checkCudaErrors (err, __FILE__, __LINE__)
void __checkCudaErrors( cudaError err, const char *file, const int line );

#endif
