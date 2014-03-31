#ifndef CUDA_FORTRAN_TYPES_
#define CUDA_FORTRAN_TYPES_

#define CUDA_USE_STREAMS
//=============================================================================
// corresponding C types for Fortran types (platform specific)
// correctness is checked in cu_testsizes_() function at runtime
typedef 	double 		fortReal8;
typedef 	int 		fortInteger;

#endif // CUDA_FORTRAN_TYPES_