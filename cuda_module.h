#ifndef CUDA_MODULE_H
#define CUDA_MODULE_H


enum module_state {
		not_initialized, 	// didn't try
		initialized,		// tried and succeeded
		not_available		// tried and failed
} state;


// corresponding c types for fortran types (platform specific)
typedef 	double 		fortReal8;
typedef 	int 		fortInteger;


/*
 * Procedures called from Fortran program
 *
 * 1. append _ underscore to function name
 * 2. do not use uppercase letters
 * 3. all arguments are passed as pointers
 */
extern "C"
{
	void cu_testsizes_(size_t *sizeInteger, size_t *sizeReal8);

	void cu_init_();

	void cu_setVars_(	fortReal8 *PI, fortReal8 *D2R, fortReal8 *A, fortReal8 *OMEGA,
						fortReal8 *G, fortReal8 *RHO0, fortReal8 *r, fortReal8 *k,
						fortReal8 *Ah, fortReal8 *missval,
						fortReal8 *gamma_new, fortReal8 *gamma_new_sponge,
						fortReal8 *new_sponge_efolding,
						fortInteger *Nx, fortInteger *Ny, fortInteger *Nt,
						fortReal8 *dt, fortReal8 *dLambda, fortReal8 *dTheta);

	void cu_setdomain_();

	void cu_advance_();

	void cu_timestep_();

	void cu_copytohost_u_(fortReal8 *h_u);

	void cu_copytohost_v_(fortReal8 *h_v);

	void cu_copytohost_eta_(fortReal8 *h_eta);

	void kernel_wrapper_(float *a, float *b, int *Np);

	void cu_finish_();
}

#endif
