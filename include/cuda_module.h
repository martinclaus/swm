#ifndef _CUDA_MODULE_H
#define _CUDA_MODULE_H

#include "cuda_fortran_types.h"

//=============================================================================
/*
 * Rules for functions called from Fortran program
 *
 * 1. append _ underscore to function name
 * 2. only lowercase letters allowed
 * 3. all arguments are passed as pointers
 * 
 * compiler does not throw any errors if wrong parameters are used
 */
extern "C"
{
	/* Initialization */
	void cu_init_();
	
	void cu_testsizes_(size_t *fortranIntegerSize, size_t *fortranReal8Size);
	
	/* Data Transfer: Host To Device */
	void cu_setconstants_(	
			fortInteger *pNx, fortInteger *pNy, fortReal8 *pDt,
			fortReal8 *pAB_C1, fortReal8 *pAB_C2, fortReal8 *A, 
			fortReal8 *dTheta, fortReal8 *dLambda
			);
	
	void cu_setfields_(
			bool *ocean_eta, bool *ocean_u, bool *ocean_v,
			fortReal8 *SWM_Coef_eta, fortReal8 *SWM_Coef_u, fortReal8 *SWM_Coef_v,
			fortReal8 *impl_eta, fortReal8 *impl_u, fortReal8 *impl_v,
			fortReal8 *G_eta, fortReal8 *G_u, fortReal8 *G_v,
			fortReal8 *SWM_eta, fortReal8 *SWM_u, fortReal8 *SWM_v,
			fortReal8 *F_eta, fortReal8 *F_x, fortReal8 *F_y,
			fortReal8 *H_u, fortReal8 *H_v,
			fortReal8 *cos_lat_v
			);

	/* Called each iteration */
	void cu_setforcing_(fortReal8* F_x, fortReal8* F_y, fortReal8* F_eta);
	
	void cu_timestep_();
      
	/* Advance Routine */
	void cu_advance_();
	
	/* used in Diag Module */
	void cu_computestreamfunction_();
	
	/* Data Transfer: Device To Host */
	void cu_copytohost_u_(fortReal8 *h_u);
	void cu_copytohost_v_(fortReal8 *h_v);
	void cu_copytohost_eta_(fortReal8 *h_eta);
	void cu_copytohost_psi_(fortReal8 *h_psi);
		
	/* Cleanup */
	void cu_finish_();
	
	/* Debug */
// 	void cu_compare_(fortReal8 *fort_eta, fortReal8 *fort_u, fortReal8 *fort_v, fortReal8 *fort_coef_eta, fortReal8 *fort_coef_u, fortReal8 *fort_coef_v);
}

#endif
