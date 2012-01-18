/*
Header file for calc_lib module of shallow water model.
 */
#ifndef FILE_CALC_LIB_MODULE_SEEN
#define FILE_CALC_LIB_MODULE_SEEN

#define CORRECT_FLOW_FOR_PSI /* Flag if non-divergent flow should be computed */

/* This should be defined in the Makefile */
/*#define CALC_LIB_ELLIPTIC_SOLVER ElSolv_SOR */

#ifdef CALC_LIB_ELLIPTIC_SOLVER
/* include headder of used Elliptic Solver Module and copy caller functions and module name */
#include CALC_LIB_ELLIPTIC_SOLVER_HEADER /*defined in Makefile if CALC_ELLIPTIC_SOLVER is defined there */
#define CALC_LIB_ELLIPTIC_SOLVER_INIT FUNC_INIT_ELSOLV
#define CALC_LIB_ELLIPTIC_SOLVER_FINISH FUNC_FINISH_ELSOLV
#define CALC_LIB_ELLIPTIC_SOLVER_MAIN FUNC_ELSOLV
#define CALC_LIB_ELLIPTIC_SOLVER_MODULE MODULE_ELSOLV
#endif

#endif
