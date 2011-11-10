/*
Header file for diagnosing module of shallow water model.
 */
#ifndef FILE_DIAG_MODULE_SEEN
#define FILE_DIAG_MODULE_SEEN

/* This should be defined in the Makefile */
/*#define DIAG_ELLIPTIC_SOLVER ElSolv_SOR */

#ifdef DIAG_ELLIPTIC_SOLVER
/* include headder of used Elliptic Solver Module and copy caller functions and module name */
#include DIAG_ELLIPTIC_SOLVER_HEADER /*defined in Makefile if DIAG_ELLIPTIC_SOLVER is defined there */
#define DIAG_ELLIPTIC_SOLVER_INIT FUNC_INIT_ELSOLV
#define DIAG_ELLIPTIC_SOLVER_FINISH FUNC_FINISH_ELSOLV
#define DIAG_ELLIPTIC_SOLVER_MAIN FUNC_ELSOLV
#define DIAG_ELLIPTIC_SOLVER_MODULE MODULE_ELSOLV
#endif

#endif

