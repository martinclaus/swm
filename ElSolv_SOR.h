/*
Header file for elliptic solver.
 */

/* use openMP parallelisation */
#define ELSOLV_SOR_PARALLEL

#ifdef ELSOLV_SOR_PARALLEL
#include "model.h"
#endif

/* define function names and module name*/
#define MODULE_ELSOLV ElSolv_SOR
#define FUNC_INIT_ELSOLV init_ElSolv_SOR
#define FUNC_FINISH_ELSOLV finish_ElSolv_SOR
#define FUNC_ELSOLV main_ElSolv_SOR

