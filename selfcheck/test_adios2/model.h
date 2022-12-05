/*
Header file for shallow water model. Meant to set different physics.
 */
#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN

/* OpenMP */
#define OMP_COLLAPSE(N) collapse(N-1)
#define OMPCHUNK 5
/*#define OMPSCHEDULE GUIDED*/
#define OMPSCHEDULE STATIC

/* Switch for Shallow Water Model */
#define SWM
#define BAROTROPIC
#define LINEARISED_STATE_OF_REST
/* Switches for timestepping (only use one at a time) */
#define SWM_TSTEP_ADAMSBASHFORTH
/* Switches for forcing */
/* #define TAU_SCALE 1e-2 */
/* Switches for damping */
/* #define LATERAL_MIXING */

/* Switch to load dynamical variables from file (mean flow) */
/* #define DYNFROMFILE */

#define H_OVERWRITE_DEF 0.

/* Possible calculation of Coriolis-Parameter */
#define CORIOLIS_NOF 0
#define CORIOLIS_FPLANE 1
#define CORIOLIS_BETAPLANE 2
#define CORIOLIS_SPHERICALGEOMETRY 3

#endif