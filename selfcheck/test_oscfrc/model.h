/*
Header file for shallow water model. Meant to set different physics.
 */
#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN

/* OpenMP */
#define OMPCHUNK Nx
#define OMPSCHEDULE GUIDED

/* Switch for Shallow Water Model */
#define SWM
#define BAROTROPIC
#define LINEARISED_STATE_OF_REST
/* Switches for timestepping (only use one at a time) */
/* #define SWM_TSTEP_ADAMSBASHFORTH */
#define SWM_TSTEP_HEAPS
/* Switches for forcing */
/* #define TAU_SCALE 1e-2 */
/* Switches for damping */
/* #define LINEAR_BOTTOM_FRICTION
#define LATERAL_MIXING
#define NEWTONIAN_COOLING */
/* Sponge layers */
/*
#define NEWTONIAN_SPONGE "NSEW"
#define VELOCITY_SPONGE "NSEW"
*/

/* Switch to load dynamical variables from file (mean flow) */
/* #define DYNFROMFILE */

#define H_OVERWRITE

/* Possible calculation of Coriolis-Parameter */
#define CORIOLIS_NOF 0
#define CORIOLIS_FPLANE 1
#define CORIOLIS_BETAPLANE 2
#define CORIOLIS_SPHERICALGEOMETRY 3

#endif
