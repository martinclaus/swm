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
/* Switches for timestepping (only use one at a time) */
!#define SWM_TSTEP_EULERFW
!#define SWM_TSTEP_HEAPS
#define SWM_TSTEP_ADAMSBASHFORTH
/* Switches for forcing */
!#define TDEP_FORCING
#define PERIODIC_FORCING_X SIN
/* Switches for Reynolds stress terms*/
!#define wo_u2_x_u
!#define wo_uv_y_u
!#define wo_uv_x_v
!#define wo_v2_y_v
/* Switches for damping */
!#define LINEAR_BOTTOM_FRICTION
!#define QUADRATIC_BOTTOM_FRICTION
#define LATERAL_MIXING
!#define oldmixing
#define NEWTONIAN_COOLING
/* Sponge layers */
#define NEWTONIAN_SPONGE_N
#define NEWTONIAN_SPONGE_S
#define VELOCITY_SPONGE_N
#define VELOCITY_SPONGE_S

/* Switch to load dynamical variables from file (instead of using SWM) */
!#define DYNFROMFILE

/* Switch for tracer module */
!#define TRACER

#ifdef ISSELFCHECK
#define LINEAR_BOTTOM_FRICTION
#define QUADRATIC_BOTTOM_FRICTION
#define LATERAL_MIXING
#define NEWTONIAN_COOLING
#define PERIODIC_FORCING
#define NEWTONIAN_SPONGE_N
#define NEWTONIAN_SPONGE_S
#define VELOCITY_SPONGE_N
#define VELOCITY_SPONGE_S
#define SWM
#define DYNFROMFILE
#define TDEP_FORCING
#define TRACER
#endif
#endif
