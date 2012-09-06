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
#define SWM_TSTEP_ADAMSBASHFORTH
!#define SWM_TSTEP_HEAPS
/* Switches for forcing */
#define PERIODIC_FORCING_X SIN
#define TAU_SCALE 9.949307700452987e+01
/* Switches for damping */
#define LINEAR_BOTTOM_FRICTION
#define LATERAL_MIXING
#define NEWTONIAN_COOLING
/* Sponge layers */
#define VELOCITY_SPONGE_N
#define VELOCITY_SPONGE_S
#define NEWTONIAN_SPONGE_N
#define NEWTONIAN_SPONGE_S

/* Switch to load dynamical variables from file (mean flow) */
!#define DYNFROMFILE

#define DIAG_START 0.

#endif
