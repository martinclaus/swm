/*
Header file for shallow water model. Meant to set different physics.
 */
#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN

/* OpenMP */
#define OMPCHUNK Nx
#define OMPSCHEDULE GUIDED

/* Switches for physics*/
!#define LINEAR_BOTTOM_FRICTION
!#define QUADRATIC_BOTTOM_FRICTION
#define LATERAL_MIXING
!#define NEWTONIAN_COOLING
!#define oldmixing
#define PERIODIC_FORCING

/* Sponge layers */
#define NEWTONIAN_SPONGE_N
#define NEWTONIAN_SPONGE_S

/* Switches for Reynolds stress terms*/
!#define wo_u2_x_u
!#define wo_uv_y_u
!#define wo_uv_x_v
!#define wo_v2_y_v

/* Switches for time dependent forcing */
!#define TDEP_FORCING

/* Switch for tracer module */
!#define TRACER

#endif
