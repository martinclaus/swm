/*
Header file for shallow water model. Meant to set different physics
 */
#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN

/* OpenMP */
#define OMPCHUNK Nx
#define OMPSCHEDULE GUIDED

/* Switches for physics*/
!#define LINEAR_BOTTOM_FRICTION
#define QUADRATIC_BOTTOM_FRICTION
#define LATERAL_MIXING
!#define oldmixing

/* Switches for Reynolds stress terms*/
!#define wo_u2_x_u
!#define wo_uv_y_u
!#define wo_uv_x_v
!#define wo_v2_y_v

#endif
