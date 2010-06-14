/*
Header file for shallow water model. Meant to set different physics
 */
#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN

/* OpenMP */
#define OMPCHUNK Nx
#define OMPSCHEDULE GUIDED

/* Switches for output control */ 
!#define writeInput

/* Switches for physics*/
!#define LINEAR_BOTTOM_FRICTION
#define QUADRATIC_BOTTOM_FRICTION
#define LATERAL_MIXING
!#define oldmixing

#endif
