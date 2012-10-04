/*
Header file for tracer module.
*/
#ifndef FILE_TRACER_SEEN
#define FILE_TRACER_SEEN

/* use openMP parallelisation */
#define TRC_PARALLEL

/* remove divergent part from velocity field */
!#define TRC_CORRECT_VELOCITY_FIELD

/* time stepping scheme for tracer equation (define only one!)*/
!#define TRC_TSTEP_LEAPFROG 
#define TRC_TSTEP_ADAMSBASHFORTH

/* properties of time stepping schemes (do not edit!) */
#ifdef TRC_TSTEP_LEAPFROG
#define TRC_NLEVEL 3
#define TRC_NLEVEL0 2
#endif
#ifdef TRC_TSTEP_ADAMSBASHFORTH
#define TRC_NLEVEL 2
#define TRC_NLEVEL0 1
#endif


#endif
