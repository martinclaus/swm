/*
Header file for tracer module.
*/
#ifndef FILE_TRACER_SEEN
#define FILE_TRACER_SEEN

/* time stepping scheme for tracer equation (define only one!)*/
#define TRC_TSTEP_ADAMSBASHFORTH

/* properties of time stepping schemes (do not edit!) */
#ifdef TRC_TSTEP_ADAMSBASHFORTH
#define TRC_NLEVEL 2
#define TRC_NLEVEL0 1
#define TRC_GLEVEL 2
#endif


#endif
