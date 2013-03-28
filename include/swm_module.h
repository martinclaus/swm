/*
 * Header file for shallow water model.
 *  */
#ifndef FILE_SWM_SEEN
#define FILE_SWM_SEEN

/* Available units of sponge layer length scale new_sponge_efolding */
#define SCU_DEGREE 1
#define SCU_RADIUS_OF_DEFORMATION 2
#define SCU_METER 3

/* Unit used to scale new_sponge_efolding */
#define SPONGE_SCALE_UNIT SCU_DEGREE

/* Sponge layer cut-off distance in units of length scale, e.g. 3 equals to 3 times length scale */
#define SPONGE_CUT_OFF 15 

/* set default value for unit of sponge layer length scale */
#ifndef SPONGE_SCALE_UNIT
#define SPONGE_SCALE_UNIT SCU_RADIUS_OF_DEFORMATION
#endif

/* Maxuimum number of forcing variables */
#define SWM_MAX_FORCING_INPUT 16
/* Default chunksize for memoryChunk member of SWM_forcingStream objects */
#define SWM_DEF_FORCING_CHUNKSIZE 100
#endif
