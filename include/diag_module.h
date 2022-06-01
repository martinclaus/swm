/*
Header file for shallow water model. Meant to set Input/Output constants
 */
#ifndef FILE_DIAG_SEEN
#define FILE_DIAG_SEEN

#define DEF_NOUT_CHUNK 10000
#define DEF_OFILENAME "out.nc"
#define DEF_OVARNAME "var"
#define DEF_DIAG_TYPE "S"
#define DEF_DIAG_FREQUENCY 1
#define DEF_DIAG_PERIOD 2.628e6
#define DEF_DIAG_PROCESS "A"
#define DEF_DIAG_UNLIMITED_TDIM .True.

/* Diagnostic Variables (must be upper case!) */
#define DVARNAME_PSI "PSI"
#define DVARNAME_CHI "VELPOT"
#define DVARNAME_U_ND "U_ND"
#define DVARNAME_V_ND "V_ND"

/* Comma separated list of diagnostic variables which are initialised
   with default values */
#define DVARNAME_LIST DVARNAME_PSI, DVARNAME_CHI, DVARNAME_U_ND, DVARNAME_V_ND

#endif
