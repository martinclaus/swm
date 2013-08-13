/*
Header file for shallow water model. Meant to set Input/Output constants
 */
#ifndef FILE_DIAG_SEEN
#define FILE_DIAG_SEEN

#define DEF_NOUT_CHUNK 10000
#define DEF_OFILENAME  "out.nc"
#define DEF_OVARNAME   "var"
#define DEF_DIAG_TYPE  "S"
#define DEF_DIAG_FREQUENCY  1
#define DEF_DIAG_PERIOD     "1M"
#define DEF_DIAG_PROCESS    "A"

#define DVARNAME_PSI "PSI"
#define DVARNAME_CHI "VELPOT"
#define DVARNAME_U_ND "U_ND"
#define DVARNAME_V_ND "V_ND"

#endif

