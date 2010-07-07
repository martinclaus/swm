/*
Header file for shallow water model. Meant to set Input/Output constants
 */
#ifndef FILE_IO_SEEN
#define FILE_IO_SEEN

/* Namelist files */
#define MODEL_NL "model.namelist"
#define OUTPUT_NL "output.namelist"

/* Switches for output control */ 
!#define writeInput
#define TAXISNAME "TIME"
#define XAXISNAME "LONGITUDE"
#define YAXISNAME "LATITUDE"
#define XUNIT     "degrees_east"
#define YUNIT     "degrees_north"
#define TUNIT     "seconds"
#define OFILEETA   "eta_out.nc"
#define OFILEU     "u_out.nc"
#define OFILEV     "v_out.nc"
#define OFILEH     "H_out.nc"
#define OFILEFX    "fx_out.nc"
#define OFILEFY    "fy_out.nc"
#define OFILEPSI   "psi_out.nc"
#define OVARNAMEETA "ETA"
#define OVARNAMEU   "U"
#define OVARNAMEV   "V"
#define OVARNAMEH   "H"
#define OVARNAMEFX  "TAUX"
#define OVARNAMEFY  "TAUY"
#define OVARNAMEPSI "PSI"

#endif
