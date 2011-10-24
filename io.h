/* TO-DO
  
  + switching on writeInput complains about OVARNAMEGAMMA_N not being known...
  --> Make sure that gamma_new is only written out if Newtonian cooling is used.
  --> Add relevant file and variable names here.

*/


/*
Header file for shallow water model. Meant to set Input/Output constants
 */
#ifndef FILE_IO_SEEN
#define FILE_IO_SEEN

/* Namelist files */
#define MODEL_NL "model.namelist"
#define OUTPUT_NL "output.namelist"
#define TRACER_NL "tracer.namelist"

/* Namlist unit identifier */
#define UNIT_MODEL_NL 17
#define UNIT_OUTPUT_NL 18
#define UNIT_TRACER_NL 19


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

/* If switched on, the diagnostics routines flush at every write_step */
#define DIAG_FLUSH

#endif
