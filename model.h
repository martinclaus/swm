/*
Header file for shallow water model. Meant to set different physics
 */
#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN

/* OpenMP */
#define OMPCHUNK Nx
#define OMPSCHEDULE GUIDED

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
