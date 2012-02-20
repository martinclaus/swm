/*
Header file for shallow water model. Meant to set Input/Output constants
 */
#ifndef FILE_IO_SEEN
#define FILE_IO_SEEN

/* Maximal length of Character arrays */
#define CHARLEN 80
#define FULLREC_STRLEN 12

/* Default values for fileHandle Type */
#define DEF_NCID  -99999
#define DEF_VARID -99999
#define DEF_TIMEDID -99999
#define DEF_TIMEVID -99999
#define DEF_NREC -99999

#define NF90_NOTIMEDIM -1

/* Define missing values in CDF Outfiles */
#define MISS_VAL_DEF 1e37

/* Namelist files */
#define MODEL_NL "model.namelist"
#define OUTPUT_NL "output.namelist"
#define TRACER_NL "model.namelist"
#define DYNFROMFILE_NL "model.namelist"

/* Namlist unit identifier */
#define UNIT_MODEL_NL 17
#define UNIT_OUTPUT_NL 18
#define UNIT_TRACER_NL 19
#define UNIT_DYNFROMFILE_NL 20


/* Switches for output control */ 
!#define WRITEINPUT
#define WRITEMEAN
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
#define OFILETRACER "C1.nc"
#define OFILEGAMMA_N "gamma_n"
#define OFILEETAMEAN "eta_mean_out.nc"
#define OFILEUMEAN "u_mean_out.nc"
#define OFILEVMEAN "v_mean_out.nc"
#define OFILEPSIMEAN "psi_mean_out.nc"
#define OFILEETA2MEAN "eta2_mean_out.nc"
#define OFILEU2MEAN "u2_mean_out.nc"
#define OFILEV2MEAN "v2_mean_out.nc"
#define OFILEPSI2MEAN "psi2_mean_out.nc"
#define OVARNAMEETA "ETA"
#define OVARNAMEU   "U"
#define OVARNAMEV   "V"
#define OVARNAMEH   "H"
#define OVARNAMEETAMEAN "ETAMEAN"
#define OVARNAMEUMEAN "UMEAN"
#define OVARNAMEVMEAN "VMEAN"
#define OVARNAMEPSIMEAN "PSIMEAN"
#define OVARNAMEETA2MEAN "ETA2MEAN"
#define OVARNAMEU2MEAN "U2MEAN"
#define OVARNAMEV2MEAN "V2MEAN"
#define OVARNAMEPSI2MEAN "PSI2MEAN"
#define OVARNAMEFX  "TAUX"
#define OVARNAMEFY  "TAUY"
#define OVARNAMEPSI "PSI"
#define OVARNAMETRACER "C"
#define OVARNAMEGAMMA_N "GAMMA_N"

/* If switched on, the diagnostics routines flush at every write_step */
#define DIAG_FLUSH

#endif
