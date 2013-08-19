/*
Header file for shallow water model. Meant to set Input/Output constants
 */
#ifndef FILE_IO_SEEN
#define FILE_IO_SEEN

/* Maximal length of Character arrays */
#define CHARLEN 80
#define FULLREC_STRLEN 12

/* Default chunk size for MemChunk */
#define DEF_NT_CHUNKSIZE 100

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
#define OUTPUT_NL MODEL_NL
#define TRACER_NL MODEL_NL
#define DYNFROMFILE_NL MODEL_NL
#define SWM_FORCING_NL MODEL_NL
#define SWM_DAMPING_NL MODEL_NL
#define DIAG_NL MODEL_NL
#define DOMAIN_NL MODEL_NL
#define CALENDAR_NL MODEL_NL

/* Namlist unit identifier */
#define UNIT_MODEL_NL 17
#define UNIT_OUTPUT_NL 18
#define UNIT_TRACER_NL 19
#define UNIT_DYNFROMFILE_NL 20
#define UNIT_SWM_FORCING_NL 21
#define UNIT_SWM_DAMPING_NL 22
#define UNIT_DIAG_NL 23
#define UNIT_DOMAIN_NL 24
#define UNIT_CALENDAR_NL 25

/* Attribute conventions from NetCDF Users' Guide (10.06.2011) */
#define NUG_ATT_UNITS "units"
#define NUG_ATT_LONG_NAME "long_name"
#define NUG_ATT_FILL "_FillValue"
#define NUG_ATT_MISS "missing_value"
#define NUG_ATT_VALMIN "valid_min"
#define NUG_ATT_VALMAX "valid_max"
#define NUG_ATT_VALRANGE "valid_range"
#define NUG_ATT_SCALE "scale_factor"
#define NUG_ATT_OFFSET "add_offset"

/* Switches for output control */
#define TAXISNAME "TIME"
#define XAXISNAME "LONGITUDE"
#define YAXISNAME "LATITUDE"
#define IDAXISNAME "ID"
#define XUNIT     "degrees_east"
#define YUNIT     "degrees_north"
#define TUNIT     "seconds since 0000-01-01 00:00:00"
#define OFILEETA   "eta_out.nc"
#define OFILEU     "u_out.nc"
#define OFILEV     "v_out.nc"
#define OFILEH     "H_out.nc"
#define OFILEFX    "fx_out.nc"
#define OFILEFY    "fy_out.nc"
#define OFILEPSI   "psi_out.nc"
#define OFILETRACER "C1.nc"
#define OFILEGAMMA_N "gamma_n.nc"
#define OFILEGAMMA_U "gamma_u.nc"
#define OFILEGAMMA_V "gamma_v.nc"
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
#define OVARNAMEGAMMA_U "GAMMA_U"
#define OVARNAMEGAMMA_V "GAMMA_V"

/* If switched on, the diagnostics routines flush at every write_step */
#define DIAG_FLUSH

#endif
