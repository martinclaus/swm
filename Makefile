FC:=gfortran
SED:=/bin/sed
O := -O3 -fopenmp
DEBUG = #-Wall #-g
FFLAGS = -cpp -ffree-line-length-none $(defSelfCheck) $(DEBUG) $(O)
libnc = -L$(HOME)/local/netcdf-3.6.3/lib -lnetcdf
includenc = -I$(HOME)/local/netcdf-3.6.3/include
libud = -L$(HOME)/local/udunits-1.12.11/lib -ludunits
includeud = -I$(HOME)/local/udunits-1.12.11/include
udunitsdef = -D'UNITSPATH="$(HOME)/local/udunits-1.12.11/etc/udunits.dat"' -D'UD_POINTER=INTEGER(ptrKind)'
include_dir = -Iinclude
include_h = -Iinclude -Ibuild
build_dir = build/
selfcheck_dir = selfcheck/
DOXYGEN=/usr/bin/doxygen

VPATH = include:src:build

# conditional modules
#cl_elsolv := ElSolv_SOR # comment this line out if you don't want to use an elliptic solver in the calc_lib module
ifneq ($(strip $(cl_elsolv)),)
  $(info Elliptic solver used by calc_lib module: $(strip $(cl_elsolv)))
  cl_elsolv.o := $(cl_elsolv:%=%.o)
  defClElSolv := $(cl_elsolv:%=-D'CALC_LIB_ELLIPTIC_SOLVER=%')
  defClElSolv += $(cl_elsolv:%=-D'CALC_LIB_ELLIPTIC_SOLVER_HEADER="%.h"')
else
  $(info Elliptic solver used by calc_lib module: none)
endif

modules = vars_module calendar_module diag_module swm_module tracer_module io_module calc_lib dynFromFile_module $(cl_elsolv) memchunk_module swm_forcing_module swm_timestep_module swm_damping_module swm_lateralmixing_module
objects = $(modules:%=$(build_dir)%.o)

.PHONY: all clean_all clean clean-doc selfcheck doc

all     : model
#clean

model   : $(modules:%=%.o) model.o
	$(FC) $(FFLAGS) -o$(build_dir)$@ $(objects) $(build_dir)model.o $(libnc) $(libud)$

model.o : model.f90 diag_module.o vars_module.o tracer_module.o swm_module.o calc_lib.o dynFromFile_module.o model.h io.h
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

vars_module.o : vars_module.f90 io.h
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

calendar_module.o : calendar_module.f90 calendar.h include/udunits.inc
	$(FC) $(FFLAGS) $(udunitsdef) -c $< $(include_dir) $(include_h) -J$(build_dir) -o$(build_dir)$@

include/udunits.inc : $(subst -I,,$(includeud))/udunits.inc include
	@echo Convert F77 udunits interface to F90
	@$(SED) -e 's/^C/!/' -e '/#define UD_POINTER/d' $< > $@

swm_module.o : swm_module.f90 vars_module.o io_module.o swm_forcing_module.o swm_timestep_module.o swm_damping_module.o
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

diag_module.o : diag_module.f90 vars_module.o io_module.o calc_lib.o tracer_module.o swm_module.o model.h
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

io_module.o : io_module.f90 io.h vars_module.o calendar_module.o
	$(FC) $(FFLAGS) -c $< $(includenc) $(include_h) -J$(build_dir) -o$(build_dir)$@

ElSolv_SOR.o : ElSolv_SOR.f90 ElSolv_SOR.h vars_module.o model.h
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

tracer_module.o : tracer_module.f90 tracer_module.h vars_module.o io_module.o calc_lib.o model.h
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

calc_lib.o : calc_lib.f90 calc_lib.h vars_module.o model.h $(cl_elsolv.o)
	$(FC) $(FFLAGS) $(defClElSolv) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

dynFromFile_module.o : dynFromFile_module.f90 vars_module.o calc_lib.o memchunk_module.o
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

memchunk_module.o : memchunk_module.f90 io_module.o calc_lib.o vars_module.o
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

swm_forcing_module.o : swm_forcing_module.f90 model.h swm_module.h io.h vars_module.o memchunk_module.o
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

swm_timestep_module.o : swm_timestep_module.f90 model.h swm_module.h io.h vars_module.o swm_damping_module.o swm_forcing_module.o swm_lateralmixing_module.o memchunk_module.o calc_lib.o
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

swm_damping_module.o : swm_damping_module.f90 model.h swm_module.h vars_module.o
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

swm_lateralmixing_module.o : swm_lateralmixing_module.f90 model.h vars_module.o
	$(FC) $(FFLAGS) -c $< $(include_h) -J$(build_dir) -o$(build_dir)$@

selfcheck : model
	sh $(selfcheck_dir)/runselfcheck.sh

doc : doc/Doxyfile doc/html doc/latex
	@cd doc && $(DOXYGEN) $(<F)
	@cd doc/latex && make


# create folder if needed

include/  :
	mkdir -p $@

doc/html :
	mkdir -p $@

doc/latex :
	mkdir -p $@

clean_all : clean clean-doc

clean :
	@rm -fv $(build_dir)*

clean-doc :
	@rm -rvf doc/latex doc/html
