FC := gfortran
O := -O3 -fopenmp
DEBUG = #-Wall -g
FFLAGS = -cpp -ffree-line-length-none $(defSelfCheck) $(DEBUG) $O
libnc = -L/bpeiler/local/netcdf-3.6.3/lib -lnetcdf
includenc = -I/bpeiler/local/netcdf-3.6.3/include

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

modules = vars_module diag_module swm_module tracer_module io_module calc_lib dynFromFile_module $(cl_elsolv) memchunk_module

.PHONY: all clean_all clean selfcheck defineSelfcheck

all     : model clean
#clean

model   : $(modules:%=%.o) model.o
	$(FC) $(FFLAGS) -o model $^ $(libnc)

model.o : model.f90 diag_module.o vars_module.o tracer_module.o swm_module.o model.h io.h
	$(FC) $(FFLAGS) -c $<

vars_module.o : vars_module.f90 io.h
	$(FC) $(FFLAGS) -c $<

swm_module.o : swm_module.f90 vars_module.o io_module.o model.h io.h
	$(FC) $(FFLAGS) -c $<

diag_module.o : diag_module.f90 vars_module.o io_module.o calc_lib.o tracer_module.o swm_module.o model.h
	$(FC) $(FFLAGS) -c $<

io_module.o : io_module.f90 io.h vars_module.o
	$(FC) $(FFLAGS) -c $< $(includenc)

ElSolv_SOR.o : ElSolv_SOR.f90 ElSolv_SOR.h vars_module.o model.h
	$(FC) $(FFLAGS) -c $<

tracer_module.o : tracer_module.f90 tracer_module.h vars_module.o io_module.o calc_lib.o model.h
	$(FC) $(FFLAGS) -c $<

calc_lib.o : calc_lib.f90 calc_lib.h vars_module.o model.h $(cl_elsolv.o)
	$(FC) $(FFLAGS) $(defClElSolv) -c $<

dynFromFile_module.o : dynFromFile_module.f90 io_module.o vars_module.o calc_lib.o memchunk_module.o
	$(FC) $(FFLAGS) -c $<

memchunk_module.o : memchunk_module.f90 io_module.o calc_lib.o
	$(FC) $(FFLAGS) -c $<

selfcheck : defineSelfcheck model

defineSelfcheck :
	$(eval defSelfCheck := -D'ISSELFCHECK')

clean_all : clean
	@rm -fv model
clean :
	@rm -fv model.o $(modules) $(modules:%=%.o) $(shell echo $(modules:%=%.mod) | tr A-Z a-z)
