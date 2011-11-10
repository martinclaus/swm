CFORTAN := gfortran
O := -O3
CFLAGS := -fopenmp -cpp
libnc = -L/home/mclaus/local/netcdf-3.6.3/lib -lnetcdf
includenc = -I/home/mclaus/local/netcdf-3.6.3/include

# conditional modules
diag_elsolv := ElSolv_SOR # comment this line out if you don't want to use an elliptic solver in the diag_module
ifneq ($(strip $(diag_elsolv)),)
  $(info Elliptic solver used by diag_module: $(strip $(diag_elsolv))) 
  diag_elsolv.o := $(diag_elsolv:%=%.o)
  defDiagElSolv := $(diag_elsolv:%=-D'DIAG_ELLIPTIC_SOLVER=%')
  defDiagElSolv += $(diag_elsolv:%=-D'DIAG_ELLIPTIC_SOLVER_HEADER="%.h"')
else
  $(info Elliptic solver used by diag_module: none) 
endif

modules = vars_module diag_module timestep_module $(diag_elsolv)

.PHONY: all clean_all clean

all     : model clean

model   : $(modules:%=%.o) model.o
	$(CFORTAN) $O $(CFLAGS) -o model $^ $(libnc)

model.o : model.f90 diag_module.o vars_module.o model.h
	$(CFORTAN) $O $(CFLAGS) -c $< $(includenc)

vars_module.o : vars_module.f90 io.h
	$(CFORTAN) $O $(CFLAGS) -c $<

timestep_module.o : timestep_module.f90
	$(CFORTAN) $O $(CFLAGS) -c $<

diag_module.o     : diag_module.f90 vars_module.o diag_module.h io.h $(diag_elsolv.o)
	$(CFORTAN) $O $(CFLAGS) $(defDiagElSolv) -c $< $(includenc)

ElSolv_SOR.o : ElSolv_SOR.f90 vars_module.o ElSolv_SOR.h model.h
	$(CFORTAN) $O $(CFLAGS) -c $<

clean_all : clean
	@rm -fv model
clean :
	@rm -fv model.o $(modules) $(modules:%=%.o) $(shell echo $(modules:%=%.mod) | tr A-Z a-z)
