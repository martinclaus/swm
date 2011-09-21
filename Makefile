# from /data1/mclaus/model_yao , 20110909
O = -O3
CFLAGS = -fopenmp -cpp
libnc = -L/home/mclaus/local/netcdf-3.6.3/lib -lnetcdf
includenc = -I/home/mclaus/local/netcdf-3.6.3/include
modules = vars_module diag_module timestep_module

# old stuff:
#O = -O3
#CFLAGS = -fopenmp -cpp
#libnc = -L/home/mclaus/local/lib -lnetcdf
#includenc = -I/home/mclaus/local/include
#modules = vars_module diag_module timestep_module

all     : model clean

model   : $(modules:%=%.o) model.o
	gfortran $O $(CFLAGS) -o model $^ $(libnc)

model.o : model.f90 diag_module.o vars_module.o
	gfortran $O $(CFLAGS) -c $< $(includenc)

vars_module.o : vars_module.f90
	gfortran $O $(CFLAGS) -c $<

timestep_module.o : timestep_module.f90
	gfortran $O $(CFLAGS) -c $<

diag_module.o     : diag_module.f90 vars_module.o timestep_module.o
	gfortran $O $(CFLAGS) -c $< $(includenc)

clean_all : clean
	rm -fv model
clean :
	rm -fv model.o $(modules) $(modules:%=%.o) $(modules:%=%.mod)
