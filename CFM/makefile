SHELL = /bin/bash
MAIN = kalman

FC       = ifort #-heap-arrays
#FC	 = gfortran
FFLAGS   = -mcmodel=medium -shared-intel -traceback
LDFLAGS  = $(FFLAGS)

OBJECTS  = kalmfilt.f90    build_global_data.f90 \
output_data.f90  mtxchi.f90 build_result_fld.f90 read_config.f90 \
global_data.f90 estimation_proc.f90 count_err.f90


$(MAIN): $(OBJECTS)
	$(FC) *.f90 -o $(MAIN) $(LDFLAGS)
$(OBJECTS): $(INCF)

