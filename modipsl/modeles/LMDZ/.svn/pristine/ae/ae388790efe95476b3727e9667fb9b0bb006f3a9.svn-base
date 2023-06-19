# This is a makefile for GNU make.
# This makefile builds "max_diff_nc".

override FFLAGS += $(addprefix -I, ${NETCDF_INC_DIR} ../NetCDF95 ../Jumble)

override LDLIBS := $(addprefix -L, ../NetCDF95 ../NR_util ../Jumble) -ljumble -lnr_util -lnetcdf95 ${LDLIBS}

name = max_diff_nc

%: %.f90
	$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

.PHONY: all clean
all: ${name}

clean:
	-rm ${name}

# Some compilers leave an object file if linking fails, so we must be
# able to process this object:
LINK.o = $(FC) $(LDFLAGS) $(TARGET_ARCH)
