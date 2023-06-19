# This is a makefile for GNU make.
# This makefile builds "NetCDF95".

# 1. Source files

sources = nf95_inq_dimid.f90 nf95_inquire_dimension.f90 nf95_inq_varid.f90 nf95_inquire_variable.f90 nf95_create.f90 nf95_def_dim.f90 nf95_redef.f90 nf95_enddef.f90 nf95_close.f90 nf95_copy_att.f90 nf95_inquire_attribute.f90 nf95_inquire.f90 netcdf95.f90 nf95_def_var.f90 nf95_gw_var.f90 nf95_put_var.f90 nf95_put_att.f90 handle_err.f90 nf95_get_att.f90 nf95_get_var.f90 find_coord.f90 nf95_open.f90

# 2. Objects and library

objects := $(sources:.f90=.o)
lib = libnetcdf95.a

# 3. Compiler-dependent part

override FFLAGS += $(if ${NETCDF_INC_DIR}, -I${NETCDF_INC_DIR})

# 4. Rules

# Extend known suffixes:
%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $<

.PHONY: all clean depend

all: ${lib}
${lib}: ${lib}(${objects})

depend depend.mk:
	makedepf90 -Wmissing -Wconfused -nosrc -u netcdf -u typesizes ${sources} >depend.mk

clean:
	rm -f ${lib} ${objects}

# Dependencies between object files and include files:
include depend.mk
