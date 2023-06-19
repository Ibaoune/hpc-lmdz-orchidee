# This is a makefile for GNU make.
# This makefile builds "libjumble.a".

# 1. Source files

VPATH = Numerical

sources := $(sort avg_mag.f90 count_lines.f90 opt_merge.f90 point.f90 compare.f90 csvread.f90 new_unit.f90 read_column.f90 jumble.f90 averge.f90 get_divisors.f90 dtridgl.f90 quadrat.f90 spherical.f90 prep_file.f90 prt_cmp.f90)

# 2. Objects and library

objects := $(sources:.f90=.o)
lib = libjumble.a

# 3. Compiler-dependent part

override FFLAGS += -I../NR_util

# 4. Rules

# Extend known suffixes:
%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $<

.PHONY: all clean depend

all: ${lib}

${lib}: ${lib}(${objects})

depend depend.mk:
	makedepf90 -Wmissing -Wconfused -I${VPATH} -nosrc -u nr_util ${sources} >depend.mk

clean:
	rm -f ${lib} ${objects}

# Dependencies between object files and include files:
include depend.mk
