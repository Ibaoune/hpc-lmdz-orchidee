#-
#- $Id: AA_make.ldef 1311 2013-06-04 08:13:54Z josefine.ghattas $
#-
#---------------------------------------------------------------------
#- Creation des elements relatifs a PARAMETER
#---------------------------------------------------------------------
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
ORCHIDEE_LIB = $(LIBDIR)/liborchidee.a
MODEL_LIB = $(LIBDIR)/libparameters.a
SXMODEL_LIB = $(MODEL_LIB)
#-
#- $Id: AA_make.gdef 4976 2020-02-25 10:06:58Z aclsce $
#-
#- Validate the correlation between the target and the environment
#-
UTIL_DIR = ../../../util
#-
#-Q-
######-Q- ada      F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- curie  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- cur_mono  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
####-Q- lxiv8    F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- occigen  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
#-
#- Global definitions for Toubkal at UM6P
LIB_MPI = MPI1
LIB_MPI_BIS = MPI1
PRISM_ARCH = X64
PRISM_NAME = toubkal
FCM_ARCH = X64_TOUBKAL
M_K = gmake
P_C = cpp
P_O = -P -C $(P_P)
F_C = mpif90 -c -cpp
F_D =
F_P = -fdefault-real-8
F_O = -DCPP_PARA -O3 $(F_D) $(F_P) -I$(MODDIR) -J$(MODDIR)
F_L = mpif90
M_M = 0
L_X = 0
L_O =
A_C = ar -rs
A_G = ar -x
C_C = mpicc -c
C_O =
C_L = mpicc
#-
NCDF_INC = /srv/software/easybuild/software/netCDF-Fortran/4.4.4-foss-2018b/include -I/srv/software/easybuild/software/netCDF/4.6.1-foss-2018b/include
NCDF_LIB = -L/srv/software/easybuild/software/netCDF-Fortran/4.4.4-foss-2018b/lib -lnetcdff -L/srv/software/easybuild/software/netCDF/4.6.1-foss-2018b/lib -lnetcdf
#-
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
#- $Id: AA_make 2582 2015-03-06 11:40:41Z josefine.ghattas $
#-
IOIPSL_LIB = $(LIBDIR)/libioipsl.a
SXIOIPSL_LIB = $(IOIPSL_LIB)
#-
PARALLEL_LIB = $(LIBDIR)/libparallel.a
SXPARALLEL_LIB = $(PARALLEL_LIB)
#-
MODS1 = constantes.f90 \
	constantes_var.f90   \
	constantes_soil.f90   \
	constantes_soil_var.f90   \
	constantes_mtc.f90   \
	control.f90        \
	pft_parameters.f90 \
	pft_parameters_var.f90 
#
OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)
#-
all:
	$(M_K) libioipsl
	$(M_K) libparallel
	$(M_K) m_all
	@echo parameter is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libioipsl:
	(cd ../../IOIPSL/src; $(M_K) -f Makefile)

libparallel:
	(cd ../src_parallel; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	$(A_C) $(ORCHIDEE_LIB) $*.o
	$(RM) $*.o

config:
	$(BINDIR)/Fparser -name PARAMETERS  $(MODS1)
	echo 'Configuration of PARAMETERS done'
 
clean:
	$(RM) $(MODEL_LIB)
 
$(MODEL_LIB)(pft_parameters.o): \
  $(MODEL_LIB)(pft_parameters_var.o) \
  $(MODEL_LIB)(constantes_mtc.o)
$(MODEL_LIB)(constantes_soil.o): \
  $(MODEL_LIB)(constantes_soil_var.o) \
  $(MODEL_LIB)(constantes.o)
$(MODEL_LIB)(control.o): \
  $(MODEL_LIB)(constantes_soil.o) \
  $(MODEL_LIB)(constantes_var.o)  \
  $(MODEL_LIB)(pft_parameters.o)

$(MODEL_LIB)(constantes.o): \
  $(MODEL_LIB)(constantes_var.o) \
  $(MODEL_LIB)(constantes_soil_var.o) \
  $(PARALLEL_LIB) \
  $(IOIPSL_LIB)

