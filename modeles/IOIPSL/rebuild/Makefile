#-
#- $Id: AA_make.ldef 53 2007-06-06 15:27:26Z bellier $
#-
#---------------------------------------------------------------------
#- Local parameters related to the "tools" of IOIPSL (rebuild, config)
#---------------------------------------------------------------------
SHELL = /bin/ksh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
P_P =
MODEL_LIB = $(LIBDIR)/libioipsl.a
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
#- $Id: AA_make 53 2007-06-06 15:27:26Z bellier $
#-
#---------------------------------------------------------------------
#- Create "rebuild" related to IOIPSL
#---------------------------------------------------------------------

SXMODEL_LIB = $(MODEL_LIB)
ifeq ($(L_X),1)
SXMODEL_LIB = $(LIBDIR)/libsxioipsl.a
endif

rbld: libioipsl $(BINDIR)/rebuild $(BINDIR)/flio_rbld
	@echo rebuild is OK

$(BINDIR)/rebuild: ./rebuild
	cp -p ./rebuild $(BINDIR)/rebuild
	chmod a+rx $(BINDIR)/rebuild

$(BINDIR)/flio_rbld: $(SXMODEL_LIB) flio_rbld.f90
	$(F_C) $(F_O) -I$(NCDF_INC) flio_rbld.f90
	$(F_L) flio_rbld.o -o $(BINDIR)/flio_rbld \
               $(L_O) $(SXMODEL_LIB) $(NCDF_LIB)
	chmod a+rx $(BINDIR)/flio_rbld
	$(RM) flio_rbld.o

libioipsl:
	(cd ../src; $(M_K) -f Makefile)

clean:
	$(RM) *.*~ Makefile~ core *.o *.mod i.*.L *.L i.*.f90 *.M
	$(RM) $(BINDIR)/rebuild $(BINDIR)/flio_rbld

