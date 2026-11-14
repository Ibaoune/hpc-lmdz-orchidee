#- $Id: AA_make.ldef 12 2010-11-05 15:42:13Z mmaipsl $
#---------------------------------------------------------------------
#-
#-
#- $Id: AA_make.gdef 4976 2020-02-25 10:06:58Z aclsce $
#-
#- Validate the correlation between the target and the environment
#-
UTIL_DIR = ../../util
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

#- $Id: AA_make 1249 2013-04-04 12:00:37Z josefine.ghattas $
all_lib : libioipsl libparallel libparameters liborglob libstomate libsechiba

libioipsl:
	(cd ../IOIPSL/src; $(M_K) -f Makefile)

libparameters :
	(cd src_parameters ; $(M_K) -f Makefile)

libparallel :
	(cd ../../modeles/ORCHIDEE/src_parallel ; $(M_K) -f Makefile)

liborglob :
	(cd src_global ; $(M_K) -f Makefile)

libstomate :
	(cd src_stomate ; $(M_K) -f Makefile)

libsechiba :
	(cd src_sechiba ; $(M_K) -f Makefile)

driver : all_lib
	(cd src_driver ; $(M_K) -f Makefile)

config : 
	(cd src_parameters; $(M_K) -f Makefile config)
	(cd src_sechiba; $(M_K) -f Makefile config)
	(cd src_stomate; $(M_K) -f Makefile config)

clean : 
	(cd src_parameters; $(M_K) -f Makefile clean)
	(cd src_parallel; $(M_K) -f Makefile clean)
	(cd src_global; $(M_K) -f Makefile clean)
	(cd src_sechiba; $(M_K) -f Makefile clean)
	(cd src_stomate; $(M_K) -f Makefile clean)
	(cd src_driver; $(M_K) -f Makefile clean)
	(rm -f ../../lib/*)

doc:
	doxygen Doxyfile_ORCHIDEE
	export TEXINPUTS="${TEXINPUTS}:${PWD}/DOC//"; export BIBINPUTS="${BIBINPUTS}:${PWD}/DOC/BIB//"; cd docs/latex ; latex refman.tex 

bib:
	export TEXINPUTS="${TEXINPUTS}:${PWD}/DOC//"; export BIBINPUTS="${BIBINPUTS}:${PWD}/DOC/BIB//"; cd docs/latex ; bibtex refman 

index:
	export TEXINPUTS="${TEXINPUTS}:${PWD}/DOC//"; export BIBINPUTS="${BIBINPUTS}:${PWD}/DOC/BIB//"; cd docs/latex ; makeindex refman 

toc:
	export TEXINPUTS="${TEXINPUTS}:${PWD}/DOC//"; export BIBINPUTS="${BIBINPUTS}:${PWD}/DOC/BIB//"; cd docs/latex ; latex refman.tex 

