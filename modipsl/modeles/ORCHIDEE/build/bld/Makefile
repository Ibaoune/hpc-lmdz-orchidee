# Automatic Makefile

FCM_BLD_CFG = /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE
FCM_BLD_TARGETS = liborglob.a libparallel.a libsechiba.a libstomate.a libparameters.a liborchidee.a liborchoasis.a

export PERL5LIB := /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/tools/FCM_V1.2/lib

export FCM_ROOTDIR = /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/build
export FCM_BINDIR = $(FCM_ROOTDIR)/bin
export FCM_BLDDIR = $(FCM_ROOTDIR)/bld
export FCM_CACHEDIR = $(FCM_ROOTDIR)/.cache
export FCM_CFGDIR = $(FCM_ROOTDIR)/cfg
export FCM_DONEDIR = $(FCM_ROOTDIR)/done
export FCM_ETCDIR = $(FCM_ROOTDIR)/etc
export FCM_FLAGSDIR = $(FCM_ROOTDIR)/flags
export FCM_INCDIR = $(FCM_ROOTDIR)/inc
export FCM_LIBDIR = $(FCM_ROOTDIR)/lib
export FCM_OBJDIR = $(FCM_ROOTDIR)/obj
export FCM_PPSRCDIR = $(FCM_ROOTDIR)/ppsrc
export FCM_SRCDIR = $(FCM_ROOTDIR)/src
export FCM_TMPDIR = $(FCM_ROOTDIR)/tmp

export FCM_ROOTPATH = $(FCM_ROOTDIR)
export FCM_BINPATH = $(FCM_BINDIR)
export FCM_BLDPATH = $(FCM_BLDDIR)
export FCM_CACHEPATH = $(FCM_CACHEDIR)
export FCM_CFGPATH = $(FCM_CFGDIR)
export FCM_DONEPATH = $(FCM_DONEDIR)
export FCM_ETCPATH = $(FCM_ETCDIR)
export FCM_FLAGSPATH = $(FCM_FLAGSDIR)
export FCM_INCPATH = $(FCM_INCDIR)
export FCM_LIBPATH = $(FCM_LIBDIR)
export FCM_OBJPATH = $(FCM_OBJDIR)
export FCM_PPSRCPATH = $(FCM_PPSRCDIR)
export FCM_SRCPATH = $(FCM_SRCDIR)
export FCM_TMPPATH = $(FCM_TMPDIR)

export AR = ar
export ARFLAGS = rs
export CC = cc
export CC_COMPILE = -c
export CC_DEFINE = -D
export CC_INCLUDE = -I
export CC_OUTPUT = -o
export CFLAGS = 
export CPPKEYS = 
export FC = /srv/software/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/bin/mpif90
export FC_COMPILE = -c
export FC_DEFINE = -D
export FC_INCLUDE = -I
export FC_OUTPUT = -o
export FFLAGS = -cpp -ffree-line-length-0 -fdefault-real-8 -DNC_DOUBLE -O3 -fcray-pointer -fopenmp -fcray-pointer -I/srv/software/easybuild/software/netCDF-Fortran/4.4.4-foss-2018b/include -I/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/../../lib -I/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/../IOIPSL/inc -I/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/../XIOS/inc
export FPPKEYS = XIOS CPP_OMP CPP_PARA 
export LD = /srv/software/easybuild/software/OpenMPI/3.1.1-GCC-7.3.0-2.30/bin/mpif90
export LDFLAGS =  -fopenmp  -O3 -L/srv/software/easybuild/software/netCDF-Fortran/4.4.4-foss-2018b/lib -lnetcdff -lnetcdf  -L/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/../../lib -L/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/../IOIPSL/lib -lioipsl -L/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/../XIOS/lib -lxios
export LD_LIBLINK = -l
export LD_LIBSEARCH = -L
export LD_OUTPUT = -o
export SHELL = /bin/bash

export FCM_VERBOSE ?= 1

vpath % $(FCM_BINPATH)
vpath %.bash $(FCM_BINPATH)
vpath %.csh $(FCM_BINPATH)
vpath %.ksh $(FCM_BINPATH)
vpath %.pl $(FCM_BINPATH)
vpath %.pm $(FCM_BINPATH)
vpath %.pro $(FCM_BINPATH)
vpath %.py $(FCM_BINPATH)
vpath %.sh $(FCM_BINPATH)
vpath %.tcl $(FCM_BINPATH)
vpath %.exe $(FCM_BINPATH)
vpath %.mk $(FCM_BLDPATH)
vpath %.done $(FCM_DONEPATH)
vpath %.etc $(FCM_DONEPATH)
vpath %.idone $(FCM_DONEPATH)
vpath %.flags $(FCM_FLAGSPATH)
vpath %.h $(FCM_INCPATH)
vpath %.h90 $(FCM_INCPATH)
vpath %.inc $(FCM_INCPATH)
vpath %.interface $(FCM_INCPATH)
vpath %.a $(FCM_LIBPATH)
vpath %.o $(FCM_OBJPATH)
vpath %.dummy $(FCM_DONEDIR)

.PHONY : all clean

all : $(FCM_BLD_TARGETS)

clean : 
	rm -rf $(FCM_BINDIR) $(FCM_LIBDIR) $(FCM_OBJDIR) $(FCM_DONEDIR)
	rm -f lib__*.a *.o

$(FCM_DONEDIR)/FCM_CP.dummy :
	touch $@

OBJECTS__stomate = lpj_constraints.o lpj_cover.o lpj_crown.o lpj_establish.o lpj_fire.o lpj_gap.o lpj_kill.o lpj_light.o lpj_pftinout.o stomate.o stomate_alloc.o stomate_data.o stomate_io.o stomate_lcchange.o stomate_litter.o stomate_lpj.o stomate_npp.o stomate_phenology.o stomate_prescribe.o stomate_resp.o stomate_season.o stomate_soilcarbon.o stomate_turnover.o stomate_vmax.o stomate_woodharvest.o

libstomate.a : $(OBJECTS__stomate)
	fcm_internal archive $@ $(^F)

OBJECTS__sechiba = chemistry.o condveg.o diffuco.o enerbil.o explicitsnow.o hydrol.o intersurf.o ioipslctrl.o qsat_moisture.o routing.o routing_simple.o routing_wrapper.o sechiba.o sechiba_io_p.o slowproc.o thermosoil.o

libsechiba.a : $(OBJECTS__sechiba)
	fcm_internal archive $@ $(^F)

OBJECTS__parameters = constantes.o constantes_mtc.o constantes_soil.o constantes_soil_var.o constantes_var.o control.o pft_parameters.o pft_parameters_var.o vertical_soil.o vertical_soil_var.o

libparameters.a : $(OBJECTS__parameters)
	fcm_internal archive $@ $(^F)

OBJECTS__parallel = ioipsl_para.o mod_orchidee_mpi_data.o mod_orchidee_mpi_transfert.o mod_orchidee_omp_data.o mod_orchidee_omp_transfert.o mod_orchidee_para.o mod_orchidee_para_var.o mod_orchidee_transfert_para.o orch_write_field.o orch_write_field_p.o timer.o xios_orchidee.o

libparallel.a : $(OBJECTS__parallel)
	fcm_internal archive $@ $(^F)

OBJECTS__orchoasis = orchoasis_tools.o

liborchoasis.a : $(OBJECTS__orchoasis)
	fcm_internal archive $@ $(^F)

OBJECTS__orchidee_ol = forcing_tools.o forcingdaily_tools.o getlandseamask.o globgrd.o readdim2.o weather.o

liborchidee_ol.a : $(OBJECTS__orchidee_ol)
	fcm_internal archive $@ $(^F)

OBJECTS__global = matrix_resolution.o grid.o grid_var.o haversine.o interpol_help.o interpweight.o interregxy.o module_llxy.o polygones.o solar.o time.o

liborglob.a : $(OBJECTS__global)
	fcm_internal archive $@ $(^F)

OBJECTS__ext_src =

libext_src.a : $(OBJECTS__ext_src)
	fcm_internal archive $@ $(^F)

export OBJECTS = $(OBJECTS__ext_src) $(OBJECTS__global) $(OBJECTS__orchidee_ol) $(OBJECTS__orchoasis) $(OBJECTS__parallel) $(OBJECTS__parameters) $(OBJECTS__sechiba) $(OBJECTS__stomate)

liborchidee.a : $(OBJECTS)
	fcm_internal archive $@ $(^F)

CC.flags :
	touch $(FCM_FLAGSDIR)/$@

CFLAGS.flags : CC.flags
	touch $(FCM_FLAGSDIR)/$@

CFLAGS__ext_src.flags CFLAGS__global.flags CFLAGS__orchidee_ol.flags CFLAGS__orchoasis.flags CFLAGS__parallel.flags CFLAGS__parameters.flags CFLAGS__sechiba.flags CFLAGS__stomate.flags : CFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS.flags :
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS__ext_src.flags CPPKEYS__global.flags CPPKEYS__orchidee_ol.flags CPPKEYS__orchoasis.flags CPPKEYS__parallel.flags CPPKEYS__parameters.flags CPPKEYS__sechiba.flags CPPKEYS__stomate.flags : CPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

FC.flags :
	touch $(FCM_FLAGSDIR)/$@

FFLAGS.flags : FC.flags
	touch $(FCM_FLAGSDIR)/$@

FFLAGS__ext_src.flags FFLAGS__global.flags FFLAGS__orchidee_ol.flags FFLAGS__orchoasis.flags FFLAGS__parallel.flags FFLAGS__parameters.flags FFLAGS__sechiba.flags FFLAGS__stomate.flags : FFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS.flags :
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__ext_src.flags FPPKEYS__global.flags FPPKEYS__orchidee_ol.flags FPPKEYS__orchoasis.flags FPPKEYS__parallel.flags FPPKEYS__parameters.flags FPPKEYS__sechiba.flags FPPKEYS__stomate.flags : FPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

LD.flags :
	touch $(FCM_FLAGSDIR)/$@

LD__ext_src.flags LD__global.flags LD__orchidee_ol.flags LD__orchoasis.flags LD__parallel.flags LD__parameters.flags LD__sechiba.flags LD__stomate.flags : LD.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS.flags :
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__ext_src.flags LDFLAGS__global.flags LDFLAGS__orchidee_ol.flags LDFLAGS__orchoasis.flags LDFLAGS__parallel.flags LDFLAGS__parameters.flags LDFLAGS__sechiba.flags LDFLAGS__stomate.flags : LDFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

include $(FCM_BLDDIR)/global.mk
include $(FCM_BLDDIR)/orchidee_ol.mk
include $(FCM_BLDDIR)/orchoasis.mk
include $(FCM_BLDDIR)/parallel.mk
include $(FCM_BLDDIR)/parameters.mk
include $(FCM_BLDDIR)/sechiba.mk
include $(FCM_BLDDIR)/stomate.mk

# EOF
