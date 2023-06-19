# Automatic Make rule for orchoasis

PPSRCDIR0__orchoasis = /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/build/ppsrc/orchoasis

SRCDIR0__orchoasis = /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/src_oasisdriver

FFLAGS__orchoasis__orchoasis_tools.flags: \
          FFLAGS__orchoasis.flags
	touch $(FCM_FLAGSDIR)/$@

orchoasis_tools.done: \
          orchoasis_tools.o \
          ioipsl_para.done \
          mod_orchidee_para.done
	touch $(FCM_DONEDIR)/$@

orchoasis_tools.o: \
          $(PPSRCDIR0__orchoasis)/orchoasis_tools.f90 \
          FFLAGS__orchoasis__orchoasis_tools.flags \
          ioipsl_para.o \
          mod_orchidee_para.o
	fcm_internal compile:F orchoasis $< $@

FFLAGS__orchoasis__driver2oasis.flags: \
          FFLAGS__orchoasis.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__orchoasis__driver2oasis.flags: \
          LDFLAGS__orchoasis.flags
	touch $(FCM_FLAGSDIR)/$@

LD__orchoasis__driver2oasis.flags: \
          LD__orchoasis.flags
	touch $(FCM_FLAGSDIR)/$@

driver2oasis.exe: \
          driver2oasis.o \
          LD__orchoasis__driver2oasis.flags \
          LDFLAGS__orchoasis__driver2oasis.flags \
          $(OBJECTS__ext_src) \
          $(OBJECTS__stomate) \
          $(OBJECTS__global) \
          $(OBJECTS__parameters) \
          $(OBJECTS__sechiba) \
          $(OBJECTS__orchoasis) \
          $(OBJECTS__orchidee_ol) \
          $(OBJECTS__parallel) \
          constantes.done \
          forcing_tools.done \
          globgrd.done \
          grid.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          timer.done
	fcm_internal load orchoasis $< $@

driver2oasis.o: \
          $(PPSRCDIR0__orchoasis)/driver2oasis.f90 \
          FFLAGS__orchoasis__driver2oasis.flags \
          constantes.o \
          forcing_tools.o \
          globgrd.o \
          grid.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          timer.o
	fcm_internal compile:F orchoasis $< $@

FFLAGS__orchoasis__orchideeoasis.flags: \
          FFLAGS__orchoasis.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__orchoasis__orchideeoasis.flags: \
          LDFLAGS__orchoasis.flags
	touch $(FCM_FLAGSDIR)/$@

LD__orchoasis__orchideeoasis.flags: \
          LD__orchoasis.flags
	touch $(FCM_FLAGSDIR)/$@

orchideeoasis.exe: \
          orchideeoasis.o \
          LD__orchoasis__orchideeoasis.flags \
          LDFLAGS__orchoasis__orchideeoasis.flags \
          $(OBJECTS__ext_src) \
          $(OBJECTS__stomate) \
          $(OBJECTS__global) \
          $(OBJECTS__parameters) \
          $(OBJECTS__sechiba) \
          $(OBJECTS__orchoasis) \
          $(OBJECTS__orchidee_ol) \
          $(OBJECTS__parallel) \
          control.done \
          globgrd.done \
          grid.done \
          ioipsl_para.done \
          ioipslctrl.done \
          mod_orchidee_para.done \
          orchoasis_tools.done \
          sechiba.done \
          timer.done
	fcm_internal load orchoasis $< $@

orchideeoasis.o: \
          $(PPSRCDIR0__orchoasis)/orchideeoasis.f90 \
          FFLAGS__orchoasis__orchideeoasis.flags \
          control.o \
          globgrd.o \
          grid.o \
          ioipsl_para.o \
          ioipslctrl.o \
          mod_orchidee_para.o \
          orchoasis_tools.o \
          sechiba.o \
          timer.o
	fcm_internal compile:F orchoasis $< $@

