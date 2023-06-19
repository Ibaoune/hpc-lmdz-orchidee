# Automatic Make rule for global

PPSRCDIR0__global = /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/build/ppsrc/global

SRCDIR0__global = /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/src_global

global.etc : \
          $(SRCDIR0__global)/AA_make \
          $(SRCDIR0__global)/AA_make.ldef \
          $(FCM_DONEDIR)/FCM_CP.dummy
	cp $^ $(FCM_ETCDIR)
	touch $(FCM_DONEDIR)/$@

FFLAGS__global__interregxy.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

interregxy.done: \
          interregxy.o \
          constantes.done \
          grid.done \
          mod_orchidee_para.done \
          module_llxy.done \
          polygones.done
	touch $(FCM_DONEDIR)/$@

interregxy.o: \
          $(PPSRCDIR0__global)/interregxy.f90 \
          FFLAGS__global__interregxy.flags \
          constantes.o \
          grid.o \
          mod_orchidee_para.o \
          module_llxy.o \
          polygones.o
	fcm_internal compile:F global $< $@

FFLAGS__global__grid_var.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

grid_var.done: \
          grid_var.o
	touch $(FCM_DONEDIR)/$@

grid_var.o: \
          $(PPSRCDIR0__global)/grid_var.f90 \
          FFLAGS__global__grid_var.flags
	fcm_internal compile:F global $< $@

FFLAGS__global__interpol_help.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

interpol_help.done: \
          interpol_help.o \
          constantes.done \
          grid.done \
          interregxy.done \
          mod_orchidee_para.done
	touch $(FCM_DONEDIR)/$@

interpol_help.o: \
          $(PPSRCDIR0__global)/interpol_help.f90 \
          FFLAGS__global__interpol_help.flags \
          constantes.o \
          grid.o \
          interregxy.o \
          mod_orchidee_para.o
	fcm_internal compile:F global $< $@

FFLAGS__global__grid.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

grid.done: \
          grid.o \
          grid_var.done \
          haversine.done \
          mod_orchidee_para_var.done \
          mod_orchidee_transfert_para.done \
          module_llxy.done
	touch $(FCM_DONEDIR)/$@

grid.o: \
          $(PPSRCDIR0__global)/grid.f90 \
          FFLAGS__global__grid.flags \
          grid_var.o \
          haversine.o \
          mod_orchidee_para_var.o \
          mod_orchidee_transfert_para.o \
          module_llxy.o
	fcm_internal compile:F global $< $@

FFLAGS__global__polygones.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

polygones.done: \
          polygones.o
	touch $(FCM_DONEDIR)/$@

polygones.o: \
          $(PPSRCDIR0__global)/polygones.f90 \
          FFLAGS__global__polygones.flags
	fcm_internal compile:F global $< $@

FFLAGS__global__gauss_jordan_method.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

matrix_resolution.done: \
          matrix_resolution.o \
          constantes.done
	touch $(FCM_DONEDIR)/$@

matrix_resolution.o: \
          $(PPSRCDIR0__global)/gauss_jordan_method.f90 \
          FFLAGS__global__gauss_jordan_method.flags \
          constantes.o
	fcm_internal compile:F global $< $@

FFLAGS__global__interpweight.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

interpweight.done: \
          interpweight.o \
          constantes_var.done \
          grid.done \
          interpol_help.done \
          ioipsl_para.done \
          mod_orchidee_para.done
	touch $(FCM_DONEDIR)/$@

interpweight.o: \
          $(PPSRCDIR0__global)/interpweight.f90 \
          FFLAGS__global__interpweight.flags \
          constantes_var.o \
          grid.o \
          interpol_help.o \
          ioipsl_para.o \
          mod_orchidee_para.o
	fcm_internal compile:F global $< $@

FFLAGS__global__haversine.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

haversine.done: \
          haversine.o \
          constantes_var.done \
          module_llxy.done
	touch $(FCM_DONEDIR)/$@

haversine.o: \
          $(PPSRCDIR0__global)/haversine.f90 \
          FFLAGS__global__haversine.flags \
          constantes_var.o \
          module_llxy.o
	fcm_internal compile:F global $< $@

FFLAGS__global__solar.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

solar.done: \
          solar.o \
          constantes.done \
          ioipsl_para.done \
          time.done
	touch $(FCM_DONEDIR)/$@

solar.o: \
          $(PPSRCDIR0__global)/solar.f90 \
          FFLAGS__global__solar.flags \
          constantes.o \
          ioipsl_para.o \
          time.o
	fcm_internal compile:F global $< $@

FFLAGS__global__time.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

time.done: \
          time.o \
          constantes_var.done \
          ioipsl_para.done \
          mod_orchidee_para_var.done
	touch $(FCM_DONEDIR)/$@

time.o: \
          $(PPSRCDIR0__global)/time.f90 \
          FFLAGS__global__time.flags \
          constantes_var.o \
          ioipsl_para.o \
          mod_orchidee_para_var.o
	fcm_internal compile:F global $< $@

FFLAGS__global__module_llxy.flags: \
          FFLAGS__global.flags
	touch $(FCM_FLAGSDIR)/$@

module_llxy.done: \
          module_llxy.o \
          constantes_var.done
	touch $(FCM_DONEDIR)/$@

module_llxy.o: \
          $(PPSRCDIR0__global)/module_llxy.f90 \
          FFLAGS__global__module_llxy.flags \
          constantes_var.o
	fcm_internal compile:F global $< $@

