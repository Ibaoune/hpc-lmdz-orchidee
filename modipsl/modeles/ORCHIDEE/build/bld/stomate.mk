# Automatic Make rule for stomate

SRCDIR0__stomate = /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/src_stomate

PPSRCDIR0__stomate = /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/Model/LMDZOR_v6.2.2/modipsl/modeles/ORCHIDEE/build/ppsrc/stomate

stomate.etc : \
          $(SRCDIR0__stomate)/AA_make \
          $(SRCDIR0__stomate)/AA_make.ldef \
          $(FCM_DONEDIR)/FCM_CP.dummy
	cp $^ $(FCM_ETCDIR)
	touch $(FCM_DONEDIR)/$@

FFLAGS__stomate__stomate_npp.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_npp.done: \
          stomate_npp.o \
          constantes.done \
          constantes_soil.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

stomate_npp.o: \
          $(PPSRCDIR0__stomate)/stomate_npp.f90 \
          FFLAGS__stomate__stomate_npp.flags \
          constantes.o \
          constantes_soil.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__lpj_light.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

lpj_light.done: \
          lpj_light.o \
          constantes.done \
          ioipsl_para.done \
          stomate_data.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

lpj_light.o: \
          $(PPSRCDIR0__stomate)/lpj_light.f90 \
          FFLAGS__stomate__lpj_light.flags \
          constantes.o \
          ioipsl_para.o \
          stomate_data.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_io.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_io.done: \
          stomate_io.o \
          constantes.done \
          constantes_soil.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

stomate_io.o: \
          $(PPSRCDIR0__stomate)/stomate_io.f90 \
          FFLAGS__stomate__stomate_io.flags \
          constantes.o \
          constantes_soil.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__lpj_establish.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

lpj_establish.done: \
          lpj_establish.o \
          constantes.done \
          grid.done \
          ioipsl_para.done \
          stomate_data.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

lpj_establish.o: \
          $(PPSRCDIR0__stomate)/lpj_establish.f90 \
          FFLAGS__stomate__lpj_establish.flags \
          constantes.o \
          grid.o \
          ioipsl_para.o \
          stomate_data.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_lcchange.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_lcchange.done: \
          stomate_lcchange.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

stomate_lcchange.o: \
          $(PPSRCDIR0__stomate)/stomate_lcchange.f90 \
          FFLAGS__stomate__stomate_lcchange.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_data.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_data.done: \
          stomate_data.o \
          constantes.done \
          pft_parameters.done \
          time.done
	touch $(FCM_DONEDIR)/$@

stomate_data.o: \
          $(PPSRCDIR0__stomate)/stomate_data.f90 \
          FFLAGS__stomate__stomate_data.flags \
          constantes.o \
          pft_parameters.o \
          time.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_season.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_season.done: \
          stomate_season.o \
          constantes.done \
          constantes_soil.done \
          grid.done \
          ioipsl_para.done \
          pft_parameters.done \
          solar.done \
          stomate_data.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

stomate_season.o: \
          $(PPSRCDIR0__stomate)/stomate_season.f90 \
          FFLAGS__stomate__stomate_season.flags \
          constantes.o \
          constantes_soil.o \
          grid.o \
          ioipsl_para.o \
          pft_parameters.o \
          solar.o \
          stomate_data.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_vmax.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_vmax.done: \
          stomate_vmax.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

stomate_vmax.o: \
          $(PPSRCDIR0__stomate)/stomate_vmax.f90 \
          FFLAGS__stomate__stomate_vmax.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate.done: \
          stomate.o \
          constantes.done \
          constantes_soil.done \
          grid.done \
          ioipsl_para.done \
          matrix_resolution.done \
          mod_orchidee_para.done \
          pft_parameters.done \
          stomate_data.done \
          stomate_io.done \
          stomate_litter.done \
          stomate_lpj.done \
          stomate_resp.done \
          stomate_season.done \
          stomate_soilcarbon.done \
          stomate_vmax.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

stomate.o: \
          $(PPSRCDIR0__stomate)/stomate.f90 \
          FFLAGS__stomate__stomate.flags \
          constantes.o \
          constantes_soil.o \
          grid.o \
          ioipsl_para.o \
          matrix_resolution.o \
          mod_orchidee_para.o \
          pft_parameters.o \
          stomate_data.o \
          stomate_io.o \
          stomate_litter.o \
          stomate_lpj.o \
          stomate_resp.o \
          stomate_season.o \
          stomate_soilcarbon.o \
          stomate_vmax.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_woodharvest.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_woodharvest.done: \
          stomate_woodharvest.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

stomate_woodharvest.o: \
          $(PPSRCDIR0__stomate)/stomate_woodharvest.f90 \
          FFLAGS__stomate__stomate_woodharvest.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_turnover.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_turnover.done: \
          stomate_turnover.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

stomate_turnover.o: \
          $(PPSRCDIR0__stomate)/stomate_turnover.f90 \
          FFLAGS__stomate__stomate_turnover.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_resp.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_resp.done: \
          stomate_resp.o \
          constantes.done \
          constantes_soil.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

stomate_resp.o: \
          $(PPSRCDIR0__stomate)/stomate_resp.f90 \
          FFLAGS__stomate__stomate_resp.flags \
          constantes.o \
          constantes_soil.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_litter.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_litter.done: \
          stomate_litter.o \
          constantes.done \
          constantes_soil.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

stomate_litter.o: \
          $(PPSRCDIR0__stomate)/stomate_litter.f90 \
          FFLAGS__stomate__stomate_litter.flags \
          constantes.o \
          constantes_soil.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__lpj_gap.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

lpj_gap.done: \
          lpj_gap.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

lpj_gap.o: \
          $(PPSRCDIR0__stomate)/lpj_gap.f90 \
          FFLAGS__stomate__lpj_gap.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__lpj_crown.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

lpj_crown.done: \
          lpj_crown.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

lpj_crown.o: \
          $(PPSRCDIR0__stomate)/lpj_crown.f90 \
          FFLAGS__stomate__lpj_crown.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__lpj_cover.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

lpj_cover.done: \
          lpj_cover.o \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

lpj_cover.o: \
          $(PPSRCDIR0__stomate)/lpj_cover.f90 \
          FFLAGS__stomate__lpj_cover.flags \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_lpj.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_lpj.done: \
          stomate_lpj.o \
          constantes.done \
          constantes_soil.done \
          grid.done \
          ioipsl_para.done \
          lpj_constraints.done \
          lpj_cover.done \
          lpj_crown.done \
          lpj_establish.done \
          lpj_fire.done \
          lpj_gap.done \
          lpj_kill.done \
          lpj_light.done \
          lpj_pftinout.done \
          pft_parameters.done \
          stomate_alloc.done \
          stomate_data.done \
          stomate_lcchange.done \
          stomate_litter.done \
          stomate_npp.done \
          stomate_phenology.done \
          stomate_prescribe.done \
          stomate_soilcarbon.done \
          stomate_turnover.done \
          stomate_vmax.done \
          stomate_woodharvest.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

stomate_lpj.o: \
          $(PPSRCDIR0__stomate)/stomate_lpj.f90 \
          FFLAGS__stomate__stomate_lpj.flags \
          constantes.o \
          constantes_soil.o \
          grid.o \
          ioipsl_para.o \
          lpj_constraints.o \
          lpj_cover.o \
          lpj_crown.o \
          lpj_establish.o \
          lpj_fire.o \
          lpj_gap.o \
          lpj_kill.o \
          lpj_light.o \
          lpj_pftinout.o \
          pft_parameters.o \
          stomate_alloc.o \
          stomate_data.o \
          stomate_lcchange.o \
          stomate_litter.o \
          stomate_npp.o \
          stomate_phenology.o \
          stomate_prescribe.o \
          stomate_soilcarbon.o \
          stomate_turnover.o \
          stomate_vmax.o \
          stomate_woodharvest.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__lpj_kill.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

lpj_kill.done: \
          lpj_kill.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

lpj_kill.o: \
          $(PPSRCDIR0__stomate)/lpj_kill.f90 \
          FFLAGS__stomate__lpj_kill.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_soilcarbon.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_soilcarbon.done: \
          stomate_soilcarbon.o \
          constantes.done \
          ioipsl_para.done \
          stomate_data.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

stomate_soilcarbon.o: \
          $(PPSRCDIR0__stomate)/stomate_soilcarbon.f90 \
          FFLAGS__stomate__stomate_soilcarbon.flags \
          constantes.o \
          ioipsl_para.o \
          stomate_data.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_prescribe.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_prescribe.done: \
          stomate_prescribe.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

stomate_prescribe.o: \
          $(PPSRCDIR0__stomate)/stomate_prescribe.f90 \
          FFLAGS__stomate__stomate_prescribe.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_phenology.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_phenology.done: \
          stomate_phenology.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

stomate_phenology.o: \
          $(PPSRCDIR0__stomate)/stomate_phenology.f90 \
          FFLAGS__stomate__stomate_phenology.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__stomate_alloc.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

stomate_alloc.done: \
          stomate_alloc.o \
          constantes.done \
          constantes_soil.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

stomate_alloc.o: \
          $(PPSRCDIR0__stomate)/stomate_alloc.f90 \
          FFLAGS__stomate__stomate_alloc.flags \
          constantes.o \
          constantes_soil.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__lpj_constraints.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

lpj_constraints.done: \
          lpj_constraints.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

lpj_constraints.o: \
          $(PPSRCDIR0__stomate)/lpj_constraints.f90 \
          FFLAGS__stomate__lpj_constraints.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__lpj_pftinout.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

lpj_pftinout.done: \
          lpj_pftinout.o \
          constantes.done \
          grid.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done
	touch $(FCM_DONEDIR)/$@

lpj_pftinout.o: \
          $(PPSRCDIR0__stomate)/lpj_pftinout.f90 \
          FFLAGS__stomate__lpj_pftinout.flags \
          constantes.o \
          grid.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o
	fcm_internal compile:F stomate $< $@

FFLAGS__stomate__lpj_fire.flags: \
          FFLAGS__stomate.flags
	touch $(FCM_FLAGSDIR)/$@

lpj_fire.done: \
          lpj_fire.o \
          constantes.done \
          ioipsl_para.done \
          pft_parameters.done \
          stomate_data.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

lpj_fire.o: \
          $(PPSRCDIR0__stomate)/lpj_fire.f90 \
          FFLAGS__stomate__lpj_fire.flags \
          constantes.o \
          ioipsl_para.o \
          pft_parameters.o \
          stomate_data.o \
          xios_orchidee.o
	fcm_internal compile:F stomate $< $@

