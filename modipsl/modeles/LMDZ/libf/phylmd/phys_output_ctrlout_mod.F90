!
! $Id: phys_output_ctrlout_mod.F90 3823 2021-02-05 09:00:13Z musat $
!
MODULE phys_output_ctrlout_mod

  USE phys_output_var_mod
  USE indice_sol_mod
  USE aero_mod

  IMPLICIT NONE
      INTEGER, PRIVATE :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition pour chaque variable du niveau d ecriture dans chaque fichier,
!! de son nom, de sa description, de son unité et du type d'écriture.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/ histmth, histday, histhf, histins /),'!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CHARACTER(len=20), dimension(nfiles) :: TEF = type_ecri_files

!!! saving lon and lat as variables for CMIP6 DataRequest
  TYPE(ctrl_out), SAVE :: o_longitude = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'io_lon', '', '', (/ ('once', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_latitude = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'io_lat', '', '', (/ ('once', i=1, 10) /))

!!! Composantes de la coordonnee sigma-hybride
!!! Ap et Bp et interfaces
  TYPE(ctrl_out), SAVE :: o_Ahyb = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'Ahyb', 'Ahyb at level interface', '', (/ ('once', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_Bhyb = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'Bhyb', 'Bhyb at level interface', '', (/ ('once', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_Ahyb_bounds = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'Ahyb_bounds', '', '', (/ ('once', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_Bhyb_bounds = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'Bhyb_bounds', '', '', (/ ('once', i=1, 10) /))
!!! Composantes de la coordonnee sigma-hybride  au milieu des couches
!!! Aps et Bps et interfaces
  TYPE(ctrl_out), SAVE :: o_Ahyb_mid = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'Ahyb_mid', 'Ahyb at the middle of the level', '', (/ ('once', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_Bhyb_mid = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'Bhyb_mid', 'Bhyb at the middle of the level', '', (/ ('once', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_Ahyb_mid_bounds = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'Ahyb_mid_bounds', '', '', (/ ('once', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_Bhyb_mid_bounds = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'Bhyb_mid_bounds', '', '', (/ ('once', i=1, 10) /))

  TYPE(ctrl_out), SAVE :: o_Alt = ctrl_out((/ 1, 1, 1, 1, 1, 1, 11, 11, 11, 11/), &
    'Alt', '', '', (/ ('', i=1, 10) /))

!!! 1D
  TYPE(ctrl_out), SAVE :: o_phis = ctrl_out((/ 1, 1, 10, 5, 1, 1, 11, 11, 11, 11/), &
    'phis', 'Surface geop.height', 'm2/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_aire = ctrl_out((/ 1, 1, 10,  10, 1, 1, 11, 11, 11, 11/), &
    'aire', 'Grid area', '-', (/ 'once', 'once', 'once', 'once', 'once', 'once', &
                                 'once', 'once', 'once', 'once' /))
  TYPE(ctrl_out), SAVE :: o_contfracATM = ctrl_out((/ 10, 1,  1, 10, 10, 10, 11, 11, 11, 11/), &
    'contfracATM', '% sfce ter+lic', '-', &
       (/ 'once', 'once', 'once', 'once', 'once', 'once', 'once', 'once', 'once', 'once' /)) 
  TYPE(ctrl_out), SAVE :: o_contfracOR = ctrl_out((/ 10, 1,  10, 10, 10, 10, 11, 11, 11, 11/), &
    'contfracOR', '% sfce terre OR', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_aireTER = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'aireTER', 'Grid area CONT', '-', (/ ('', i=1, 10) /))

!!! 2D
  TYPE(ctrl_out), SAVE :: o_sza = ctrl_out((/ 1, 1, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'sza', 'Solar zenithal angle', 'degrees', (/ ('', i=1, 10) /))

! Marine

  TYPE(ctrl_out), SAVE :: o_alt_tropo = ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'alt_tropo','Tropopause pressure','hPa',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_prop_hc = ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_prop_hc','Proportion of high clouds',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_prop_hist = &
  ctrl_out((/1,1,1,1,1,1,10,10,10,10/),&
  'map_prop_hist','Proportion of high ice semi-transp clouds',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_emis_hc = &
  ctrl_out((/1,1,1,1,1,1,10,10,10,10/),&
  'map_emis_hc','Emissivity of high clouds',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_iwp_hc = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_iwp_hc','Ice water path of high clouds','g/m2',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_deltaz_hc = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_deltaz_hc','geom thickness of high clouds','m',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_pcld_hc = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_pcld_hc','cloud pressure of high clouds','hPa',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

   TYPE(ctrl_out), SAVE :: o_map_tcld_hc = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_tcld_hc','cloud temperature of high clouds','K',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_emis_hist = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_emis_hist','Emissivity of high ice st clouds',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_iwp_hist = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_iwp_hist','Ice water path of high ice st clouds','g/m2',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_deltaz_hist = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_deltaz_hist','geom thickness of high ice st clouds','m',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_rad_hist = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_rad_hist','ice crystals radius in high ice st clouds','µm',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))


 TYPE(ctrl_out), SAVE :: o_map_emis_Cb = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_emis_Cb','Emissivity of high Cb clouds',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

 TYPE(ctrl_out), SAVE :: o_map_pcld_Cb = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_pcld_Cb','cloud pressure of high Cb clouds','hPa',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

 TYPE(ctrl_out), SAVE :: o_map_tcld_Cb = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_tcld_Cb','cloud temperature of high Cb clouds','K',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))


 TYPE(ctrl_out), SAVE :: o_map_emis_Anv = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_emis_Anv','Emissivity of high Anv clouds',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

 TYPE(ctrl_out), SAVE :: o_map_pcld_Anv = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_pcld_Anv','cloud pressure of high Anv clouds','hPa',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_tcld_Anv = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_tcld_Anv','cloud temperature of high Anv clouds','K',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_emis_ThCi = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_emis_ThCi','Emissivity of high ThCi clouds',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_pcld_ThCi = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_pcld_ThCi','cloud pressure of high ThCi clouds','hPa',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_tcld_ThCi = &
  ctrl_out((/10,10,1,10,10,10,10,10,10,10/),&
  'map_tcld_ThCi','cloud temperature of high ThCi clouds','K',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

   TYPE(ctrl_out), SAVE :: o_map_ntot = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_ntot','total AIRS cloud fraction',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_hc = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_hc','high clouds AIRS cloud fraction',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_hist = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_hist','high clouds ice st AIRS cloud fraction',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE :: o_map_Cb = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_Cb','high clouds Cb AIRS cloud fraction',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

 TYPE(ctrl_out), SAVE :: o_map_ThCi = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_ThCi','high clouds ThCi AIRS cloud fraction',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

 TYPE(ctrl_out), SAVE :: o_map_Anv = &
  ctrl_out((/1,1,1,1,1,10,10,10,10,10/),&
  'map_Anv','high clouds Anv AIRS cloud fraction',' ',&
   (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)",&
      "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

! Fin Marine

  TYPE(ctrl_out), SAVE :: o_flat = ctrl_out((/ 5, 1, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'flat', 'Latent heat flux', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ptstar = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'ptstar', 'Air Surface Temperature', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_pt0 = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'pt0', 'Standard Air Surface Temperature', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_slp = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'slp', 'Sea Level Pressure', 'Pa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tsol = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'tsol', 'Surface Temperature', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_t2m = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    't2m', 'Temperature 2m', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_t2m_min = ctrl_out((/ 20, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    't2m_min', 'Temp 2m min', 'K', &
      (/ "t_min(X)", "t_min(X)", "t_min(X)", "t_min(X)", "t_min(X)", & 
         "t_min(X)", "t_min(X)", "t_min(X)", "t_min(X)", "t_min(X)" /))
  TYPE(ctrl_out), SAVE :: o_t2m_max = ctrl_out((/ 20, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    't2m_max', 'Temp 2m max', 'K', &
      (/ "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", &
         "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)" /))

  TYPE(ctrl_out), SAVE :: o_t2m_min_mon = ctrl_out((/ 1, 20, 20, 20, 20, 20, 20, 20, 20, 20 /), &
    't2m_min_mon', 'Monthly average min 2m temperature', 'K', &
      (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", & 
         "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))
  TYPE(ctrl_out), SAVE :: o_t2m_max_mon = ctrl_out((/ 1, 20, 20, 20, 20, 20, 20, 20, 20, 20 /), &
    't2m_max_mon', 'Monthly average max 2m temperature', 'K', &
      (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", &
         "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_t2m_srf = (/ &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    't2m_ter', "Temp 2m "//clnsurf(1), "K", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    't2m_lic', "Temp 2m "//clnsurf(2), "K", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    't2m_oce', "Temp 2m "//clnsurf(3), "K", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    't2m_sic', "Temp 2m "//clnsurf(4), "K", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE :: o_nt2mout = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'nt2mout', 'Nbt2m out of range complete computation', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_nq2mout = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'nq2mout', 'Nbq2m out of range complete computation', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_nu2mout = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'nu2mout', 'Nbu2m out of range complete computation', '-', (/ ('', i=1, 10) /))

  TYPE(ctrl_out), SAVE :: o_nt2moutfg = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'nt2moutfg', 'Nbt2m out of range complete/fgRi1 computation', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_nq2moutfg = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'nq2moutfg', 'Nbq2m out of range complete/fgRi1 computation', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_nu2moutfg = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'nu2moutfg', 'Nbu2m out of range complete/fgRi1 computation', '-', (/ ('', i=1, 10) /))

  TYPE(ctrl_out), SAVE :: o_gusts = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'gusts', 'surface gustiness', 'm2/s2', (/ ('', i=1, 10) /))

  TYPE(ctrl_out), SAVE :: o_wind10m = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'wind10m', '10-m wind speed', 'm/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wind100m = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wind100m', '100-m wind speed', 'm/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_loadfactor_wind_onshore = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'woncfr', 'Onshore Wind Capacity factor', 'kW/kW_installed', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_loadfactor_wind_offshore = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wofcfr', 'Offshore Wind Capacity factor', 'kW/kW_installed', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wind10max = ctrl_out((/ 10, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wind10max', '10m wind speed max', 'm/s', &
    (/ "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", & 
       "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)" /))

  TYPE(ctrl_out), SAVE :: o_sicf = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sicf', 'Sea-ice fraction', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_q2m = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'q2m', 'Specific humidity 2m', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ustar = ctrl_out((/ 1, 1, 10, 5, 10, 10, 11, 11, 11, 11/), &
    'ustar', 'Friction velocity', 'm/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_u10m = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'u10m', 'Vent zonal 10m', 'm/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_v10m = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'v10m', 'Vent meridien 10m', 'm/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_psol = ctrl_out((/ 1, 1, 1, 5, 10, 10, 11, 11, 11, 11/), &
    'psol', 'Surface Pressure', 'Pa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_qsurf = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'qsurf', 'Surface Air humidity', 'kg/kg', (/ ('', i=1, 10) /))

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_ustar_srf     = (/ &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'ustar_ter', &
      "Friction velocity "//clnsurf(1),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'ustar_lic', &
      "Friction velocity "//clnsurf(2),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'ustar_oce', &
      "Friction velocity "//clnsurf(3),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'ustar_sic', &
      "Friction velocity "//clnsurf(4),"m/s", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(5) :: o_wstar         = (/ &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'wstar_ter', &
      "Friction velocity "//clnsurf(1),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'wstar_lic', &
      "Friction velocity "//clnsurf(2),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'wstar_oce', &
      "Friction velocity "//clnsurf(3),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'wstar_sic', &
      "Friction velocity "//clnsurf(4),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 5, 5, 10, 10, 10, 10, 11, 11, 11, 11/),'wstar', &
      "w* convective velocity "//clnsurf(4),"m/s", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_u10m_srf     = (/ &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'u10m_ter', &
      "Vent Zonal 10m "//clnsurf(1),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'u10m_lic', &
      "Vent Zonal 10m "//clnsurf(2),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'u10m_oce', &
      "Vent Zonal 10m "//clnsurf(3),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'u10m_sic', &
      "Vent Zonal 10m "//clnsurf(4),"m/s", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_v10m_srf     = (/ &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'v10m_ter', &
      "Vent meredien 10m "//clnsurf(1),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'v10m_lic', &
      "Vent meredien 10m "//clnsurf(2),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'v10m_oce', &
      "Vent meredien 10m "//clnsurf(3),"m/s", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'v10m_sic', &
      "Vent meredien 10m "//clnsurf(4),"m/s", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE :: o_qsol = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'qsol', 'Soil watter content', 'mm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ndayrain = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ndayrain', 'Number of dayrain(liq+sol)', '-', &
      (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)" /))
  TYPE(ctrl_out), SAVE :: o_rain_fall = ctrl_out((/ 1, 1, 1, 10, 5, 10, 11, 11, 11, 11/), &
    'rain_fall', 'Precip Totale liq', 'kg/(s*m2)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rain_con = ctrl_out((/ 7, 7, 7, 10, 7, 10, 11, 11, 11, 11/), &
    'rain_con', 'Precip liq conv.', 'kg/(s*m2)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_precip = ctrl_out((/ 1, 1, 1, 10, 5, 10, 11, 11, 11, 11/), &
    'precip', 'Precip Totale liq+sol', 'kg/(s*m2)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_plul = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'plul', 'Large-scale Precip.', 'kg/(s*m2)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_plun = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'plun', 'Numerical Precip.', 'kg/(s*m2)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_pluc = ctrl_out((/ 1, 1, 1, 10, 5, 10, 11, 11, 11, 11/), &
    'pluc', 'Convective Precip.', 'kg/(s*m2)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_snow = ctrl_out((/ 1, 1, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'snow', 'Snow fall', 'kg/(s*m2)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_evap = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'evap', 'Evaporat', 'kg/(s*m2)', (/ ('', i=1, 10) /))

  TYPE(ctrl_out), SAVE :: o_sens_prec_liq_oce = ctrl_out((/ 5, 5, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'sens_rain_oce', 'Sensible heat flux of liquid prec. over ocean', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sens_prec_liq_sic = ctrl_out((/ 5, 5, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'sens_rain_sic', 'Sensible heat flux of liquid prec. over seaice', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sens_prec_sol_oce = ctrl_out((/ 5, 5, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'sens_snow_oce', 'Sensible heat flux of solid prec. over ocean', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sens_prec_sol_sic = ctrl_out((/ 5, 5, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'sens_snow_sic', 'Sensible heat flux of solid prec. over seaice', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lat_prec_liq_oce = ctrl_out((/ 5, 5, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'lat_rain_oce', 'Latent heat flux of liquid prec. over ocean', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lat_prec_liq_sic = ctrl_out((/ 5, 5, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'lat_rain_sic', 'Latent heat flux of liquid prec. over seaice', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lat_prec_sol_oce = ctrl_out((/ 5, 5, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'lat_snow_oce', 'Latent heat flux of solid prec. over ocean', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lat_prec_sol_sic = ctrl_out((/ 5, 5, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'lat_snow_sic', 'Latent heat flux of solid prec. over seaice', 'W/m2', (/ ('', i=1, 10) /))


  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_evap_srf     = (/ &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'evap_ter', &
      "evaporation at surface "//clnsurf(1),"kg/(s*m2)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'evap_lic', &
      "evaporation at surface "//clnsurf(2),"kg/(s*m2)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'evap_oce', &
      "evaporation at surface "//clnsurf(3),"kg/(s*m2)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'evap_sic', &
      "evaporation at surface "//clnsurf(4),"kg/(s*m2)", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE :: o_msnow = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'msnow', 'Surface snow amount', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_fsnow = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'fsnow', 'Surface snow area fraction', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tops = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tops', 'Solar rad. at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tops0 = ctrl_out((/ 1, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tops0', 'CS Solar rad. at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_topl = ctrl_out((/ 1, 1, 10, 5, 10, 10, 11, 11, 11, 11/), &
    'topl', 'IR rad. at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_topl0 = ctrl_out((/ 1, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'topl0', 'IR rad. at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWupTOA = ctrl_out((/ 1, 4, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'SWupTOA', 'SWup at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWupTOAclr = ctrl_out((/ 1, 4, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'SWupTOAclr', 'SWup clear sky at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWupTOAcleanclr = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'SWupTOAcleanclr', 'SWup clear sky clean (no aerosol) at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWdnTOA = ctrl_out((/ 1, 4, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'SWdnTOA', 'SWdn at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWdnTOAclr = ctrl_out((/ 1, 4, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'SWdnTOAclr', 'SWdn clear sky at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_nettop = ctrl_out((/ 1, 4, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'nettop', 'Net dn radiatif flux at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWup200 = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'SWup200', 'SWup at 200mb', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWup200clr = ctrl_out((/ 10, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'SWup200clr', 'SWup clear sky at 200mb', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWdn200 = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'SWdn200', 'SWdn at 200mb', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWdn200clr = ctrl_out((/ 10, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'SWdn200clr', 'SWdn clear sky at 200mb', 'W/m2', (/ ('', i=1, 10) /))

  ! arajouter
  !  type(ctrl_out),save :: o_LWupTOA     = ctrl_out((/ 1, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'LWupTOA', &
  !    (/ ('', i=1, 10) /))
  !  type(ctrl_out),save :: o_LWupTOAclr  = ctrl_out((/ 1, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'LWupTOAclr', &
  !    (/ ('', i=1, 10) /))
  !  type(ctrl_out),save :: o_LWdnTOA     = ctrl_out((/ 1, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'LWdnTOA', &
  !    (/ ('', i=1, 10) /))
  !  type(ctrl_out),save :: o_LWdnTOAclr  = ctrl_out((/ 1, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'LWdnTOAclr', &
  !    (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWup200 = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'LWup200', 'LWup at 200mb', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWup200clr = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'LWup200clr', 'LWup clear sky at 200mb', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWdn200 = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'LWdn200', 'LWdn at 200mb', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWdn200clr = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'LWdn200clr', 'LWdn clear sky at 200mb', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sols = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sols', 'Solar rad. at surf.', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sols0 = ctrl_out((/ 1, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sols0', 'Solar rad. at surf.', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_soll = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'soll', 'IR rad. at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_soll0 = ctrl_out((/ 1, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'soll0', 'IR rad. at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_radsol = ctrl_out((/ 1, 7, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'radsol', 'Rayonnement au sol', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWupSFC = ctrl_out((/ 1, 4, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'SWupSFC', 'SWup at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWupSFCclr = ctrl_out((/ 1, 4, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'SWupSFCclr', 'SWup clear sky at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWupSFCcleanclr = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'SWupSFCcleanclr', 'SWup clear sky clean (no aerosol) at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWdnSFC = ctrl_out((/ 1, 1, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'SWdnSFC', 'SWdn at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWdnSFCclr = ctrl_out((/ 1, 4, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'SWdnSFCclr', 'SWdn clear sky at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SWdnSFCcleanclr = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'SWdnSFCcleanclr', 'SWdn clear sky clean (no aerosol) at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWupSFC = ctrl_out((/ 1, 4, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'LWupSFC', 'Upwd. IR rad. at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWupSFCclr = ctrl_out((/ 1, 4, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'LWupSFCclr', 'CS Upwd. IR rad. at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWdnSFC = ctrl_out((/ 1, 4, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'LWdnSFC', 'Down. IR rad. at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWdnSFCclr = ctrl_out((/ 1, 4, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'LWdnSFCclr', 'Down. CS IR rad. at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWupTOAcleanclr = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'LWupTOAcleanclr', 'Upward CS clean (no aerosol) IR rad. at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWdnSFCcleanclr = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'LWdnSFCcleanclr', 'Downward CS clean (no aerosol) IR rad. at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_bils = ctrl_out((/ 1, 2, 10, 5, 10, 10, 11, 11, 11, 11/), &
    'bils', 'Surf. total heat flux', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_bils_tke = ctrl_out((/ 1, 2, 10, 5, 10, 10, 11, 11, 11, 11/), &
    'bils_tke', 'Surf. total heat flux', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_bils_diss = ctrl_out((/ 1, 2, 10, 5, 10, 10, 11, 11, 11, 11/), &
    'bils_diss', 'Surf. total heat flux', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_bils_ec = ctrl_out((/ 1, 2, 10, 5, 10, 10, 11, 11, 11, 11/), &
    'bils_ec', 'Surf. total heat flux correction', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_bils_ech = ctrl_out((/ 1, 2, 10, 5, 10, 10, 11, 11, 11, 11/), &
    'bils_ech', 'Surf. total heat flux correction', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_bils_kinetic = ctrl_out((/ 1, 2, 10, 5, 10, 10, 11, 11, 11, 11/), &
    'bils_kinetic', 'Surf. total heat flux', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_bils_enthalp = ctrl_out((/ 1, 2, 10, 5, 10, 10, 11, 11, 11, 11/), &
    'bils_enthalp', 'Surf. total heat flux', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_bils_latent = ctrl_out((/ 1, 2, 10, 5, 10, 10, 11, 11, 11, 11/), &
    'bils_latent', 'Surf. total heat flux', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sens = ctrl_out((/ 1, 1, 10, 10, 5, 10, 11, 11, 11, 11/), &
    'sens', 'Sensible heat flux', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_fder = ctrl_out((/ 1, 2, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'fder', 'Heat flux derivation', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ffonte = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ffonte', 'Thermal flux for snow melting', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_fqcalving = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'fqcalving', 'Ice Calving', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_fqfonte = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'fqfonte', 'Land ice melt', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_mrroli = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'mrroli', 'Runoff flux over land ice', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_runofflic = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'runofflic', 'Land ice melt to ocean', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_taux = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'taux', 'Zonal wind stress', 'Pa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tauy = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tauy', 'Meridional wind stress', 'Pa', (/ ('', i=1, 10) /))

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_taux_srf = (/           &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'taux_ter',             &
      "Zonal wind stress"//clnsurf(1), "Pa", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'taux_lic',             &
      "Zonal wind stress"//clnsurf(2), "Pa", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'taux_oce',             &
      "Zonal wind stress"//clnsurf(3), "Pa", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'taux_sic',             &
      "Zonal wind stress"//clnsurf(4), "Pa", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_tauy_srf     = (/             &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'tauy_ter',                   &
      "Meridional wind stress "//clnsurf(1),"Pa", (/ ('', i=1, 10) /)),  &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'tauy_lic',                   &
      "Meridional wind stress "//clnsurf(2),"Pa", (/ ('', i=1, 10) /)),  &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'tauy_oce',                   &
      "Meridional wind stress "//clnsurf(3),"Pa", (/ ('', i=1, 10) /)),  &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'tauy_sic',                   &
      "Meridional wind stress "//clnsurf(4),"Pa", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_pourc_srf    = (/ &
      ctrl_out((/ 1, 7, 10, 10, 10, 10, 11, 11, 11, 11/),'pourc_ter',      &
      "% "//clnsurf(1),"%", (/ ('', i=1, 10) /)),            &
      ctrl_out((/ 1, 7, 10, 10, 10, 10, 11, 11, 11, 11/),'pourc_lic',      &
      "% "//clnsurf(2),"%", (/ ('', i=1, 10) /)),            &
      ctrl_out((/ 1, 7, 10, 10, 10, 10, 11, 11, 11, 11/),'pourc_oce',      &
      "% "//clnsurf(3),"%", (/ ('', i=1, 10) /)),            &
      ctrl_out((/ 1, 7, 10, 10, 10, 10, 11, 11, 11, 11/),'pourc_sic',      &
      "% "//clnsurf(4),"%", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_fract_srf    = (/ &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'fract_ter',      &
      "Fraction "//clnsurf(1),"1", (/ ('', i=1, 10) /)),     &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'fract_lic',      &
      "Fraction "//clnsurf(2),"1", (/ ('', i=1, 10) /)),     &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'fract_oce',      &
      "Fraction "//clnsurf(3),"1", (/ ('', i=1, 10) /)),     &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'fract_sic',      &
      "Fraction "//clnsurf(4),"1", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_tsol_srf     = (/ &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'tsol_ter',       &
      "Temperature "//clnsurf(1),"K", (/ ('', i=1, 10) /)),  &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'tsol_lic',       &
      "Temperature "//clnsurf(2),"K", (/ ('', i=1, 10) /)),  &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'tsol_oce',       &
      "Temperature "//clnsurf(3),"K", (/ ('', i=1, 10) /)),  &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'tsol_sic',       &
      "Temperature "//clnsurf(4),"K", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_evappot_srf  = (/ &
      ctrl_out((/ 1, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'evappot_ter',    &
      "Temperature"//clnsurf(1),"K", (/ ('', i=1, 10) /)),   &
      ctrl_out((/ 4, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'evappot_lic',    &
      "Temperature"//clnsurf(2),"K", (/ ('', i=1, 10) /)),   &
      ctrl_out((/ 4, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'evappot_oce',    &
      "Temperature"//clnsurf(3),"K", (/ ('', i=1, 10) /)),   &
      ctrl_out((/ 4, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'evappot_sic',    &
      "Temperature"//clnsurf(4),"K", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_sens_srf     = (/          &
      ctrl_out((/ 1, 6, 10, 7, 10, 10, 11, 11, 11, 11/),'sens_ter',                 &
      "Sensible heat flux "//clnsurf(1),"W/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 7, 10, 10, 11, 11, 11, 11/),'sens_lic',                 &
      "Sensible heat flux "//clnsurf(2),"W/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 7, 10, 10, 11, 11, 11, 11/),'sens_oce',                 &
      "Sensible heat flux "//clnsurf(3),"W/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 7, 10, 10, 11, 11, 11, 11/),'sens_sic',                 &
      "Sensible heat flux "//clnsurf(4),"W/m2", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_lat_srf      = (/        &
      ctrl_out((/ 1, 6, 10, 7, 10, 10, 11, 11, 11, 11/),'lat_ter',                &
      "Latent heat flux "//clnsurf(1),"W/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 7, 10, 10, 11, 11, 11, 11/),'lat_lic',                &
      "Latent heat flux "//clnsurf(2),"W/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 7, 10, 10, 11, 11, 11, 11/),'lat_oce',                &
      "Latent heat flux "//clnsurf(3),"W/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 6, 10, 7, 10, 10, 11, 11, 11, 11/),'lat_sic',                &
      "Latent heat flux "//clnsurf(4),"W/m2", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_flw_srf      = (/ &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'flw_ter',       &
      "LW "//clnsurf(1),"W/m2", (/ ('', i=1, 10) /)),        &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'flw_lic',       &
      "LW "//clnsurf(2),"W/m2", (/ ('', i=1, 10) /)),        &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'flw_oce',       &
      "LW "//clnsurf(3),"W/m2", (/ ('', i=1, 10) /)),        &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'flw_sic',       &
      "LW "//clnsurf(4),"W/m2", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_fsw_srf      = (/ &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'fsw_ter',       &
      "SW "//clnsurf(1),"W/m2", (/ ('', i=1, 10) /)),        &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'fsw_lic',       &
      "SW "//clnsurf(2),"W/m2", (/ ('', i=1, 10) /)),        &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'fsw_oce',       &
      "SW "//clnsurf(3),"W/m2", (/ ('', i=1, 10) /)),        &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'fsw_sic',       &
      "SW "//clnsurf(4),"W/m2", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_wbils_srf    = (/ &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wbils_ter',     &
      "Bilan sol "//clnsurf(1),"W/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wbils_lic',     &
      "Bilan sol "//clnsurf(2),"W/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wbils_oce',     &
      "Bilan sol "//clnsurf(3),"W/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wbils_sic',     &
      "Bilan sol "//clnsurf(4),"W/m2", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_wbilo_srf    = (/      &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wbilo_ter',          &
      "Bilan eau "//clnsurf(1),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wbilo_lic',          &
      "Bilan eau "//clnsurf(2),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wbilo_oce',          &
      "Bilan eau "//clnsurf(3),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wbilo_sic',          &
      "Bilan eau "//clnsurf(4),"kg/(m2*s)", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_wevap_srf    = (/      &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wevap_ter',          &
      "Evap eau "//clnsurf(1),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wevap_lic',          &
      "Evap eau "//clnsurf(2),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wevap_oce',          &
      "Evap eau "//clnsurf(3),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wevap_sic',          &
      "Evap eau "//clnsurf(4),"kg/(m2*s)", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_wrain_srf    = (/      &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wrain_ter',          &
      "Pluie eau "//clnsurf(1),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wrain_lic',          &
      "Pluie eau "//clnsurf(2),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wrain_oce',          &
      "Pluie eau "//clnsurf(3),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wrain_sic',          &
      "Pluie eau "//clnsurf(4),"kg/(m2*s)", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_wsnow_srf    = (/      &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wsnow_ter',          &
      "Neige eau "//clnsurf(1),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wsnow_lic',          &
      "Neige eau "//clnsurf(2),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wsnow_oce',          &
      "Neige eau "//clnsurf(3),"kg/(m2*s)", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'wsnow_sic',          &
      "Neige eau "//clnsurf(4),"kg/(m2*s)", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE :: o_cdrm = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cdrm', 'Momentum drag coef.', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cdrh = ctrl_out((/ 1, 10, 10, 7, 10, 10, 11, 11, 11, 11/), &
    'cdrh', 'Heat drag coef.', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldl = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldl', 'Low-level cloudiness', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldm = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldm', 'Mid-level cloudiness', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldh = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldh', 'High-level cloudiness', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldt = ctrl_out((/ 1, 1, 2, 10, 5, 10, 11, 11, 11, 11/), &
    'cldt', 'Total cloudiness', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_JrNt = ctrl_out((/ 1, 1, 10, 7, 10, 10, 11, 11, 11, 11/), &
    'JrNt', '1 if Day 0 if Night', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldhjn = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldhjn', 'High-level cloudiness Day', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldmjn = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &      
    'cldmjn', 'Mid-level cloudiness day', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldljn = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &       
    'cldljn', 'Low-level cloudiness day', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldtjn = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &     
    'cldtjn', 'Total cloudiness day', '-', (/ ('', i=1, 10) /))
 
  TYPE(ctrl_out), SAVE :: o_cldq = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldq', 'Cloud water path', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lwp = ctrl_out((/ 1, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'lwp', 'Cloud liquid water path', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_iwp = ctrl_out((/ 1, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'iwp', 'Cloud ice water path', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ue = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ue', 'Zonal dry static energy transport', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ve = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    've', 'Merid dry static energy transport', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_uq = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'uq', 'Zonal humidity transport', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_vq = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'vq', 'Merid humidity transport', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_uwat = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'uwat', 'Zonal total water transport', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_vwat = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'vwat', 'Merid total water transport', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cape = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cape', 'Conv avlbl pot ener', 'J/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_pbase = ctrl_out((/ 1, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'pbase', 'Cld base pressure', 'Pa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ptop = ctrl_out((/ 1, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ptop', 'Cld top pressure', 'Pa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_fbase = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'fbase', 'Cld base mass flux', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_plcl = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'plcl', 'Lifting Condensation Level', 'hPa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_plfc = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'plfc', 'Level of Free Convection', 'hPa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wbeff = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wbeff', 'Conv. updraft velocity at LFC (<100)', 'm/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_convoccur = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'convoccur', 'Convective occurence', '', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_prw = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'prw', 'Precipitable water', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_prlw = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'prlw', 'Precipitable liquid water', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_prsw = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'prsw', 'Precipitable solid water', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_s_pblh = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    's_pblh', 'Boundary Layer Height', 'm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_s_pblt = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    's_pblt', 't at Boundary Layer Height', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_s_lcl = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    's_lcl', 'Condensation level', 'm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_s_therm = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    's_therm', 'Exces du thermique', 'K', (/ ('', i=1, 10) /))
  !IM : Les champs suivants (s_capCL, s_oliqCL, s_cteiCL, s_trmb1, s_trmb2, s_trmb3) ne sont pas definis dans HBTM.F
  ! type(ctrl_out),save :: o_s_capCL      = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'s_capCL', &
!    (/ ('', i=1, 10) /))
  ! type(ctrl_out),save :: o_s_oliqCL     = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'s_oliqCL', &
!    (/ ('', i=1, 10) /))
  ! type(ctrl_out),save :: o_s_cteiCL     = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'s_cteiCL', &
!    (/ ('', i=1, 10) /))
  ! type(ctrl_out),save :: o_s_trmb1      = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'s_trmb1', &
!    (/ ('', i=1, 10) /))
  ! type(ctrl_out),save :: o_s_trmb2      = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'s_trmb2', &
!    (/ ('', i=1, 10) /))
  ! type(ctrl_out),save :: o_s_trmb3      = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'s_trmb3', &
    !(/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_slab_bils = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'slab_bils', 'flux atmos - slab ponderes foce', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_slab_bilg = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'slab_bilg', 'flux glace - slab ponderes fsic', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_slab_qflux = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'slab_qflux', 'Correction flux slab', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tslab = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tslab', 'Temperature ocean slab', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tslab1 = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11 /), &
    'tslab1', 'Temperature ocean slab', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tslab2 = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11 /), &
    'tslab2', 'Temperature ocean slab', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_slab_tice = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'slab_tice', 'Temperature banquise slab', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_slab_sic = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'seaice', 'Epaisseur banquise slab', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_slab_hdiff = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'slab_hdiff', 'Horizontal diffusion', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_slab_gm = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11 /), &
    'slab_gm', 'GM eddy advection', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_slab_ekman = ctrl_out((/ 1, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'slab_ekman', 'Ekman heat transport', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ale_bl = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'ale_bl', 'ALE BL', 'm2/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_alp_bl = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'alp_bl', 'ALP BL', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ale_wk = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'ale_wk', 'ALE WK', 'm2/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_alp_wk = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'alp_wk', 'ALP WK', 'W/m2', (/ ('', i=1, 10) /))
!!!
!nrlmd+jyg<
  type(ctrl_out),save :: o_dtvdf_x        = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtvdf_x', ' dtvdf off_wake','K/s', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_dtvdf_w        = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtvdf_w', ' dtvdf within_wake','K/s', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_dqvdf_x        = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqvdf_x', ' dqvdf off_wake','kg/kg/s', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_dqvdf_w        = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqvdf_w', ' dqvdf within_wake','kg/kg/s', (/ ('', i=1, 10) /))
!!
  type(ctrl_out),save :: o_sens_x        = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'sens_x', 'sens off_wake', 'W/m2', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_sens_w        = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'sens_w', 'sens within_wake', 'W/m2', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_flat_x        = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'flat_x', 'flat off_wake', 'W/m2', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_flat_w        = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'flat_w', 'flat within_wake', 'W/m2', (/ ('', i=1, 10) /))
!!
  type(ctrl_out),save :: o_delta_tsurf    = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'delta_tsurf', 'Temperature difference (w-x)', 'K', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_cdragh_x       = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'cdragh_x', 'cdragh off-wake', '', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_cdragh_w       = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'cdragh_w', 'cdragh within-wake', '', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_cdragm_x       = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'cdragm_x', 'cdragm off-wake', '', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_cdragm_w       = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'cdragm_w', 'cdrgam within-wake', '', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_kh             = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'kh', 'Kh', 'kg/s/m2', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_kh_x           = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'kh_x', 'Kh off-wake', 'kg/s/m2', (/ ('', i=1, 10) /))
  type(ctrl_out),save :: o_kh_w           = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
'kh_w', 'Kh within-wake', 'kg/s/m2', (/ ('', i=1, 10) /))
!>nrlmd+jyg
!!!
  TYPE(ctrl_out), SAVE :: o_ale = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'ale', 'ALE', 'm2/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_alp = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'alp', 'ALP', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cin = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'cin', 'Convective INhibition', 'm2/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wape = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'wape', '', 'm2/s2', (/ ('', i=1, 10) /))

!!! nrlmd le 10/04/2012

!-------Spectre de thermiques de type 2 au LCL
  TYPE(ctrl_out), SAVE :: o_n2 = ctrl_out((/ 1, 6, 6, 6, 10, 10, 11, 11, 11, 11/), &
    'n2', 'Nombre de panaches de type 2', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_s2 = ctrl_out((/ 1, 6, 6, 6, 10, 10, 11, 11, 11, 11/), &
    's2', 'Surface moyenne des panaches de type 2', 'm2', (/ ('', i=1, 10) /))
              
!-------Déclenchement stochastique
  TYPE(ctrl_out), SAVE :: o_proba_notrig = ctrl_out((/ 1, 6, 6, 6, 10, 10, 11, 11, 11, 11/), &
    'proba_notrig', 'Probabilite de non-declenchement', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_random_notrig = ctrl_out((/ 1, 6, 6, 6, 10, 10, 11, 11, 11, 11/), &
    'random_notrig', 'Tirage aleatoire de non-declenchement', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ale_bl_stat = ctrl_out((/ 1, 6, 6, 6, 10, 10, 11, 11, 11, 11/), &
    'ale_bl_stat', 'ALE_BL_STAT', 'm2/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ale_bl_trig = ctrl_out((/ 1, 6, 6, 6, 10, 10, 11, 11, 11, 11/), &
    'ale_bl_trig', 'ALE_BL_STAT + Condition S>Sthreshold', 'm2/s2', (/ ('', i=1, 10) /))

!-------Fermeture statistique
  TYPE(ctrl_out), SAVE :: o_alp_bl_det = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'alp_bl_det', 'ALP_BL_DET', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_alp_bl_fluct_m = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'alp_bl_fluct_m', 'ALP_BL_FLUCT_M', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_alp_bl_fluct_tke = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'alp_bl_fluct_tke', 'ALP_BL_FLUCT_TKE', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_alp_bl_conv = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'alp_bl_conv', 'ALP_BL_CONV', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_alp_bl_stat = ctrl_out((/ 1, 1, 1, 10, 10, 10, 11, 11, 11, 11/), &
    'alp_bl_stat', 'ALP_BL_STAT', 'W/m2', (/ ('', i=1, 10) /))

!!! fin nrlmd le 10/04/2012

  ! Champs interpolles sur des niveaux de pression ??? a faire correctement

  TYPE(ctrl_out), SAVE, DIMENSION(7) :: o_uSTDlevs     = (/                    &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'u850', "Zonal wind 850hPa", "m/s",     &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'u700', "Zonal wind 700hPa", "m/s",     &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'u500', "Zonal wind 500hPa", "m/s",     &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'u200', "Zonal wind 200hPa", "m/s",     &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'u100', "Zonal wind 100hPa", "m/s",     &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'u50', "Zonal wind 50hPa", "m/s",     &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'u10', "Zonal wind 10hPa", "m/s",     &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(7) :: o_vSTDlevs     = (/                     &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'v850', "Meridional wind 850hPa", "m/s", &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)),  &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'v700', "Meridional wind 700hPa", "m/s", &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)),  &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'v500', "Meridional wind 500hPa", "m/s", &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)),  &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'v200', "Meridional wind 200hPa", "m/s", &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)),  &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'v100', "Meridional wind 100hPa", "m/s", &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)),  &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'v50', "Meridional wind 50hPa", "m/s",  &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)),  &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'v10', "Meridional wind 10hPa", "m/s",  &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(7) :: o_wSTDlevs     = (/                    &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'w850', "Vertical wind 850hPa", "Pa/s", &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'w700', "Vertical wind 700hPa", "Pa/s", &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'w500', "Vertical wind 500hPa", "Pa/s", &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'w200', "Vertical wind 200hPa", "Pa/s", &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'w100', "Vertical wind 100hPa", "Pa/s", &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'w50', "Vertical wind 50hPa", "Pa/s",  &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'w10', "Vertical wind 10hPa", "Pa/s",  &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(7) :: o_tSTDlevs     = (/                    &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'t850', "Temperature 850hPa", "K",      &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'t700', "Temperature 700hPa", "K",      &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'t500', "Temperature 500hPa", "K",      &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'t200', "Temperature 200hPa", "K",      &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'t100', "Temperature 100hPa", "K",      &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'t50',  "Temperature 50hPa", "K",      &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'t10',  "Temperature 10hPa", "K",      &
      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(7) :: o_qSTDlevs     = (/ &
       ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'q850', &
       "Specific humidity 850hPa", "kg/kg", &
       (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
       ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'q700', &
       "Specific humidity 700hPa", "kg/kg", &
       (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
       ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'q500', &
       "Specific humidity 500hPa", "kg/kg", &
       (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
       ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'q200', &
       "Specific humidity 200hPa", "kg/kg", &
       (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
       ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'q100', &
       "Specific humidity 100hPa", "kg/kg", &
       (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
       ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'q50', &
       "Specific humidity 50hPa", "kg/kg", &
       (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
       ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'q10', &
       "Specific humidity 10hPa", "kg/kg", &
       (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(7) :: o_zSTDlevs   = (/                           &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'z850', "Geopotential height 850hPa",        &
      "m", (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'z700', "Geopotential height 700hPa",        &
      "m", (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'z500', "Geopotential height 500hPa",        &
      "m", (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'z200', "Geopotential height 200hPa",        &
      "m", (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'z100', "Geopotential height 100hPa",        &
      "m", (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'z50', "Geopotential height 50hPa",         &
      "m", (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)), &
      ctrl_out((/ 1, 7, 7, 10, 10, 10, 11, 11, 11, 11/),'z10', "Geopotential height 10hPa",         &
      "m", (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /)) /)

  TYPE(ctrl_out), SAVE :: o_t_oce_sic = ctrl_out((/ 1, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    't_oce_sic', 'Temp mixte oce-sic', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_weakinv = ctrl_out((/ 10, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'weakinv', 'Weak inversion', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dthmin = ctrl_out((/ 10, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dthmin', 'dTheta mini', 'K/m', (/ ('', i=1, 10) /))

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_u10_srf      = (/ &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'u10_ter', "", "", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'u10_lic', "", "", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'u10_oce', "", "", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'u10_sic', "", "", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_v10_srf      = (/ &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'v10_ter', "", "", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'v10_lic', "", "", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'v10_oce', "", "", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'v10_sic', "", "", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE :: o_cldtau = ctrl_out((/ 10, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldtau', 'Cloud optical thickness', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldemi = ctrl_out((/ 10, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldemi', 'Cloud optical emissivity', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rh2m = ctrl_out((/ 5, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rh2m', 'Relative humidity at 2m', '%', (/ ('', i=1, 10) /))
!  TYPE(ctrl_out), SAVE :: o_rh2m_min = ctrl_out((/ 10, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
!    'rh2m_min', 'Min Relative humidity at 2m', '%',                        &
!      (/ 't_min(X)', 't_min(X)', 't_min(X)', 't_min(X)', 't_min(X)', & 
!         't_min(X)', 't_min(X)', 't_min(X)', 't_min(X)', 't_min(X)' /))
!  TYPE(ctrl_out), SAVE :: o_rh2m_max = ctrl_out((/ 10, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
!    'rh2m_max', 'Max Relative humidity at 2m', '%',                         &
!      (/ 't_max(X)', 't_max(X)', 't_max(X)', 't_max(X)', 't_max(X)', &
!         't_max(X)', 't_max(X)', 't_max(X)', 't_max(X)', 't_max(X)' /))
  TYPE(ctrl_out), SAVE :: o_qsat2m = ctrl_out((/ 10, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'qsat2m', 'Saturant humidity at 2m', '%', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tpot = ctrl_out((/ 10, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tpot', 'Surface air potential temperature', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tpote = ctrl_out((/ 10, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tpote', 'Surface air equivalent potential temperature', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tke = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tke ', 'TKE', 'm2/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tke_max = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tke_max', 'TKE max', 'm2/s2',                                  &
      (/ 't_max(X)', 't_max(X)', 't_max(X)', 't_max(X)', 't_max(X)', &
         't_max(X)', 't_max(X)', 't_max(X)', 't_max(X)', 't_max(X)' /))

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_tke_srf      = (/             &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'tke_ter',       &
      "Max Turb. Kinetic Energy "//clnsurf(1),"m2/s2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'tke_lic',       &
      "Max Turb. Kinetic Energy "//clnsurf(2),"m2/s2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'tke_oce',       &
      "Max Turb. Kinetic Energy "//clnsurf(3),"m2/s2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'tke_sic',       &
      "Max Turb. Kinetic Energy "//clnsurf(4),"m2/s2", (/ ('', i=1, 10) /)) /)
!FC
!  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_treedrg_srf      = (/             &
!      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'treedrg_ter',       &
!      "Drag from trees "//clnsurf(1),"-", (/ ('', i=1, 10) /)), &
!      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'treedrg_lic',       &
!      "Drag from trees "//clnsurf(2),"-", (/ ('', i=1, 10) /)), &
!      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'treedrg_oce',       &
!      "Drag from trees "//clnsurf(3),"-", (/ ('', i=1, 10) /)), &
!      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'treedrg_sic',       &
!      "Drag from trees "//clnsurf(4),"-", (/ ('', i=1, 10) /)) /)
!FC

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_l_mixmin      = (/             &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'l_mixmin_ter',       &
      "PBL mixing length "//clnsurf(1),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'l_mixmin_lic',       &
      "PBL mixing length "//clnsurf(2),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'l_mixmin_oce',       &
      "PBL mixing length "//clnsurf(3),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'l_mixmin_sic',       &
      "PBL mixing length "//clnsurf(4),"m", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_l_mix      = (/             &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'l_mix_ter',       &
      "min PBL mixing length "//clnsurf(1),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'l_mix_lic',       &
      "min PBL mixing length "//clnsurf(2),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'l_mix_oce',       &
      "min PBL mixing length "//clnsurf(3),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'l_mix_sic',       &
      "min PBL mixing length "//clnsurf(4),"m", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_tke_max_srf  = (/                          &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'tke_max_ter',                &
      "Max Turb. Kinetic Energy "//clnsurf(1),"-",                                   &
      (/ "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", &
         "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)" /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'tke_max_lic',                &
      "Max Turb. Kinetic Energy "//clnsurf(2),"-",                                   &
      (/ "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", &
         "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)" /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'tke_max_oce',                &
      "Max Turb. Kinetic Energy "//clnsurf(3),"-",                                   &
      (/ "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", &
         "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)" /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'tke_max_sic',                &
      "Max Turb. Kinetic Energy "//clnsurf(4),"-",                                   &
      (/ "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", &
         "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)" /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_dltpbltke_srf      = (/             &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'dltpbltke_ter',       &
      "TKE difference (w - x) "//clnsurf(1),"-", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'dltpbltke_lic',       &
      "TKE difference (w - x) "//clnsurf(2),"-", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'dltpbltke_oce',       &
      "TKE difference (w - x) "//clnsurf(3),"-", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 4, 10, 10, 10, 10, 11, 11, 11, 11/),'dltpbltke_sic',       &
      "TKE difference (w - x) "//clnsurf(4),"-", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE :: o_kz = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'kz', 'Kz melange', 'm2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_kz_max = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'kz_max', 'Kz melange max', 'm2/s',                                  &
      (/ 't_max(X)', 't_max(X)', 't_max(X)', 't_max(X)', 't_max(X)', &
         't_max(X)', "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)" /))
  TYPE(ctrl_out), SAVE :: o_SWnetOR = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'SWnetOR', 'Sfce net SW radiation OR', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_LWdownOR = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'LWdownOR', 'Sfce incident LW radiation OR', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_snowl = ctrl_out((/ 10, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'snowl', 'Solid Large-scale Precip.', 'kg/(m2*s)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cape_max = ctrl_out((/ 10, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cape_max', 'CAPE max.', 'J/kg',                                       &
      (/ "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", &
         "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)", "t_max(X)" /))
  TYPE(ctrl_out), SAVE :: o_solldown = ctrl_out((/ 10, 1, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'solldown', 'Down. IR rad. at surface', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtsvdfo = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtsvdfo', 'Boundary-layer dTs(o)', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtsvdft = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtsvdft', 'Boundary-layer dTs(t)', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtsvdfg = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtsvdfg', 'Boundary-layer dTs(g)', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtsvdfi = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtsvdfi', 'Boundary-layer dTs(g)', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_z0m = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'z0m', 'roughness length, momentum', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_z0h = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'z0h', 'roughness length, enthalpy', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_topswad = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'topswad', 'ADE at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_topswad0 = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'topswad0', 'ADE clear-sky at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_topswai = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'topswai', 'AIE at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_solswad = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'solswad', 'ADE at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_solswad0 = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'solswad0', 'ADE clear-sky at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_solswai = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'solswai', 'AIE at SFR', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_toplwad = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'toplwad', 'LW-ADE at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_toplwad0 = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'toplwad0', 'LW-ADE clear-sky at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_toplwai = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'toplwai', 'LW-AIE at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sollwad = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sollwad', 'LW-ADE at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sollwad0 = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sollwad0', 'LW-ADE clear-sky at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sollwai = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sollwai', 'LW-AIE at SFR', 'W/m2', (/ ('', i=1, 10) /))

  TYPE(ctrl_out),SAVE,DIMENSION(naero_tot) :: o_tausumaero =                              &
       (/ ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(1),     &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(1),"1", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(2),        &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(2),"2", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(3),        &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(3),"3", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(4),        &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(4),"4", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(5),        &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(5),"5", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(6),        &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(6),"6", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(7),        &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(7),"7", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(8),        &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(8),"8", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(9),        &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(9),"9", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(10),       &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(10),"10", (/ ('', i=1, 10) /)),   &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(11),       &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(11),"11", (/ ('', i=1, 10) /)),   &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(12),       &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(12),"12", (/ ('', i=1, 10) /)),   &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(13),       &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(13),"13", (/ ('', i=1, 10) /)),   &
       ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'od550_'//name_aero_tau(14),       &
       "Aerosol Optical depth at 550 nm "//name_aero_tau(14),"14", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out),SAVE,DIMENSION(naero_tot-1) :: o_drytausumaero =                              &
       (/ ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(1),     &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(1),"1", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(2),        &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(2),"2", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(3),        &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(3),"3", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(4),        &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(4),"4", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(5),        &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(5),"5", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(6),        &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(6),"6", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(7),        &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(7),"7", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(8),        &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(8),"8", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(9),        &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(9),"9", (/ ('', i=1, 10) /)),     &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(10),       &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(10),"10", (/ ('', i=1, 10) /)),   &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(11),       &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(11),"11", (/ ('', i=1, 10) /)),   &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(12),       &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(12),"12", (/ ('', i=1, 10) /)),   &
       ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'dryod550_'//name_aero_tau(13),       &
       "Dry aerosol Optical depth at 550 nm "//name_aero_tau(13),"13", (/ ('', i=1, 10) /)) /)
!
  TYPE(ctrl_out), SAVE :: o_tausumaero_lw = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'od_10um_STRAT', 'Stratospheric Aerosol Optical depth at 10 um ', '1', (/ ('', i=1, 10) /))
!
  TYPE(ctrl_out), SAVE :: o_od443aer = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'od443aer', 'Total aerosol optical depth at 440nm', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_od550aer = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'od550aer', 'Total aerosol optical depth at 550nm', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dryod550aer = ctrl_out((/ 11, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dryod550aer', 'Total dry aerosol optical depth at 550nm', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_od865aer = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'od865aer', 'Total aerosol optical depth at 870nm', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_abs550aer = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'abs550aer', 'Absorption aerosol optical depth at 550nm', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_od550lt1aer = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'od550lt1aer', 'Fine mode optical depth', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sconcso4 = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sconcso4', 'Surface Concentration of Sulfate ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sconcno3 = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sconcno3', 'Surface Concentration of Nitrate ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sconcoa = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sconcoa', 'Surface Concentration of Organic Aerosol ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sconcbc = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sconcbc', 'Surface Concentration of Black Carbon ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sconcss = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sconcss', 'Surface Concentration of Sea Salt ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sconcdust = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sconcdust', 'Surface Concentration of Dust ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_concso4 = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'concso4', 'Concentration of Sulfate ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_concno3 = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'concno3', 'Concentration of Nitrate ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_concoa = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'concoa', 'Concentration of Organic Aerosol ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_concbc = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'concbc', 'Concentration of Black Carbon ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_concss = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'concss', 'Concentration of Sea Salt ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_concdust = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'concdust', 'Concentration of Dust ', 'kg/m3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_loadso4 = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'loadso4', 'Column Load of Sulfate ', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_loadoa = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'loadoa', 'Column Load of Organic Aerosol ', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_loadbc = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'loadbc', 'Column Load of Black Carbon ', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_loadss = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'loadss', 'Column Load of Sea Salt ', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_loaddust = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'loaddust', 'Column Load of Dust ', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_loadno3 = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'loadno3', 'Column Load of Nitrate ', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swtoaas_nat = ctrl_out((/ 11, 11, 1, 11, 11, 11, 11, 11, 11, 11/), &
    'swtoaas_nat', 'Natural aerosol radiative forcing all-sky at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swsrfas_nat = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swsrfas_nat', 'Natural aerosol radiative forcing all-sky at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swtoacs_nat = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swtoacs_nat', 'Natural aerosol radiative forcing clear-sky at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swsrfcs_nat = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swsrfcs_nat', 'Natural aerosol radiative forcing clear-sky at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swtoaas_ant = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swtoaas_ant', 'Anthropogenic aerosol radiative forcing all-sky at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swsrfas_ant = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swsrfas_ant', 'Anthropogenic aerosol radiative forcing all-sky at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swtoacs_ant = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swtoacs_ant', 'Anthropogenic aerosol radiative forcing clear-sky at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swsrfcs_ant = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swsrfcs_ant', 'Anthropogenic aerosol radiative forcing clear-sky at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swtoacf_nat = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swtoacf_nat', 'Natural aerosol impact on cloud radiative forcing at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swsrfcf_nat = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swsrfcf_nat', 'Natural aerosol impact on cloud radiative forcing  at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swtoacf_ant = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swtoacf_ant', 'Anthropogenic aerosol impact on cloud radiative forcing at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swsrfcf_ant = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swsrfcf_ant', 'Anthropogenic aerosol impact on cloud radiative forcing at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swtoacf_zero = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swtoacf_zero', 'Cloud radiative forcing (allsky-clearsky fluxes) at TOA', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_swsrfcf_zero = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11/), &
    'swsrfcf_zero', 'Cloud radiative forcing (allsky-clearsky fluxes) at SRF', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldncl = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldncl', 'CDNC at top of liquid water cloud', 'm-3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_reffclwtop = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'reffclwtop', 'Droplet effective radius at top of liquid water cloud', 'm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldnvi = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldnvi', 'Column Integrated Cloud Droplet Number', 'm-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lcc = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'lcc', 'Cloud liquid fraction at top of cloud', '1', (/ ('', i=1, 10) /))

!--tropopause pressure
  TYPE(ctrl_out), SAVE :: o_p_tropopause = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'p_tropopause', 'Tropopause pressure', 'Pa', (/ ('', i=1, 10) /))
!--tropopause height
  TYPE(ctrl_out), SAVE :: o_z_tropopause = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'z_tropopause', 'Tropopause height', 'm', (/ ('', i=1, 10) /))
!--tropopause temperature
  TYPE(ctrl_out), SAVE :: o_t_tropopause = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    't_tropopause', 'Tropopause temperature', 'K', (/ ('', i=1, 10) /))
!--Added ThL
  TYPE(ctrl_out), SAVE :: o_col_O3_strato = ctrl_out((/2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'colO3_strat','Ozone stratospheric column', 'DU', (/('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_col_O3_tropo = ctrl_out((/2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'colO3_trop','Ozone tropospheric column', 'DU', (/('', i=1, 10) /))
!--end add ThL

!---CO2 fluxes for interactive CO2 configuration
  TYPE(ctrl_out), SAVE :: o_flx_co2_ff = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'flx_co2_ff', 'CO2 flux from fossil fuel and cement', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_flx_co2_bb = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'flx_co2_bb', 'CO2 flux from biomass burning', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_flx_co2_ocean = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'flx_co2_ocean', 'CO2 flux from the ocean', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_flx_co2_land = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'flx_co2_land', 'CO2 flux from the land', '1', (/ ('', i=1, 10) /))

#ifdef CPP_StratAer
!--extinction coefficient
  TYPE(ctrl_out), SAVE :: o_ext_strat_550 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'ext_strat_550', 'Strat. aerosol extinction coefficient at 550 nm', '1/m', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ext_strat_1020 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'ext_strat_1020', 'Strat. aerosol extinction coefficient at 1020 nm', '1/m', (/ ('', i=1, 10) /))
!--strat aerosol optical depth
  TYPE(ctrl_out), SAVE :: o_tau_strat_550 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'OD550_strat_only', 'Stratospheric Aerosol Optical depth at 550 nm ', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tau_strat_1020 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'OD1020_strat_only', 'Stratospheric Aerosol Optical depth at 1020 nm ', '1', (/ ('', i=1, 10) /))
!--chemistry
  TYPE(ctrl_out), SAVE :: o_R2SO4 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'R2SO4', 'H2SO4 mass fraction in aerosol', '%', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_OCS_lifetime = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'OCS_lifetime', 'OCS lifetime', 's', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_SO2_lifetime = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'SO2_lifetime', 'SO2 lifetime', 's', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_f_r_wet = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'f_r_wet', 'Conversion factor dry to wet aerosol radius', '-', (/ ('', i=1, 10) /))
!--budget  3D
  TYPE(ctrl_out), SAVE :: o_budg_3D_nucl = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_3D_nucl', 'H2SO4 nucleation mass flux', 'kg(S)/m2/layer/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_3D_cond_evap = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_3D_cond_evap', 'H2SO4 condensation/evaporation mass flux', 'kg(S)/m2/layer/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_3D_ocs_to_so2 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_3D_ocs_to_so2', 'OCS mass flux converted to SO2', 'kg(S)/m2/layer/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_3D_so2_to_h2so4 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_3D_so2_to_h2so4', 'SO2 mass flux converted to H2SO4', 'kg(S)/m2/layer/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_3D_backgr_ocs = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_3D_backgr_ocs', 'OCS background tendency', 'kg(S)/m2/layer/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_3D_backgr_so2 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_3D_backgr_so2', 'SO2 background tendency', 'kg(S)/m2/layer/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_vsed_aer = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'vsed_aer', 'Strat. aerosol sedimentation velocity (mass-weighted)', 'm/s', (/ ('', i=1, 10) /))
!--budget  2D
  TYPE(ctrl_out), SAVE :: o_budg_dep_dry_ocs = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_dep_dry_ocs',   'OCS dry deposition flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_dep_wet_ocs = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_dep_wet_ocs',   'OCS wet deposition flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_dep_dry_so2 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_dep_dry_so2',   'SO2 dry deposition flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_dep_wet_so2 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_dep_wet_so2',   'SO2 wet deposition flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_dep_dry_h2so4 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_dep_dry_h2so4', 'H2SO4 dry deposition flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_dep_wet_h2so4 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_dep_wet_h2so4', 'H2SO4 wet deposition flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_dep_dry_part = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_dep_dry_part', 'particle dry deposition flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_dep_wet_part = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_dep_wet_part', 'particle wet deposition flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_emi_ocs = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_emi_ocs', 'OCS emission flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_emi_so2 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_emi_so2', 'SO2 emission flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_emi_h2so4 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_emi_h2so4', 'H2SO4 emission flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_emi_part = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_emi_part', 'Particle emission flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_ocs_to_so2 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_ocs_to_so2', 'OCS to SO2 flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_so2_to_h2so4 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_so2_to_h2so4', 'SO2 to H2SO4 flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_h2so4_to_part = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_h2so4_to_part', 'H2SO4 to part flux', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_budg_sed_part = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'budg_sed_part', 'Ground sedimentation flux of strat. particles', 'kg(S)/m2/s', (/ ('', i=1, 10) /))
!--surface PM25 due to strat aerosol
  TYPE(ctrl_out), SAVE :: o_surf_PM25_sulf = ctrl_out((/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 1/), &
    'surf_PM25_sulf', 'Sulfate PM2.5 concentration at the surface', 'ug/m3', (/ ('', i=1, 10) /))
#endif

!!!!!!!!!!!!!!!!!!!!!! 3D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  TYPE(ctrl_out), SAVE :: o_ec550aer = ctrl_out((/ 2, 6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ec550aer', 'Extinction at 550nm', 'm^-1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lwcon = ctrl_out((/ 2, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'lwcon', 'Cloud liquid water content', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_iwcon = ctrl_out((/ 2, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'iwcon', 'Cloud ice water content', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_temp = ctrl_out((/ 2, 3, 4, 10, 10, 10, 11, 11, 11, 11/), &
    'temp', 'Air temperature', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_heat_volc = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'heat_volc', 'SW heating rate due to volcano', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cool_volc = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cool_volc', 'LW cooling rate due to volcano', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_theta = ctrl_out((/ 2, 3, 4, 10, 10, 10, 11, 11, 11, 11/), &
    'theta', 'Potential air temperature', 'K', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ovap = ctrl_out((/ 2, 3, 4, 10, 10, 10, 11, 11, 11, 11/), &
    'ovap', 'Specific humidity', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ovapinit = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ovapinit', 'Specific humidity (begin of timestep)', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_oliq = ctrl_out((/ 2, 3, 4, 10, 10, 10, 11, 11, 11, 11/), &
    'oliq', 'Liquid water', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ocond = ctrl_out((/ 2, 3, 4, 10, 10, 10, 11, 11, 11, 11/), &
    'ocond', 'Condensed water', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wvapp = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wvapp', '', '', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_geop = ctrl_out((/ 2, 3, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'geop', 'Geopotential height', 'm2/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_vitu = ctrl_out((/ 2, 3, 4, 6, 10, 10, 11, 11, 11, 11/), &
    'vitu', 'Zonal wind', 'm/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_vitv = ctrl_out((/ 2, 3, 4, 6, 10, 10, 11, 11, 11, 11/), &
    'vitv', 'Meridional wind', 'm/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_vitw = ctrl_out((/ 2, 3, 10, 6, 10, 10, 11, 11, 11, 11/), &
    'vitw', 'Vertical wind', 'Pa/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_pres = ctrl_out((/ 2, 3, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'pres', 'Air pressure', 'Pa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_paprs = ctrl_out((/ 2, 3, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'paprs', 'Air pressure Inter-Couches', 'Pa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_mass = ctrl_out((/ 2, 3, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'mass', 'Masse Couches', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_zfull = ctrl_out((/ 2, 3, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'zfull', 'Altitude of full pressure levels', 'm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_zhalf = ctrl_out((/ 2, 3, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'zhalf', 'Altitude of half pressure levels', 'm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rneb = ctrl_out((/ 2, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rneb', 'Cloud fraction', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rnebjn = ctrl_out((/ 2, 5, 10, 10, 10, 10, 11, 11,11, 11/), &      
    'rnebjn', 'Cloud fraction in day', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rnebcon = ctrl_out((/ 2, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rnebcon', 'Convective Cloud Fraction', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rnebls = ctrl_out((/ 2, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rnebls', 'LS Cloud fraction', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rneblsvol = ctrl_out((/ 2, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rneblsvol', 'LS Cloud fraction by volume', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rhum = ctrl_out((/ 2, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rhum', 'Relative humidity', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ozone = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ozone', 'Ozone mole fraction', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ozone_light = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ozone_daylight', 'Daylight ozone mole fraction', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_upwd = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'upwd', 'saturated updraft', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_epmax_diag = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'epmax', 'epmax en fn cape', 'su', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ep = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ep', 'ep', 'su', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_duphy = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'duphy', 'Physics du', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtphy = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtphy', 'Physics dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqphy = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqphy', 'Physics dQ', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqphy2d = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqphy2d', 'Physics dQ', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqlphy = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqlphy', 'Physics dQL', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqlphy2d = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqlphy2d', 'Physics dQL', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqsphy = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqsphy', 'Physics dQS', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqsphy2d = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqsphy2d', 'Physics dQS', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_pr_con_l = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'pr_con_l', 'Convective precipitation lic', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_pr_con_i = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'pr_con_i', 'Convective precipitation ice', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_pr_lsc_l = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'pr_lsc_l', 'Large scale precipitation lic', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_pr_lsc_i = ctrl_out((/ 2, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'pr_lsc_i', 'Large scale precipitation ice', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_re = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    're', 'Cloud droplet effective radius', 'um', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_fl = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'fl', 'Denominator of Cloud droplet effective radius', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_scdnc = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'scdnc', 'Cloud droplet number concentration', 'm-3', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_reffclws = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'reffclws', 'Stratiform Cloud Droplet Effective Radius (aerosol diags.)', 'm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_reffclwc = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'reffclwc', 'Convective Cloud Droplet Effective Radius (aerosol diags.)', 'm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lcc3d = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'lcc3d', 'Cloud liquid fraction', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lcc3dcon = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'lcc3dcon', 'Convective cloud liquid fraction', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lcc3dstra = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'lcc3dstra', 'Stratiform cloud liquid fraction', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_icc3dcon = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'icc3dcon', 'Convective cloud ice fraction', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_icc3dstra = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'icc3dstra', 'Stratiform cloud ice fraction', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldicemxrat = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldicemxrat', 'Cloud Ice Mixing Ratio', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cldwatmxrat = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cldwatmxrat', 'Cloud Water Mixing Ratio', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_solbnd = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'solbnd', 'Top-of-Atmosphere Solar Insolation for each band', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_stratomask = ctrl_out((/ 2,  6, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'stratomask', 'Stratospheric fraction', '1', (/ ('', i=1, 10) /))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_albe_srf     = (/ &
      ctrl_out((/ 3, 7, 10, 7, 10, 10, 11, 11, 11, 11/),'albe_ter', "Albedo VIS surf. "//clnsurf(1),"-", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 7, 10, 7, 10, 10, 11, 11, 11, 11/),'albe_lic', "Albedo VIS surf. "//clnsurf(2),"-", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 7, 10, 7, 10, 10, 11, 11, 11, 11/),'albe_oce', "Albedo VIS surf. "//clnsurf(3),"-", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 7, 10, 7, 10, 10, 11, 11, 11, 11/),'albe_sic', "Albedo VIS surf. "//clnsurf(4),"-", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_ages_srf     = (/ &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'ages_ter', "Snow age", "day", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'ages_lic', "Snow age", "day", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'ages_oce',"Snow age", "day", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'ages_sic',"Snow age", "day", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_snow_srf     = (/ &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'snow_ter', "Snow", "kg/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'snow_lic', "Snow", "kg/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'snow_oce',"Snow", "kg/m2", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 10, 10, 10, 10, 10, 11, 11, 11, 11/),'snow_sic',"Snow", "kg/m2", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_z0m_srf     = (/ &
      ctrl_out((/ 3, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'z0m_ter', "Surface roughness "//clnsurf(1),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'z0m_lic', "Surface roughness "//clnsurf(2),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'z0m_oce', "Surface roughness "//clnsurf(3),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'z0m_sic', "Surface roughness "//clnsurf(4),"m", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE, DIMENSION(4) :: o_z0h_srf     = (/ &
      ctrl_out((/ 3, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'z0h_ter', "Surface roughness "//clnsurf(1),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'z0h_lic', "Surface roughness "//clnsurf(2),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'z0h_oce', "Surface roughness "//clnsurf(3),"m", (/ ('', i=1, 10) /)), &
      ctrl_out((/ 3, 6, 10, 10, 10, 10, 11, 11, 11, 11/),'z0h_sic', "Surface roughness "//clnsurf(4),"m", (/ ('', i=1, 10) /)) /)

  TYPE(ctrl_out), SAVE :: o_alb1 = ctrl_out((/ 3, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'alb1', 'Surface VIS albedo', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_alb2 = ctrl_out((/ 3, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'alb2', 'Surface Near IR albedo', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_clwcon = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'clwcon', 'Convective Cloud Liquid water content', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_Mipsh = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'Mipsh', 'mass flux shed from adiab. ascents', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_Ma = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'Ma', 'undilute adiab updraft mass flux', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dnwd = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dnwd', 'saturated downdraft', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dnwd0 = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dnwd0', 'unsat. downdraft', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_mc = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'mc', 'Convective mass flux', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ftime_deepcv = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ftime_deepcv', 'Fraction of time deep convection Occurs', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ftime_con = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ftime_con', 'Fraction of time convection Occurs', ' ', (/ ('', i=1, 10) /))
!!jyg    'ftime_con', 'Fraction of time convection Occurs', ' ',                 &
!!jyg      (/ 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', & 
!!jyg         'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)', 'inst(X)' /))
  TYPE(ctrl_out), SAVE :: o_dtdyn = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtdyn', 'Dynamics dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqdyn = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqdyn', 'Dynamics dQ', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqdyn2d = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqdyn2d', 'Dynamics dQ', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqldyn = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqldyn', 'Dynamics dQL', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqldyn2d = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqldyn2d', 'Dynamics dQL', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqsdyn = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqsdyn', 'Dynamics dQS', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqsdyn2d = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqsdyn2d', 'Dynamics dQS', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dudyn = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dudyn', 'Dynamics dU', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dvdyn = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dvdyn', 'Dynamics dV', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtcon = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtcon', 'Convection dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ducon = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ducon', 'Convection du', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dvcon = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dvcon', 'Convection dv', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqcon = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqcon', 'Convection dQ', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqcon2d = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqcon2d', 'Convection dQ', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtwak = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtwak', 'Wake dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqwak = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqwak', 'Wake dQ', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqwak2d = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqwak2d', 'Wake dQ', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cv_gen = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'cv_gen', 'Cumulonimbus genesis', '1/(m2 s)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wake_h = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wake_h', 'wake_h', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wake_dens = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wake_dens', 'number of wakes per m2', '1/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wake_s = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wake_s', 'wake_s', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wake_deltat = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wake_deltat', 'wake_deltat', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wake_deltaq = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wake_deltaq', 'wake_deltaq', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wake_omg = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'wake_omg', 'wake_omg', 'Pa/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wdtrainA = ctrl_out((/ 4, 5, 10,  4, 10, 10, 11, 11, 11, 11 /), &
    'wdtrainA', 'precipitation from AA', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wdtrainS = ctrl_out((/ 4, 5, 10,  4, 10, 10, 11, 11, 11, 11 /), &
    'wdtrainS', 'precipitation from shedding of AA', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_wdtrainM = ctrl_out((/ 4, 5, 10,  4, 10, 10, 11, 11, 11, 11 /), &
    'wdtrainM', 'precipitation from mixture', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_Vprecip = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'Vprecip', 'precipitation vertical profile', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_qtaa = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'qtaa', 'specific total water in adiabatic ascents', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_clwaa = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'Clwaa', 'specific condensed water in adiabatic ascents', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ftd = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ftd', 'tend temp due aux descentes precip', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_fqd = ctrl_out((/ 4, 5, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'fqd', 'tend vap eau due aux descentes precip', '-', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtlsc = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtlsc', 'Condensation dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtlschr = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtlschr', 'Large-scale condensational heating rate', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqlsc = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqlsc', 'Condensation dQ', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqlsc2d = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqlsc2d', 'Condensation dQ', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_beta_prec = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'beta_prec', 'LS Conversion rate to prec', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtvdf = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtvdf', 'Boundary-layer dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtdis = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtdis', 'TKE dissipation dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqvdf = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqvdf', 'Boundary-layer dQ', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqvdf2d = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqvdf2d', 'Boundary-layer dQ', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dteva = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dteva', 'Reevaporation dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqeva = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqeva', 'Reevaporation dQ', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqeva2d = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqeva2d', 'Reevaporation dQ', '(kg/m2)/s', (/ ('', i=1, 10) /))

!!!!!!!!!!!!!!!! Specifique thermiques
  TYPE(ctrl_out), SAVE :: o_dqlscth = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqlscth', 'dQ therm.', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqlscth2d = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqlscth2d', 'dQ therm.', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqlscst = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqlscst', 'dQ strat.', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqlscst2d = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqlscst2d', 'dQ strat.', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtlscth = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtlscth', 'dQ therm.', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtlscst = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtlscst', 'dQ strat.', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_plulth = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'plulth', 'Rainfall therm.', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_plulst = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'plulst', 'Rainfall strat.', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lmaxth = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'lmaxth', "Upper level thermals", "", (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ptconvth = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ptconvth', 'POINTS CONVECTIFS therm.', ' ', (/ ('', i=1, 10) /))
!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE(ctrl_out), SAVE :: o_ptconv = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ptconv', 'POINTS CONVECTIFS', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ratqs = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ratqs', 'RATQS', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtthe = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtthe', 'Thermal dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_duthe = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'duthe', 'Thermal du', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dvthe = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dvthe', 'Thermal dv', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_f_th = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'f_th', 'Thermal plume mass flux', 'kg/(m2*s)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_e_th = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'e_th', 'Thermal plume entrainment', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_w_th = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'w_th', 'Thermal plume vertical velocity', 'm/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_lambda_th = ctrl_out((/ 10, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'lambda_th', 'Thermal plume vertical velocity', 'm/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ftime_th = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ftime_th', 'Fraction of time Shallow convection occurs', ' ', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_q_th = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'q_th', 'Thermal plume total humidity', 'kg/kg', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_a_th = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'a_th', "Thermal plume fraction", "", (/ ('', i=1, 10) /))

  TYPE(ctrl_out), SAVE :: o_cloudth_sth = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    's_th', "Thermal plume saturation deficit", "kg/kg", (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cloudth_senv = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    's_env', "Environment saturation deficit", "kg/kg", (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cloudth_sigmath = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sigma_th', "Thermal plume gauss variance", "kg/kg", (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_cloudth_sigmaenv = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'sigma_env', "Environment gauss variance", "kg/kg", (/ ('', i=1, 10) /))

  TYPE(ctrl_out), SAVE :: o_d_th = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'd_th', 'Thermal plume detrainment', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_f0_th = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'f0_th', 'Thermal closure mass flux', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_zmax_th = ctrl_out((/ 4,  4,  4,  5, 10, 10, 11, 11, 11, 11/), &
    'zmax_th', 'Thermal plume height', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqthe = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqthe', 'Thermal dQ', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqthe2d = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqthe2d', 'Thermal dQ', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtajs = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtajs', 'Dry adjust. dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqajs = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqajs', 'Dry adjust. dQ', '(kg/kg)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqajs2d = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqajs2d', 'Dry adjust. dQ', '(kg/m2)/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtswr = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtswr', 'SW radiation dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtsw0 = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtsw0', 'CS SW radiation dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtlwr = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtlwr', 'LW radiation dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtlw0 = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtlw0', 'CS LW radiation dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtec = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtec', 'Cinetic dissip dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_duvdf = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'duvdf', 'Boundary-layer dU', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dvvdf = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dvvdf', 'Boundary-layer dV', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_duoro = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'duoro', 'Orography dU', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dvoro = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dvoro', 'Orography dV', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dulif = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dulif', 'Orography dU', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dvlif = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dvlif', 'Orography dV', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_du_gwd_hines = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'du_gwd_hines', 'Hines GWD dU', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dv_gwd_hines = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dv_gwd_hines', 'Hines GWD dV', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_du_gwd_front = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'du_gwd_front', 'Fronts GWD dU', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dv_gwd_front = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dv_gwd_front', 'Fronts GWD dV', 'm/s2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_east_gwstress = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'east_gwstress', 'Eastward GW Stress', 'Pa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_west_gwstress = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'west_gwstress', 'Westward GW Stress', 'Pa', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtoro = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtoro', 'Orography dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dtlif = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dtlif', 'Orography dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dthin = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dthin', 'Hines GWD dT', 'K/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dqch4 = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dqch4', 'H2O due to CH4 oxidation & photolysis', '(kg/kg)/s', (/ ('', i=1, 10) /))

  type(ctrl_out), save:: o_du_gwd_rando &
       = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), 'du_gwd_rando', &
       "Random gravity waves dU/dt", "m/s2", (/ ('', i=1, 10) /))
  type(ctrl_out), save:: o_dv_gwd_rando &
       = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), 'dv_gwd_rando', &
       "Random gravity waves dV/dt", "m/s2", (/ ('', i=1, 10) /))
  type(ctrl_out), save:: o_ustr_gwd_hines &
       = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), 'ustr_gwd_hines', &
       "zonal wind stress Hines gravity waves", "Pa", (/ ('', i=1, 10) /))
  type(ctrl_out), save:: o_vstr_gwd_hines &
       = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), 'vstr_gwd_hines', &
       "meridional wind stress Hines gravity waves", "Pa", (/ ('', i=1, 10) /))
  type(ctrl_out), save:: o_ustr_gwd_front &
       = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), 'ustr_gwd_front', &
       "zonal wind stress fronts gravity waves", "Pa", (/ ('', i=1, 10) /))
  type(ctrl_out), save:: o_vstr_gwd_front &
       = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), 'vstr_gwd_front', &
       "meridional wind stress fronts gravity waves", "Pa", (/ ('', i=1, 10) /))
  type(ctrl_out), save:: o_ustr_gwd_rando &
       = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), 'ustr_gwd_rando', &
       "zonal wind stress random gravity waves", "Pa", (/ ('', i=1, 10) /))
  type(ctrl_out), save:: o_vstr_gwd_rando &
       = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), 'vstr_gwd_rando', &
       "meridional wind stress random gravity waves", "Pa", (/ ('', i=1, 10) /))

  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_trac(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_trac_cum(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_vdf(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_the(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_con(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_lessi_impa(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_lessi_nucl(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_insc(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_bcscav(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_evapls(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_ls(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_trsp(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_sscav(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_sat(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_uscav(:)
  TYPE(ctrl_out), SAVE, ALLOCATABLE :: o_dtr_dry(:)

  TYPE(ctrl_out), SAVE :: o_rsu = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsu', 'SW upward radiation', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsd = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsd', 'SW downward radiation', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rlu = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rlu', 'LW upward radiation', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rld = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rld', 'LW downward radiation', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsucs = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsucs', 'SW CS upward radiation', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsucsaf = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsucsaf', 'SW CS clean (no aerosol) upward radiation', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsdcs = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsdcs', 'SW CS downward radiation', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsdcsaf = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsdcsaf', 'SW CS clean (no aerosol) downward radiation', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rlucs = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rlucs', 'LW CS upward radiation', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rldcs = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rldcs', 'LW CS downward radiation', 'W m-2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tnt = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tnt', 'Tendency of air temperature', 'K s-1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tntc = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tntc', 'Tendency of air temperature due to Moist Convection', 'K s-1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tntr = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tntr', 'Air temperature tendency due to Radiative heating', 'K s-1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tntscpbl = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/),                  &
    'tntscpbl', 'Air temperature tendency due to St cloud and precipitation and BL mixing', &
      'K s-1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tnhus = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tnhus', 'Tendency of specific humidity', 's-1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tnhusc = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tnhusc', 'Tendency of specific humidity due to convection', 's-1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_tnhusscpbl = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'tnhusscpbl', 'Tendency of Specific humidity due to ST cl, precip and BL mixing', 's-1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_evu = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'evu', 'Eddy viscosity coefficient for Momentum Variables', 'm2 s-1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_h2o = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'h2o', 'Mass Fraction of Water', '1', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_mcd = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'mcd', 'Downdraft COnvective Mass Flux', 'kg/(m2*s)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_dmc = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'dmc', 'Deep COnvective Mass Flux', 'kg/(m2*s)', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ref_liq = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ref_liq', 'Effective radius of convective cloud liquid water particle', 'm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_ref_ice = ctrl_out((/ 4, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'ref_ice', 'Effective radius of startiform cloud ice particle', 'm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsut4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsut4co2', 'TOA Out SW in 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rlut4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rlut4co2', 'TOA Out LW in 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsutcs4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsutcs4co2', 'TOA Out CS SW in 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rlutcs4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rlutcs4co2', 'TOA Out CS LW in 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsu4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsu4co2', 'Upwelling SW 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rlu4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rlu4co2', 'Upwelling LW 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsucs4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsucs4co2', 'Upwelling CS SW 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rlucs4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rlucs4co2', 'Upwelling CS LW 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsd4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsd4co2', 'Downwelling SW 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rld4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rld4co2', 'Downwelling LW 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rsdcs4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rsdcs4co2', 'Downwelling CS SW 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_rldcs4co2 = ctrl_out((/ 5, 10, 10, 10, 10, 10, 11, 11, 11, 11/), &
    'rldcs4co2', 'Downwelling CS LW 4xCO2 atmosphere', 'W/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_snowsrf = ctrl_out((/ 1, 1, 10, 1, 10, 10, 11, 11, 11, 11/), &
    'snowsrf', 'Snow mass at surface', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_qsnow = ctrl_out((/ 1, 1, 10, 1, 10, 10, 11, 11, 11, 11/), &
    'qsnow', 'Water contained in snow', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_snowhgt = ctrl_out((/ 1, 1, 10, 1, 10, 10, 11, 11, 11, 11/), &
    'snowhgt', 'Snow height at surface', 'm', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_toice = ctrl_out((/ 1, 1, 10, 1, 10, 10, 11, 11, 11, 11/), &
    'to_ice', 'Snow passed to ice model', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_sissnow = ctrl_out((/ 1, 1, 10, 1, 10, 10, 11, 11, 11, 11/), &
    'sissnow', 'Snow in snow model', 'kg/m2', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_runoff = ctrl_out((/ 1, 1, 10, 1, 10, 10, 11, 11, 11, 11/), &
    'runoff', 'Run-off rate land ice', 'kg/m2/s', (/ ('', i=1, 10) /))
  TYPE(ctrl_out), SAVE :: o_albslw3 = ctrl_out((/ 1, 1, 1, 1, 10, 10, 11, 11, 11, 11/), &
    'albslw3', 'Surface albedo LW3', '-', (/ ('', i=1, 10) /))

!!!!!!!!!!!!! Sorties niveaux standards de pression NMC 
  TYPE(ctrl_out), SAVE :: o_tnondef = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'tnondef', 'Undefined value of T', 'K', (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_ta = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'ta', 'Air temperature', 'K', (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_zg  = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'zg', 'Geopotential height', 'm', (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_hus = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'hus', 'Specific humidity', '1', (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_hur = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'hur', 'Relative humidity', '%', (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_ua = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'ua', 'Eastward wind', 'm s-1', (/ "inst(X)", "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_va = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'va', 'Northward wind', 'm s-1', (/ ('', i=1, 10)/))
  TYPE(ctrl_out), SAVE :: o_wap = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'wap', 'Lagrangian tendency of air pressure', 'Pa s-1', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_psbg = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'psbg', 'Pressure sfce below ground', '%', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_tro3 = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'tro3', 'Ozone mole fraction', '1e-9', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_tro3_daylight = ctrl_out((/ 11, 11, 11, 11, 11, 11, 5, 5, 5, 11/), &
       'tro3_daylight', 'Daylight ozone mole fraction', '1e-9', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_uxv = ctrl_out((/ 11, 11, 11, 11, 11, 11, 6, 6, 6, 11/), &
       'uv', 'uv', 'm2/s2', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_vxq = ctrl_out((/ 11, 11, 11, 11, 11, 11, 6, 6, 6, 11/), &
       'vxq', 'vxq', 'm/s * (kg/kg)', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_vxT = ctrl_out((/ 11, 11, 11, 11, 11, 11, 6, 6, 6, 11/), &
       'vT', 'vT', 'mK/s', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_wxq = ctrl_out((/ 11, 11, 11, 11, 11, 11, 6, 6, 6, 11/), &
       'wq', 'wq', '(Pa/s)*(kg/kg)', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_vxphi = ctrl_out((/ 11, 11, 11, 11, 11, 11, 6, 6, 6, 11/), &
       'vphi', 'vphi', 'm2/s', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_wxT = ctrl_out((/ 11, 11, 11, 11, 11, 11, 6, 6, 6, 11/), &
       'wT', 'wT', '"K*Pa/s', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /)) 
  TYPE(ctrl_out), SAVE :: o_uxu = ctrl_out((/ 11, 11, 11, 11, 11, 11, 6, 6, 6, 11/), &
       'u2', 'u2', 'm2/s2', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
  TYPE(ctrl_out), SAVE :: o_vxv = ctrl_out((/ 11, 11, 11, 11, 11, 11, 6, 6, 6, 11/), &
       'v2', 'v2', 'm2/s2', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))
   TYPE(ctrl_out), SAVE :: o_TxT = ctrl_out((/ 11, 11, 11, 11, 11, 11, 6, 6, 6, 11/), &
       'T2', 'T2', 'K2', (/ "inst(X)", "inst(X)", "inst(X)", &
       "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)", "inst(X)"  /))

#ifdef CPP_Dust
#include "Dust/spla_output_dat.h"
#endif

END MODULE phys_output_ctrlout_mod
