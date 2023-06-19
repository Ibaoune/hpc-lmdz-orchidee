SUBROUTINE inistats(ierr)

  USE vertical_layers_mod, ONLY: ap,bp,preff,presnivs,pseudoalt

  IMPLICIT NONE

  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
  include "comconst.h"
  include "statto.h"
  include "netcdf.inc"

  INTEGER, INTENT (OUT) :: ierr
  INTEGER :: nid
  INTEGER :: l, nsteppd
  REAL, DIMENSION (llm) :: sig_s
  INTEGER :: idim_lat, idim_lon, idim_llm, idim_llmp1, idim_time
  REAL, DIMENSION (istime) :: lt
  INTEGER :: nvarid

  WRITE (*, *)
  WRITE (*, *) '                        || STATS ||'
  WRITE (*, *)
  WRITE (*, *) 'daysec', daysec
  WRITE (*, *) 'dtphys', dtphys
  nsteppd = nint(daysec/dtphys)
  WRITE (*, *) 'nsteppd=', nsteppd
  IF (abs(float(nsteppd)-daysec/dtphys)>1.E-8*daysec) &
    STOP 'Dans Instat:  1jour .ne. n pas physiques'

  IF (mod(nsteppd,istime)/=0) STOP &
    'Dans Instat:  1jour .ne. n*istime pas physiques'

  istats = nsteppd/istime
  WRITE (*, *) 'istats=', istats
  WRITE (*, *) 'Storing ', istime, 'times per day'
  WRITE (*, *) 'thus every ', istats, 'physical timestep '
  WRITE (*, *)

  DO l = 1, llm
    sig_s(l) = ((ap(l)+ap(l+1))/preff+bp(l)+bp(l+1))/2.
  END DO

  ierr = nf_create('stats.nc', nf_clobber, nid)
  IF (ierr/=nf_noerr) THEN
    WRITE (*, *) nf_strerror(ierr)
    STOP ''
  END IF

  ierr = nf_def_dim(nid, 'latitude', jjp1, idim_lat)
  ierr = nf_def_dim(nid, 'longitude', iip1, idim_lon)
  ierr = nf_def_dim(nid, 'altitude', llm, idim_llm)
  ierr = nf_def_dim(nid, 'llmp1', llm+1, idim_llmp1)
  ierr = nf_def_dim(nid, 'Time', nf_unlimited, idim_time)

  ierr = nf_enddef(nid)
  CALL def_var_stats(nid, 'Time', 'Time', 'hours since 0000-00-0 00:00:00', &
    1, idim_time, nvarid, ierr)
  ! Time is initialised later by mkstats subroutine

  CALL def_var_stats(nid, 'latitude', 'latitude', 'degrees_north', 1, &
    idim_lat, nvarid, ierr)
#ifdef NC_DOUBLE
  ierr = nf_put_var_double(nid, nvarid, rlatu/pi*180)
#else
  ierr = nf_put_var_real(nid, nvarid, rlatu/pi*180)
#endif
  CALL def_var_stats(nid, 'longitude', 'East longitude', 'degrees_east', 1, &
    idim_lon, nvarid, ierr)
#ifdef NC_DOUBLE
  ierr = nf_put_var_double(nid, nvarid, rlonv/pi*180)
#else
  ierr = nf_put_var_real(nid, nvarid, rlonv/pi*180)
#endif

  ! Niveaux verticaux, aps et bps
  ierr = nf_redef(nid)
  ! presnivs
#ifdef NC_DOUBLE
  ierr = nf_def_var(nid, 'presnivs', nf_double, 1, idim_llm, nvarid)
#else
  ierr = nf_def_var(nid, 'presnivs', nf_float, 1, idim_llm, nvarid)
#endif
  ierr = nf_put_att_text(nid, nvarid, 'long_name', 15, 'Vertical levels')
  ierr = nf_put_att_text(nid, nvarid, 'units', 2, 'Pa')
  ierr = nf_put_att_text(nid, nvarid, 'positive', 4, 'down')
  ierr = nf_enddef(nid)
#ifdef NC_DOUBLE
  ierr = nf_put_var_double(nid, nvarid, presnivs(1:llm))
#else
  ierr = nf_put_var_real(nid, nvarid, presnivs(1:llm))
#endif
  ! Pseudo alts
#ifdef NC_DOUBLE
  ierr = nf_def_var(nid, 'altitude', nf_double, 1, idim_llm, nvarid)
#else
  ierr = nf_def_var(nid, 'altitude', nf_float, 1, idim_llm, nvarid)
#endif
  ierr = nf_put_att_text(nid, nvarid, 'long_name', 8, 'altitude')
  ierr = nf_put_att_text(nid, nvarid, 'units', 2, 'km')
  ierr = nf_put_att_text(nid, nvarid, 'positive', 2, 'up')
  ierr = nf_enddef(nid)
#ifdef NC_DOUBLE
  ierr = nf_put_var_double(nid, nvarid, pseudoalt)
#else
  ierr = nf_put_var_real(nid, nvarid, pseudoalt)
#endif
  ! call def_var_stats(nid,"aps","hybrid pressure at midlayers"," ",
  ! &            1,idim_llm,nvarid,ierr)
  ! #ifdef NC_DOUBLE
  ! ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aps)
  ! #else
  ! ierr = NF_PUT_VAR_REAL (nid,nvarid,aps)
  ! #endif

  ! call def_var_stats(nid,"bps","hybrid sigma at midlayers"," ",
  ! &            1,idim_llm,nvarid,ierr)
  ! #ifdef NC_DOUBLE
  ! ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bps)
  ! #else
  ! ierr = NF_PUT_VAR_REAL (nid,nvarid,bps)
  ! #endif

  ierr = nf_close(nid)

END SUBROUTINE inistats

