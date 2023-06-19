!
! $Id $
!
      COMMON/coefils/jfiltnu,jfiltsu,jfiltnv,jfiltsv,sddu(iim),sddv(iim)&
     & ,unsddu(iim),unsddv(iim),coefilu(iim,jjm),coefilv(iim,jjm),      &
     & modfrstu(jjm),modfrstv(jjm),eignfnu(iim,iim),eignfnv(iim,iim)    &
     & ,coefilu2(iim,jjm),coefilv2(iim,jjm)
!c
      INTEGER jfiltnu ! index of the last lat line filtered in NH (U grid)
      INTEGER jfiltsu ! index of the first lat line filtered in SH (U grid)
      INTEGER jfiltnv ! index of the last lat line filtered in NH (V grid)
      INTEGER jfiltsv ! index of the first lat line filtered in SH (V grid)
      INTEGER modfrstu ! number of retained (ie: unfiltered) modes on U grid
      INTEGER modfrstv ! number of retained (ie: unfiltered) modes on V grid
      REAL    sddu,sddv,unsddu,unsddv,coefilu,coefilv,eignfnu,eignfnv
      REAL    coefilu2,coefilv2
