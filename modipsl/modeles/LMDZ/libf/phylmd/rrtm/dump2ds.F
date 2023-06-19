      SUBROUTINE dump2ds(im,jm,z,nom_z)
C Copyright (C) 2005 Centre National de la Recherche Scientifique
c ==================================================================
c Perform a scatter plot print of big matrices using regular
c intervals between min and max matrix coefficient values.
c ==================================================================
c adapted from LMD3 by Alain Lahellec and retranscipted for LMDZ5
c NAN et INF ajoute aux plots                           Pat fin 2006
c ==================================================================
c Comme dump2d sauf que le signe est préservé, la valeur zero
c identifiee par un blanc.
c detection des Infty (= ou -) et NaN (?)
c ==================================================================
      IMPLICIT NONE
      INTEGER im,jm
      REAL z(im,jm),az
      CHARACTER*32 jform,jline*1000
      CHARACTER*1 iform(32)
      CHARACTER*16 F1000
      logical zinf,znan,zsign
      CHARACTER (len=*) :: nom_z
      INTEGER*4 icheck(2)
      EQUIVALENCE(az,icheck)
      EQUIVALENCE(iform,jform)
      DATA jform/'-@zyxwvutsrqpon NOPQRSTUVWXYZ*+?'/
      INTEGER i,j,k,imin,imax,jmin,jmax,kzero,kchar(im)
      REAL zmin,zmax,zllu,zllm
      write(F1000,'(''(4x,'',I3,''(1H-))'')')im+3
      DO 10001 i=1,200
      jline(1+(i-1)*5:5*i)='.    '
10001 CONTINUE
10002 zmin=z(1,1)
      imin=1
      jmin=1
      zmax=z(1,1)
      imax=1
      jmax=1
      kzero=0
      DO 10003 j=1,jm
      DO 10005 i=1,im
      IF(.NOT.( z(i,j).GT.zmax))GOTO 10007
      zmax=z(i,j)
      imax=i
      jmax=j
10007 IF(.NOT.( z(i,j).LT.zmin))GOTO 10009
      zmin=z(i,j)
      imin=i
      jmin=j
10009 IF(.NOT.( z(i,j).eq.0.))GOTO 10011
      kzero=kzero+1
10011 CONTINUE
10005 CONTINUE
10006 CONTINUE
10003 CONTINUE
10004 zsign=(sign(1.,zmin)*sign(1.,zmax).gt.0.)
      WRITE(*,*)'>>> dump2ds: ',trim(nom_z)
      PRINT*,'>>>  ',kzero,' zero values  <<<'
      IF(.NOT.( ZMin.lt.0.))GOTO 10013
      PRINT2000,' [MIN:',zmin,'(',imin,',',jmin,') (-@zyxwvutsrqpon NOPQ
     *RSTUVWXYZ*+) ',zmax,'(',imax,',',jmax,'):MAX]'
      GOTO 10014
10013 PRINT2000,' [MIN:',zmin,'(',imin,',',jmin,') ( NOPQRSTUVWXYZ*+) ',
     *zmax,'(',imax,',',jmax,'):MAX]'
10014 CONTINUE
2000  Format(a,1pg11.4,a1,i3,a1,i3,a,1pg11.4,a1,i3,a1,i3,a)
      IF(.NOT.( zmax.GT.zmin))GOTO 10015
      zllm=max(abs(zmax),abs(zmin))
      zllu=min(abs(zmax),abs(zmin))
      if(im.ge.100)WRITE(*,'(104x,900i1)')(mod(i/100,10),i=100,im)
      WRITE(*,'(14x,1000i1)')(mod(i/10,10),i=10,im)
      WRITE(*,'(5x,1009i1)')(mod(i,10),i=1,im)
      write(*,F1000)
      zinf=.false.
      znan=.false.
      DO 10017 j=1,jm
      DO 10019 i=1,im
      az=abs(z(i,j))
      IF(.NOT.( az.eq.0.))GOTO 10021
      kchar(i)=16
      GOTO 10022
10021 IF(.NOT.((az.ne.0.and.icheck(1).eq.0.and.icheck(2).eq.2146435072))
     *)GOTO 10023
      kchar(i)=31
      zinf=.true.
      GOTO 10022
10023 IF(.NOT.((az.ne.0.and.icheck(1).eq.0.and.icheck(2).eq.2146959360))
     *)GOTO 10024
      kchar(i)=32
      znan=.true.
      GOTO 10022
10024 IF(.NOT.( zsign))GOTO 10025
      kchar(i)=NINT(13.*(az-zllu)/(zllm-zllu)+17)
      GOTO 10026
10025 kchar(i)=NINT(13.*az/zllm+17)
10026 CONTINUE
10022 IF(.NOT.( z(i,j).lt.0.))GOTO 10027
      kchar(i)=32-kchar(i)
10027 CONTINUE
10019 CONTINUE
10020 WRITE(*,'(1x,i3,''|'',1000a)')j+1,(iform(kchar(i)),i=1,im),'|','| 
     *'
10017 CONTINUE
10018 write(*,F1000)
      WRITE(*,'(5x,1000i1)')(mod(i,10),i=1,im)
      WRITE(*,'(14x,1000i1)')(mod(i/10,10),i=10,im)
      if(im.ge.100)WRITE(*,'(104x,900i1)')(mod(i/100,10),i=100,im)
      GOTO 10016
10015 print*,'>>> ZERO MAP  <<<'
10016 if(zinf)print*,' *** Infty value(s) (+ or -) in map ***'
      if(znan)print*,' *** NaN value(s) (?) in map ***'
      print*
      RETURN
      END
