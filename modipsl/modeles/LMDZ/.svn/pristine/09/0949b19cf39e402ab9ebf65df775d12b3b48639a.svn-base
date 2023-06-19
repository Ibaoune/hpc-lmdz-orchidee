!IM for NMC files
!      real twriteSTD(klon,nlevSTD,nfiles)
!      real qwriteSTD(klon,nlevSTD,nfiles)
!      real rhwriteSTD(klon,nlevSTD,nfiles)
!      real phiwriteSTD(klon,nlevSTD,nfiles)
!      real uwriteSTD(klon,nlevSTD,nfiles)
!      real vwriteSTD(klon,nlevSTD,nfiles)
!      real wwriteSTD(klon,nlevSTD,nfiles)

      real twriteSTD3(klon,nlevSTD3)
      real qwriteSTD3(klon,nlevSTD3)
      real rhwriteSTD3(klon,nlevSTD3)
      real phiwriteSTD3(klon,nlevSTD3)
      real uwriteSTD3(klon,nlevSTD3)
      real vwriteSTD3(klon,nlevSTD3)
      real wwriteSTD3(klon,nlevSTD3)

      real tnondefSTD8(klon,nlevSTD8)
      real twriteSTD8(klon,nlevSTD8)
      real qwriteSTD8(klon,nlevSTD8)
      real rhwriteSTD8(klon,nlevSTD8)
      real phiwriteSTD8(klon,nlevSTD8)
      real uwriteSTD8(klon,nlevSTD8)
      real vwriteSTD8(klon,nlevSTD8)
      real wwriteSTD8(klon,nlevSTD8)

      real, save :: rlevSTD(nlevSTD)
      DATA rlevSTD/100000., 92500., 85000., 70000., &
      60000., 50000., 40000., 30000., 25000., 20000., &
      15000., 10000., 7000., 5000., 3000., 2000., 1000./
!$OMP THREADPRIVATE(rlevstd)

      CHARACTER*4, SAVE :: clevSTD(nlevSTD)
      DATA clevSTD/'1000','925 ','850 ','700 ','600 ', &
      '500 ','400 ','300 ','250 ','200 ','150 ','100 ', &
      '70  ','50  ','30  ','20  ','10  '/
!$OMP THREADPRIVATE(clevSTD)

      real, save :: rlevSTD3(nlevSTD3)
      DATA rlevSTD3/85000., 50000., 25000./
!$OMP THREADPRIVATE(rlevSTD3)

      real, save :: rlevSTD8(nlevSTD8)
      DATA rlevSTD8/100000., 85000., 70000., 50000., 25000., 10000., &
           5000., 1000./
!$OMP THREADPRIVATE(rlevSTD8) 
!
      REAL geo500(klon)

! nout : niveau de output des variables a une pression donnee
      logical oknondef(klon,nlevSTD,nout)
!
! les produits uvSTD, vqSTD, .., T2STD sont calcules
! a partir des valeurs instantannees toutes les 6 h
! qui sont moyennees sur le mois

      REAL zx_tmp_fiNC(klon,nlevSTD)

!     REAL missing_val
      REAL, SAVE :: freq_moyNMC(nout)
!$OMP THREADPRIVATE(freq_moyNMC)
