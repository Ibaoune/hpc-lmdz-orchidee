      module Mod_SISVAT_dzS

      use    Mod_Real


      implicit none


      integer, SAVE          ,ALLOCATABLE ,dimension(:)   ::  islpSV
      integer, SAVE          ,ALLOCATABLE ,dimension(:)   ::  isnpSV
      integer, SAVE          ,ALLOCATABLE ,dimension(:)   ::  islmSV

      real(kind=real8), SAVE                              ::  Implic,Explic
      real(kind=real8), SAVE                              ::  OcndSV        ! Swab Ocean / Soil Ratio 
      real(kind=real8), SAVE ,ALLOCATABLE ,dimension(:)   ::  dzmiSV        ! dz_(i-1/2)
      real(kind=real8), SAVE ,ALLOCATABLE ,dimension(:)   ::  dzi_SV        ! dz_(i-1)/(dz_(i)+dz_(i-1))
      real(kind=real8), SAVE ,ALLOCATABLE ,dimension(:)   ::  dziiSV        ! dz_(i)  /(dz_(i)+dz_(i-1))
      real(kind=real8), SAVE ,ALLOCATABLE ,dimension(:)   ::  dtz_SV        ! dt / dz
      real(kind=real8), SAVE ,ALLOCATABLE ,dimension(:)   ::  dz78SV        ! 7/8 (dz)
      real(kind=real8), SAVE ,ALLOCATABLE ,dimension(:)   ::  dz34SV        ! 3/4 (dz)
      real(kind=real8), SAVE ,ALLOCATABLE ,dimension(:)   ::  dz_8SV        ! 1/8 (dz)
      real(kind=real8), SAVE ,ALLOCATABLE ,dimension(:)   ::  dzAvSV        ! 1/8dz_(-1)+3/4dz+1/8dz_(+1)
      real(kind=real8), SAVE ,ALLOCATABLE ,dimension(:,:) ::  RF__SV        ! 1/8dz_(-1)+3/4dz+1/8dz_(+1)

      end module Mod_SISVAT_dzS
