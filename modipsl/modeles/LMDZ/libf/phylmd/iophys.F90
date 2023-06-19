!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface pour ecrire en netcdf avec les routines d'enseignement
! iotd de Frederic Hourdin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine iophys_ecrit(nom,lllm,titre,unite,px)

      USE mod_phys_lmdz_para, ONLY: klon_omp, is_mpi_root
      USE mod_phys_lmdz_transfert_para, ONLY: gather
      USE mod_grid_phy_lmdz, ONLY: klon_glo, nbp_lon, nbp_lat, grid1dto2d_glo
      IMPLICIT NONE



!  Ecriture de variables diagnostiques au choix dans la physique 
!  dans un fichier NetCDF nomme  'diagfi'. Ces variables peuvent etre
!  3d (ex : temperature), 2d (ex : temperature de surface), ou
!  0d (pour un scalaire qui ne depend que du temps : ex : la longitude
!  solaire)
!  (ou encore 1d, dans le cas de testphys1d, pour sortir une colonne)
!  La periode d'ecriture est donnee par 
!  "ecritphy " regle dans le fichier de controle de run :  run.def
!
!    writediagfi peut etre appele de n'importe quelle subroutine
!    de la physique, plusieurs fois. L'initialisation et la creation du
!    fichier se fait au tout premier appel.
!
! WARNING : les variables dynamique (u,v,t,q,ps)
!  sauvees par writediagfi avec une
! date donnee sont legerement differentes que dans le fichier histoire car 
! on ne leur a pas encore ajoute de la dissipation et de la physique !!!
! IL est  RECOMMANDE d'ajouter les tendance physique a ces variables
! avant l'ecriture dans diagfi (cf. physiq.F)
!  
! Modifs: Aug.2010 Ehouarn: enforce outputs to be real*4
!
!  parametres (input) :
!  ----------
!      unit : unite logique du fichier de sortie (toujours la meme)
!      nom  : nom de la variable a sortir (chaine de caracteres)
!      titre: titre de la variable (chaine de caracteres)
!      unite : unite de la variable (chaine de caracteres)
!      px : variable a sortir (real 0, 1, 2, ou 3d)
!
!=================================================================


! Arguments on input:
      integer lllm
      character (len=*) :: nom,titre,unite
      integer imjmax
      parameter (imjmax=100000)
      real px(klon_omp,lllm)
      real xglo(klon_glo,lllm)
      real zx(nbp_lon,nbp_lat,lllm)


      CALL Gather(px,xglo)
!$OMP MASTER
      IF (is_mpi_root) THEN       
        CALL Grid1Dto2D_glo(xglo,zx)
        call iotd_ecrit(nom,lllm,titre,unite,zx)
      ENDIF
!$OMP END MASTER

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Version avec reindexation pour appeler depuis les routines internes
! à la sous surface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine iophys_ecrit_index(nom,lllm,titre,unite,knon,knindex,px)

    USE mod_phys_lmdz_para, ONLY: klon_omp
    USE dimphy, ONLY : klon
    USE mod_grid_phy_lmdz, ONLY: klon_glo
    IMPLICIT NONE

! This subroutine returns the sea surface temperature already read from limit.nc
!

! Arguments on input:
    INTEGER lllm
    CHARACTER (len=*) :: nom,titre,unite
    REAL px(klon_omp,lllm)
    INTEGER, INTENT(IN)                  :: knon     ! nomber of points on compressed grid
    INTEGER, DIMENSION(klon), INTENT(IN) :: knindex  ! grid point number for compressed grid
    REAL, DIMENSION(klon,lllm) :: varout

    INTEGER :: i,l

    IF (klon/=klon_omp) THEN
      print*,'klon, klon_omp',klon,klon_omp
      CALL abort_physic('iophys_ecrit','probleme de dimension parallele',1)
    ENDIF

    varout(1:klon,1:lllm)=0.
    DO l = 1, lllm
    DO i = 1, knon
       varout(knindex(i),l) = px(i,l)
    END DO
    END DO
    CALL iophys_ecrit(nom,lllm,titre,unite,varout)

  END SUBROUTINE iophys_ecrit_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE iophys_ini
      USE mod_phys_lmdz_para, ONLY: is_mpi_root
      USE vertical_layers_mod, ONLY: presnivs
      USE regular_lonlat_mod, ONLY: lon_reg, lat_reg
      USE dimphy, ONLY: klev
      USE mod_grid_phy_lmdz, ONLY: klon_glo

      IMPLICIT NONE

!=======================================================================
!
!   Auteur:  L. Fairhead  ,  P. Le Van, Y. Wanherdrick, F. Forget
!   -------
!
!   Objet:
!   ------
!
!   'Initialize' the diagfi.nc file: write down dimensions as well
!   as time-independent fields (e.g: geopotential, mesh area, ...)
!
!=======================================================================
!-----------------------------------------------------------------------
!   Declarations:
!   -------------

real pi
INTEGER nlat_eff

!   Arguments:
!   ----------

!$OMP MASTER
    IF (is_mpi_root) THEN       

! Bidouille pour gerer le fait que lat_reg contient deux latitudes
! en version uni-dimensionnelle (chose qui pourrait être résolue
! par ailleurs)
IF (klon_glo==1) THEN
   nlat_eff=1
ELSE
   nlat_eff=size(lat_reg)
ENDIF
pi=2.*asin(1.)
call iotd_ini('phys.nc   ', &
size(lon_reg),nlat_eff,klev,lon_reg(:)*180./pi,lat_reg*180./pi,presnivs)
    ENDIF
!$OMP END MASTER

      END

#ifdef und
      SUBROUTINE gr_fi_ecrit(nfield,nlon,iim,jjmp1,fi,ecrit)
      IMPLICIT none

!=======================================================================
      INTEGER nfield,nlon,iim,jjmp1, jjm
      REAL fi(nlon,nfield), ecrit(iim*jjmp1,nfield)

      INTEGER i, n, ig

      jjm = jjmp1 - 1
      DO n = 1, nfield
         DO i=1,iim
            ecrit(i,n) = fi(1,n)
            ecrit(i+jjm*iim,n) = fi(nlon,n)
         ENDDO
         DO ig = 1, nlon - 2
           ecrit(iim+ig,n) = fi(1+ig,n)
         ENDDO
      ENDDO
      RETURN
      END

#endif
