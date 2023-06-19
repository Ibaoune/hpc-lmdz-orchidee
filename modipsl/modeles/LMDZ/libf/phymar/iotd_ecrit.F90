      subroutine iotd_ecrit(nom,llm,titre,unite,px)


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
 
      implicit none

! Commons

#include "netcdf.inc"
#include "iotd.h"


! Arguments on input:
      integer llm
      character (len=*) :: nom,titre,unite
      integer imjmax
      parameter (imjmax=100000)
      real px(imjmax*llm)

! Local variables:

      real*4 date
      real*4 zx(imjmax*llm)


      integer ierr,ndim,dim_cc(4)
      integer iq
      integer i,j,l

      integer zitau
      character firstnom*20
      SAVE firstnom
      SAVE zitau
      SAVE date
      data firstnom /'1234567890'/
      data zitau /0/

! Ajouts
      integer, save :: ntime=0
      integer :: idim,varid
      character (len =50):: fichnom
      integer, dimension(4) :: id
      integer, dimension(4) :: edges,corner
      

!***************************************************************
! Initialisation of 'firstnom' and create/open the "diagfi.nc" NetCDF file
! ------------------------------------------------------------------------
! (Au tout premier appel de la subroutine durant le run.)


!--------------------------------------------------------
! Write the variables to output file if it's time to do so
!--------------------------------------------------------


! Compute/write/extend 'Time' coordinate (date given in days)
! (done every "first call" (at given time level) to writediagfi)
! Note: date is incremented as 1 step ahead of physics time
!--------------------------------------------------------

        zx(1:imax*jmax*llm)=px(1:imax*jmax*llm)
        if (firstnom =='1234567890') then
            firstnom=nom
        endif

!      print*,'nom ',nom,firstnom

!! Quand on tombe sur la premiere variable on ajoute un pas de temps
        if (nom.eq.firstnom) then
        ! We have identified a "first call" (at given date)

           ntime=ntime+1 ! increment # of stored time steps

!!          print*,'ntime ',ntime
           date=ntime
!          date= float (zitau +1)/float (day_step)

           ! compute corresponding date (in days and fractions thereof)
           ! Get NetCDF ID of 'Time' variable

           ierr=NF_SYNC(nid)

           ierr= NF_INQ_VARID(nid,"Time",varid)
           ! Write (append) the new date to the 'Time' array


           ierr= NF_PUT_VARA_REAL(nid,varid,ntime,1,date)

!          print*,'date ',date,ierr,nid
!        print*,'IOTD Date ,varid,nid,ntime,date',varid,nid,ntime,date

           if (ierr.ne.NF_NOERR) then
              write(*,*) "***** PUT_VAR matter in writediagfi_nc"
              write(*,*) "***** with time"
              write(*,*) 'ierr=', ierr   
           endif

!          write(6,*)'WRITEDIAGFI: date= ', date
        end if ! of if (nom.eq.firstnom)


!Case of a 3D variable
!---------------------
        if (llm==lmax) then
           ndim=4
           corner(1)=1
           corner(2)=1
           corner(3)=1
           corner(4)=ntime
           edges(1)=imax
           edges(2)=jmax
           edges(3)=llm
           edges(4)=1
           dim_cc=dim_coord


!Case of a 2D variable
!---------------------

        else if (llm==1) then
           ndim=3
           corner(1)=1
           corner(2)=1
           corner(3)=ntime
           corner(4)=1
           edges(1)=imax
           edges(2)=jmax
           edges(3)=1
           edges(4)=1
           dim_cc(1:2)=dim_coord(1:2)
           dim_cc(3)=dim_coord(4)

        endif ! of if llm=1 ou llm

! AU premier pas de temps, on crée les variables
!-----------------------------------------------

      if (ntime==1) then
          ierr = NF_REDEF (nid)
          ierr = NF_DEF_VAR(nid,nom,NF_FLOAT,ndim,dim_cc,varid)
          print*,'DEF ',nom,nid,varid
          ierr = NF_ENDDEF(nid)
      else
         ierr= NF_INQ_VARID(nid,nom,varid)
          print*,'INQ ',nom,nid,varid
! Commandes pour recuperer automatiquement les coordonnees
!             ierr= NF_INQ_DIMID(nid,"longitude",id(1))
      endif


      ierr= NF_PUT_VARA_REAL(nid,varid,corner,edges,zx)

      if (ierr.ne.NF_NOERR) then
           write(*,*) "***** PUT_VAR problem in writediagfi"
           write(*,*) "***** with ",nom
           write(*,*) 'ierr=', ierr
      endif 


      end
