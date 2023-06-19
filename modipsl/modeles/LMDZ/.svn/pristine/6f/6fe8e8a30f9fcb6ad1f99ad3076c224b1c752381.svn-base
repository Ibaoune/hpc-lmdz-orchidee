  subroutine global_mean(field,airephy,laire,mfield)
!
! I.Musat: 05.2011
! calcul moyenne globale d'un champ pondere par l'aire de la maille
! (laire=.TRUE.) ou somme globale du champ (laire=.FALSE.)
!
  USE dimphy
  USE mod_phys_lmdz_para, only: is_sequential
  USE mod_phys_lmdz_transfert_para, only: reduce_sum
  use mod_phys_lmdz_mpi_data, only: is_mpi_root
  USE ioipsl
  implicit none

  real,dimension(klon),intent(in) :: field
  real,dimension(klon),intent(in) :: airephy
  LOGICAL, intent(in) :: laire
  REAL, intent(out) :: mfield
  REAL :: airetot     ! Total area the earth
  REAL :: sumtmp
  INTEGER :: i

  if (is_sequential) then

   airetot = 0.
   sumtmp = 0.
   DO i=1, klon
    airetot = airetot + airephy(i)
    sumtmp = sumtmp + field(i)
   END DO
   if (laire) THEN
    if(airetot.NE.0.) THEN
     mfield=sumtmp/airetot
    endif
   else
    mfield=sumtmp
   endif

  else

   CALL reduce_sum(SUM(airephy),airetot)
   CALL reduce_sum(SUM(field),sumtmp)

!$OMP MASTER
  if (is_mpi_root) THEN
  if (laire) THEN
!  print*,'gmean airetot=',airetot
   if(airetot.NE.0.) THEN
    mfield=sumtmp/airetot
!  else
!   mfield=sumtmp
   endif
  else
   mfield=sumtmp
  endif

! print*,'gmean sumtmp mfield=',sumtmp,mfield

  endif !(is_mpi_root) THEN
!$OMP END MASTER

  endif

  end subroutine global_mean
