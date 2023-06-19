MODULE ITEMPORAL_SPLITTING
   USE, INTRINSIC :: ISO_C_BINDING
   USE TEMPORAL_SPLITTING_INTERFACE
   TYPE xios_temporal_splitting
      INTEGER(kind = C_INTPTR_T) :: daddr
   END TYPE xios_temporal_splitting
   CONTAINS ! Fonctions disponibles pour les utilisateurs.
   SUBROUTINE xios_get_temporal_splitting_handle(idt,ret)
      IMPLICIT NONE
      CHARACTER(len = *), INTENT(IN) :: idt
      TYPE(xios_temporal_splitting) , INTENT(OUT):: ret
      CALL cxios_temporal_splitting_handle_create(ret%daddr, idt, len(idt))
   END SUBROUTINE xios_get_temporal_splitting_handle
   LOGICAL FUNCTION xios_is_valid_temporal_splitting(idt)
      IMPLICIT NONE
      CHARACTER(len = *) , INTENT(IN) :: idt
      LOGICAL (kind = 1) :: val
      CALL cxios_temporal_splitting_valid_id(val, idt, len(idt))
      xios_is_valid_temporal_splitting = val
   END FUNCTION xios_is_valid_temporal_splitting
END MODULE ITEMPORAL_SPLITTING
