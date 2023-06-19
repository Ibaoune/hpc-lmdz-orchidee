! * ************************************************************************** *
! * Interface auto generated - do not modify *
! * ************************************************************************** *
MODULE ireduce_axis_to_axis_attr
  USE, INTRINSIC :: ISO_C_BINDING
  USE ireduce_axis_to_axis
  USE reduce_axis_to_axis_interface_attr
CONTAINS
  SUBROUTINE xios_set_reduce_axis_to_axis_attr &
    ( reduce_axis_to_axis_id, operation )
    IMPLICIT NONE
      TYPE(xios_reduce_axis_to_axis) :: reduce_axis_to_axis_hdl
      CHARACTER(LEN=*), INTENT(IN) ::reduce_axis_to_axis_id
      CHARACTER(len = *) , OPTIONAL, INTENT(IN) :: operation
      CALL xios_get_reduce_axis_to_axis_handle &
      (reduce_axis_to_axis_id,reduce_axis_to_axis_hdl)
      CALL xios_set_reduce_axis_to_axis_attr_hdl_ &
      ( reduce_axis_to_axis_hdl, operation )
  END SUBROUTINE xios_set_reduce_axis_to_axis_attr
  SUBROUTINE xios_set_reduce_axis_to_axis_attr_hdl &
    ( reduce_axis_to_axis_hdl, operation )
    IMPLICIT NONE
      TYPE(xios_reduce_axis_to_axis) , INTENT(IN) :: reduce_axis_to_axis_hdl
      CHARACTER(len = *) , OPTIONAL, INTENT(IN) :: operation
      CALL xios_set_reduce_axis_to_axis_attr_hdl_ &
      ( reduce_axis_to_axis_hdl, operation )
  END SUBROUTINE xios_set_reduce_axis_to_axis_attr_hdl
  SUBROUTINE xios_set_reduce_axis_to_axis_attr_hdl_ &
    ( reduce_axis_to_axis_hdl, operation_ )
    IMPLICIT NONE
      TYPE(xios_reduce_axis_to_axis) , INTENT(IN) :: reduce_axis_to_axis_hdl
      CHARACTER(len = *) , OPTIONAL, INTENT(IN) :: operation_
      IF (PRESENT(operation_)) THEN
        CALL cxios_set_reduce_axis_to_axis_operation &
      (reduce_axis_to_axis_hdl%daddr, operation_, len(operation_))
      ENDIF
  END SUBROUTINE xios_set_reduce_axis_to_axis_attr_hdl_
  SUBROUTINE xios_get_reduce_axis_to_axis_attr &
    ( reduce_axis_to_axis_id, operation )
    IMPLICIT NONE
      TYPE(xios_reduce_axis_to_axis) :: reduce_axis_to_axis_hdl
      CHARACTER(LEN=*), INTENT(IN) ::reduce_axis_to_axis_id
      CHARACTER(len = *) , OPTIONAL, INTENT(OUT) :: operation
      CALL xios_get_reduce_axis_to_axis_handle &
      (reduce_axis_to_axis_id,reduce_axis_to_axis_hdl)
      CALL xios_get_reduce_axis_to_axis_attr_hdl_ &
      ( reduce_axis_to_axis_hdl, operation )
  END SUBROUTINE xios_get_reduce_axis_to_axis_attr
  SUBROUTINE xios_get_reduce_axis_to_axis_attr_hdl &
    ( reduce_axis_to_axis_hdl, operation )
    IMPLICIT NONE
      TYPE(xios_reduce_axis_to_axis) , INTENT(IN) :: reduce_axis_to_axis_hdl
      CHARACTER(len = *) , OPTIONAL, INTENT(OUT) :: operation
      CALL xios_get_reduce_axis_to_axis_attr_hdl_ &
      ( reduce_axis_to_axis_hdl, operation )
  END SUBROUTINE xios_get_reduce_axis_to_axis_attr_hdl
  SUBROUTINE xios_get_reduce_axis_to_axis_attr_hdl_ &
    ( reduce_axis_to_axis_hdl, operation_ )
    IMPLICIT NONE
      TYPE(xios_reduce_axis_to_axis) , INTENT(IN) :: reduce_axis_to_axis_hdl
      CHARACTER(len = *) , OPTIONAL, INTENT(OUT) :: operation_
      IF (PRESENT(operation_)) THEN
        CALL cxios_get_reduce_axis_to_axis_operation &
      (reduce_axis_to_axis_hdl%daddr, operation_, len(operation_))
      ENDIF
  END SUBROUTINE xios_get_reduce_axis_to_axis_attr_hdl_
  SUBROUTINE xios_is_defined_reduce_axis_to_axis_attr &
    ( reduce_axis_to_axis_id, operation )
    IMPLICIT NONE
      TYPE(xios_reduce_axis_to_axis) :: reduce_axis_to_axis_hdl
      CHARACTER(LEN=*), INTENT(IN) ::reduce_axis_to_axis_id
      LOGICAL, OPTIONAL, INTENT(OUT) :: operation
      LOGICAL(KIND=C_BOOL) :: operation_tmp
      CALL xios_get_reduce_axis_to_axis_handle &
      (reduce_axis_to_axis_id,reduce_axis_to_axis_hdl)
      CALL xios_is_defined_reduce_axis_to_axis_attr_hdl_ &
      ( reduce_axis_to_axis_hdl, operation )
  END SUBROUTINE xios_is_defined_reduce_axis_to_axis_attr
  SUBROUTINE xios_is_defined_reduce_axis_to_axis_attr_hdl &
    ( reduce_axis_to_axis_hdl, operation )
    IMPLICIT NONE
      TYPE(xios_reduce_axis_to_axis) , INTENT(IN) :: reduce_axis_to_axis_hdl
      LOGICAL, OPTIONAL, INTENT(OUT) :: operation
      LOGICAL(KIND=C_BOOL) :: operation_tmp
      CALL xios_is_defined_reduce_axis_to_axis_attr_hdl_ &
      ( reduce_axis_to_axis_hdl, operation )
  END SUBROUTINE xios_is_defined_reduce_axis_to_axis_attr_hdl
  SUBROUTINE xios_is_defined_reduce_axis_to_axis_attr_hdl_ &
    ( reduce_axis_to_axis_hdl, operation_ )
    IMPLICIT NONE
      TYPE(xios_reduce_axis_to_axis) , INTENT(IN) :: reduce_axis_to_axis_hdl
      LOGICAL, OPTIONAL, INTENT(OUT) :: operation_
      LOGICAL(KIND=C_BOOL) :: operation__tmp
      IF (PRESENT(operation_)) THEN
        operation__tmp = cxios_is_defined_reduce_axis_to_axis_operation &
      (reduce_axis_to_axis_hdl%daddr)
        operation_ = operation__tmp
      ENDIF
  END SUBROUTINE xios_is_defined_reduce_axis_to_axis_attr_hdl_
END MODULE ireduce_axis_to_axis_attr
