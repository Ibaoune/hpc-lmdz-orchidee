! * ************************************************************************** *
! * Interface auto generated - do not modify *
! * ************************************************************************** *
MODULE iduplicate_scalar_to_axis_attr
  USE, INTRINSIC :: ISO_C_BINDING
  USE iduplicate_scalar_to_axis
  USE duplicate_scalar_to_axis_interface_attr
CONTAINS
  SUBROUTINE xios_set_duplicate_scalar_to_axis_attr &
    ( duplicate_scalar_to_axis_id )
    IMPLICIT NONE
      TYPE(xios_duplicate_scalar_to_axis) :: duplicate_scalar_to_axis_hdl
      CHARACTER(LEN=*), INTENT(IN) ::duplicate_scalar_to_axis_id
      CALL xios_get_duplicate_scalar_to_axis_handle &
      (duplicate_scalar_to_axis_id,duplicate_scalar_to_axis_hdl)
      CALL xios_set_duplicate_scalar_to_axis_attr_hdl_ &
      ( duplicate_scalar_to_axis_hdl )
  END SUBROUTINE xios_set_duplicate_scalar_to_axis_attr
  SUBROUTINE xios_set_duplicate_scalar_to_axis_attr_hdl &
    ( duplicate_scalar_to_axis_hdl )
    IMPLICIT NONE
      TYPE(xios_duplicate_scalar_to_axis) , INTENT(IN) :: duplicate_scalar_to_axis_hdl
      CALL xios_set_duplicate_scalar_to_axis_attr_hdl_ &
      ( duplicate_scalar_to_axis_hdl )
  END SUBROUTINE xios_set_duplicate_scalar_to_axis_attr_hdl
  SUBROUTINE xios_set_duplicate_scalar_to_axis_attr_hdl_ &
    ( duplicate_scalar_to_axis_hdl )
    IMPLICIT NONE
      TYPE(xios_duplicate_scalar_to_axis) , INTENT(IN) :: duplicate_scalar_to_axis_hdl
  END SUBROUTINE xios_set_duplicate_scalar_to_axis_attr_hdl_
  SUBROUTINE xios_get_duplicate_scalar_to_axis_attr &
    ( duplicate_scalar_to_axis_id )
    IMPLICIT NONE
      TYPE(xios_duplicate_scalar_to_axis) :: duplicate_scalar_to_axis_hdl
      CHARACTER(LEN=*), INTENT(IN) ::duplicate_scalar_to_axis_id
      CALL xios_get_duplicate_scalar_to_axis_handle &
      (duplicate_scalar_to_axis_id,duplicate_scalar_to_axis_hdl)
      CALL xios_get_duplicate_scalar_to_axis_attr_hdl_ &
      ( duplicate_scalar_to_axis_hdl )
  END SUBROUTINE xios_get_duplicate_scalar_to_axis_attr
  SUBROUTINE xios_get_duplicate_scalar_to_axis_attr_hdl &
    ( duplicate_scalar_to_axis_hdl )
    IMPLICIT NONE
      TYPE(xios_duplicate_scalar_to_axis) , INTENT(IN) :: duplicate_scalar_to_axis_hdl
      CALL xios_get_duplicate_scalar_to_axis_attr_hdl_ &
      ( duplicate_scalar_to_axis_hdl )
  END SUBROUTINE xios_get_duplicate_scalar_to_axis_attr_hdl
  SUBROUTINE xios_get_duplicate_scalar_to_axis_attr_hdl_ &
    ( duplicate_scalar_to_axis_hdl )
    IMPLICIT NONE
      TYPE(xios_duplicate_scalar_to_axis) , INTENT(IN) :: duplicate_scalar_to_axis_hdl
  END SUBROUTINE xios_get_duplicate_scalar_to_axis_attr_hdl_
  SUBROUTINE xios_is_defined_duplicate_scalar_to_axis_attr &
    ( duplicate_scalar_to_axis_id )
    IMPLICIT NONE
      TYPE(xios_duplicate_scalar_to_axis) :: duplicate_scalar_to_axis_hdl
      CHARACTER(LEN=*), INTENT(IN) ::duplicate_scalar_to_axis_id
      CALL xios_get_duplicate_scalar_to_axis_handle &
      (duplicate_scalar_to_axis_id,duplicate_scalar_to_axis_hdl)
      CALL xios_is_defined_duplicate_scalar_to_axis_attr_hdl_ &
      ( duplicate_scalar_to_axis_hdl )
  END SUBROUTINE xios_is_defined_duplicate_scalar_to_axis_attr
  SUBROUTINE xios_is_defined_duplicate_scalar_to_axis_attr_hdl &
    ( duplicate_scalar_to_axis_hdl )
    IMPLICIT NONE
      TYPE(xios_duplicate_scalar_to_axis) , INTENT(IN) :: duplicate_scalar_to_axis_hdl
      CALL xios_is_defined_duplicate_scalar_to_axis_attr_hdl_ &
      ( duplicate_scalar_to_axis_hdl )
  END SUBROUTINE xios_is_defined_duplicate_scalar_to_axis_attr_hdl
  SUBROUTINE xios_is_defined_duplicate_scalar_to_axis_attr_hdl_ &
    ( duplicate_scalar_to_axis_hdl )
    IMPLICIT NONE
      TYPE(xios_duplicate_scalar_to_axis) , INTENT(IN) :: duplicate_scalar_to_axis_hdl
  END SUBROUTINE xios_is_defined_duplicate_scalar_to_axis_attr_hdl_
END MODULE iduplicate_scalar_to_axis_attr
