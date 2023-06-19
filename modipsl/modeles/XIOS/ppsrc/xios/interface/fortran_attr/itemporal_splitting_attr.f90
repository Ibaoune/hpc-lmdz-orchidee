! * ************************************************************************** *
! * Interface auto generated - do not modify *
! * ************************************************************************** *
MODULE itemporal_splitting_attr
  USE, INTRINSIC :: ISO_C_BINDING
  USE itemporal_splitting
  USE temporal_splitting_interface_attr
CONTAINS
  SUBROUTINE xios_set_temporal_splitting_attr &
    ( temporal_splitting_id )
    IMPLICIT NONE
      TYPE(xios_temporal_splitting) :: temporal_splitting_hdl
      CHARACTER(LEN=*), INTENT(IN) ::temporal_splitting_id
      CALL xios_get_temporal_splitting_handle &
      (temporal_splitting_id,temporal_splitting_hdl)
      CALL xios_set_temporal_splitting_attr_hdl_ &
      ( temporal_splitting_hdl )
  END SUBROUTINE xios_set_temporal_splitting_attr
  SUBROUTINE xios_set_temporal_splitting_attr_hdl &
    ( temporal_splitting_hdl )
    IMPLICIT NONE
      TYPE(xios_temporal_splitting) , INTENT(IN) :: temporal_splitting_hdl
      CALL xios_set_temporal_splitting_attr_hdl_ &
      ( temporal_splitting_hdl )
  END SUBROUTINE xios_set_temporal_splitting_attr_hdl
  SUBROUTINE xios_set_temporal_splitting_attr_hdl_ &
    ( temporal_splitting_hdl )
    IMPLICIT NONE
      TYPE(xios_temporal_splitting) , INTENT(IN) :: temporal_splitting_hdl
  END SUBROUTINE xios_set_temporal_splitting_attr_hdl_
  SUBROUTINE xios_get_temporal_splitting_attr &
    ( temporal_splitting_id )
    IMPLICIT NONE
      TYPE(xios_temporal_splitting) :: temporal_splitting_hdl
      CHARACTER(LEN=*), INTENT(IN) ::temporal_splitting_id
      CALL xios_get_temporal_splitting_handle &
      (temporal_splitting_id,temporal_splitting_hdl)
      CALL xios_get_temporal_splitting_attr_hdl_ &
      ( temporal_splitting_hdl )
  END SUBROUTINE xios_get_temporal_splitting_attr
  SUBROUTINE xios_get_temporal_splitting_attr_hdl &
    ( temporal_splitting_hdl )
    IMPLICIT NONE
      TYPE(xios_temporal_splitting) , INTENT(IN) :: temporal_splitting_hdl
      CALL xios_get_temporal_splitting_attr_hdl_ &
      ( temporal_splitting_hdl )
  END SUBROUTINE xios_get_temporal_splitting_attr_hdl
  SUBROUTINE xios_get_temporal_splitting_attr_hdl_ &
    ( temporal_splitting_hdl )
    IMPLICIT NONE
      TYPE(xios_temporal_splitting) , INTENT(IN) :: temporal_splitting_hdl
  END SUBROUTINE xios_get_temporal_splitting_attr_hdl_
  SUBROUTINE xios_is_defined_temporal_splitting_attr &
    ( temporal_splitting_id )
    IMPLICIT NONE
      TYPE(xios_temporal_splitting) :: temporal_splitting_hdl
      CHARACTER(LEN=*), INTENT(IN) ::temporal_splitting_id
      CALL xios_get_temporal_splitting_handle &
      (temporal_splitting_id,temporal_splitting_hdl)
      CALL xios_is_defined_temporal_splitting_attr_hdl_ &
      ( temporal_splitting_hdl )
  END SUBROUTINE xios_is_defined_temporal_splitting_attr
  SUBROUTINE xios_is_defined_temporal_splitting_attr_hdl &
    ( temporal_splitting_hdl )
    IMPLICIT NONE
      TYPE(xios_temporal_splitting) , INTENT(IN) :: temporal_splitting_hdl
      CALL xios_is_defined_temporal_splitting_attr_hdl_ &
      ( temporal_splitting_hdl )
  END SUBROUTINE xios_is_defined_temporal_splitting_attr_hdl
  SUBROUTINE xios_is_defined_temporal_splitting_attr_hdl_ &
    ( temporal_splitting_hdl )
    IMPLICIT NONE
      TYPE(xios_temporal_splitting) , INTENT(IN) :: temporal_splitting_hdl
  END SUBROUTINE xios_is_defined_temporal_splitting_attr_hdl_
END MODULE itemporal_splitting_attr
