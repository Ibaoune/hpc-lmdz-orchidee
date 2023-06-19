SUBROUTINE print_debug_phys (i,debug_lev,text)

USE dimphy, ONLY: klev
USE phys_local_var_mod, ONLY: u_seri, v_seri, t_seri, q_seri, ql_seri
USE geometry_mod, ONLY: longitude_deg, latitude_deg
IMPLICIT NONE
integer i,debug_lev
CHARACTER*(*) text


integer k

print*,'PLANTAGE POUR LE POINT i=',i,text
print*,'l    u, v, T, q, ql'
DO k = 1, klev
   write(*,'(i3,2f8.4,3f14.4,2e14.2)') k,longitude_deg(i),latitude_deg(i), &
   u_seri(i,k),v_seri(i,k),t_seri(i,k),q_seri(i,k),ql_seri(i,k)
ENDDO

RETURN
END
