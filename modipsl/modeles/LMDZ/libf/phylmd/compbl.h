      !
      ! $Header$
      !
!jyg+al1<
!!      integer iflag_pbl,iflag_pbl_split
!!      common/compbl/iflag_pbl,iflag_pbl_split
!!FC      integer iflag_pbl, iflag_pbl_split, iflag_order2_sollw
!FC      common/compbl/iflag_pbl, iflag_pbl_split, iflag_order2_sollw
      integer iflag_pbl, iflag_pbl_split, iflag_order2_sollw, ifl_pbltree
      common/compbl/iflag_pbl, iflag_pbl_split, iflag_order2_sollw, ifl_pbltree
!>jyg+al1
!$OMP THREADPRIVATE(/compbl/)
