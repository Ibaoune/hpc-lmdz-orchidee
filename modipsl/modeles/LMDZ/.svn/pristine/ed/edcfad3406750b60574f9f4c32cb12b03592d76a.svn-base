MODULE YOMCOAPHY

USE PARKIND1  ,ONLY : JPIM
USE GRIDPOINT_BUFFERS , ONLY : gridpoint_buffer

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOEPHY* - SWITCHES RELATED TO DIABATIC PROCESSES
!     -----------------------------------------------------------------

!        * E.C.M.W.F. PHYSICS PACKAGE *

INTEGER(KIND=JPIM) :: NPHYINT      ! NPHYINT=0 -> Physics and dynamics at the same resolution
                                   ! NPHYINT=1 -> Physics grid coarser than dynamics grid
				   ! NPHYINT=2 -> Physics grid finer than dynamics grid
TYPE(gridpoint_buffer) :: PHYS_GPPBUF
CHARACTER (LEN = 32) ::   CPTABLEFIL
CHARACTER (LEN = 256) ::   CPTABLEDIR

!$OMP THREADPRIVATE(cptabledir,cptablefil,nphyint,phys_gppbuf)
END MODULE YOMCOAPHY
