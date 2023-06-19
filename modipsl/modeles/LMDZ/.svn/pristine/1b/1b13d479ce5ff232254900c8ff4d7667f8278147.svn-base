      MODULE indice_sol_mod

            IMPLICIT NONE

            INTEGER, PARAMETER :: nbsrf = 4 ! nombre de sous-fractions pour une maille
            INTEGER, PARAMETER :: is_ter=1, is_lic=2, is_oce = 3, is_sic=4
                              ! terre ! ocean ! glacier continental ! glace de mer

            INTEGER, PARAMETER :: is_ave=nbsrf+1 ! valeur moyenne sur l'ensemble des surfaces
            REAL, PARAMETER :: epsfra=1.0E-05

            CHARACTER(len=3), DIMENSION(nbsrf), PARAMETER :: clnsurf = (/'ter', 'lic', 'oce', 'sic'/)
!FC 
           INTEGER, SAVE    :: nvm_orch ! Nombre de type de vegetation ds ORCHIDEE                 

      END MODULE indice_sol_mod

