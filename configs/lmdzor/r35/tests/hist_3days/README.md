# Test Couplé HIST 3 Jours (`hist_3days`)

## 1. Objectif
Tester sur Toubkal la chaîne couplée complète LMDZ–ORCHIDEE–STOMATE–XIOS sur le domaine zoomé 35 km ($384 \times 220 \times 79$) sur une durée de 3 jours consécutifs.

## 2. Configuration (`config.card`)
- **Période :** 1980-01-01 au 1980-01-03 (`PeriodLength = 1D`, 3 pas journaliers).
- **Restarts utilisés (date 1979-12-31, archivés sur Lustre) :**
  - ATM : `/srv/lustre01/.../ARCHIVE/Simulations/LMDZ/r35/hist/ATM/Restart/sim35km-HIST_19791231_restart.nc` et `restartphy.nc`.
  - SRF : `/srv/lustre01/.../ARCHIVE/Simulations/LMDZ/r35/hist/SRF/Restart/sim35km-HIST_19791231_sechiba_rest.nc`.
  - SBG : `/srv/lustre01/.../ARCHIVE/Simulations/LMDZ/r35/hist/SBG/Restart/sim35km-HIST_19791231_stomate_rest.nc`.
- **Conditions aux limites (année 1980) :** `limit.nc`, `climoz_LMDZ.nc`, aérosols INCA 1980, PFTmap 1980.
- **Ressources :** 14 nœuds, 111 MPI, partition `compute`.

## 3. Lancement
Le script Slurm `Job_test-3days-HIST-r35` a été généré via `ins_job`.
