# Test Couplé SSP2-4.5 3 Jours (`ssp245_3days`)

## 1. Objectif
Tester et valider sur le cluster HPC Toubkal (UM6P) la transition entre la simulation historique (terminée au 2014-12-31) et le scénario climatique futur moyen SSP2-4.5 sur 3 jours consécutifs (2015-01-01 au 2015-01-03).

## 2. Configuration Scientifique et Technique
- **Période :** 2015-01-01 au 2015-01-03 (`PeriodLength = 1D`, `PeriodNb = 3`).
- **Restarts utilisés (date 2014-12-31, fin de HIST, via `simulations/lmdzor_restart_stage/sim35km-HIST-after2002`) :**
  - ATM : `sim35km-HIST-after2002_20141231_restart.nc` (433 Mo) et `restartphy.nc` (1,66 Go)
  - SRF : `sim35km-HIST-after2002_20141231_sechiba_rest.nc` (591 Mo)
  - SBG : `sim35km-HIST-after2002_20141231_stomate_rest.nc` (1,52 Go)
- **Forçages atmosphériques SSP2-4.5 :**
  - GES : `CO2_CMIP6_2015_2500.txt`, `CH4`, `N2O`, `CFC` sous `ATM/GHG/CMIP6/ssp/ssp245/`
  - Aérosols : `aerosols2015_from_inca.nc` (4,8 Go) sous `ATM/AEROSOLS/CMIP6/Scenarios/ssp245/384x220/L79/`
  - Conditions aux limites : `ELI-384x220x79-ssp245-T01_2015_limit.nc` (1,9 Go) et `climoz_LMDZ.nc` (313 Mo)
  - Solaire : `solarforcing_daily_3.2_6bands_2015.nc`
- **Ressources Slurm :**
  - 28 nœuds exclusifs (partition `compute`).
  - 111 rangs MPI : 110 rangs LMDZ/ORCHIDEE (14 threads OpenMP par rang, 4 rangs par nœud = 56 cœurs) + 1 rang serveur XIOS dédié.
  - Walltime : `01:15:00` (75 minutes).
  - Lanceur Slurm natif : `srun --mpi=pmix` avec distribution `multi-prog` et affinité `taskset`.

## 3. Forçages de Végétation et Distinction avec la Simulation Historique
- **Forçages nominaux utilisés dans ce test (`COMP/orchidee.card`) :**
  - `PFTmap` : `${R_IN}/SRF/PFTMAPS/CMIP6/ESA-LUH2v2/scenarios/ssp245/PFTmap_2015.nc` (62 Mo).
  - `woodharvest` : `${R_IN}/SRF/WOODHARVEST/LUH2v2/scenarios/ssp245/woodharvest_2015.nc` (8,3 Mo).
- **Différence avec la simulation historique archivée :**
  - La simulation archivée de production (`experiments/ssp245/`, `sim35km-HIST-future-ssp245`) pointait par erreur vers `scenarios/ssp585/PFTmap_${year}.nc` et `scenarios/ssp585/woodharvest_${year}.nc`.
  - Pour préserver la traçabilité historique stricte, `experiments/ssp245/COMP/orchidee.card` est conservé inchangé avec ses chemins d'origine `ssp585`.
  - En revanche, le présent test nominal `tests/ssp245_3days/` utilise les véritables forçages nominaux du scénario SSP245.

## 4. Résultats de Validation (Job `8235853`)
La simulation des 3 jours a été exécutée et validée avec succès le 2026-11-06 :
- **Job ID :** `8235853`
- **Statut Slurm :** `COMPLETED` (Code de sortie `0:0`)
- **Durée totale du job :** 36 minutes 03 secondes (début 17:34:14, fin 18:10:17)
- **Détail des 3 périodes au sein du même job :**
  - **Jour 1 (2015-01-01) :** Étape Slurm `8235853.0`, calcul modèle 14 min 43 s (`RealCpuTime = 882.70 s`), sortie `COMPLETED 0:0`.
  - **Jour 2 (2015-01-02) :** Étape Slurm `8235853.1`, calcul modèle 6 min 19 s (`RealCpuTime = 378.66 s`), sortie `COMPLETED 0:0`.
  - **Jour 3 (2015-01-03) :** Étape Slurm `8235853.2`, calcul modèle 6 min 17 s (`RealCpuTime = 377.03 s`), sortie `COMPLETED 0:0`.
- **État final `run.card` :** `CumulPeriod = 4`, `PeriodDateBegin = 2015-01-04`, `PeriodState = Completed`.

## 5. Sorties et Restarts Archivés (`IGCM_OUT/`)
Toutes les sorties scientifiques et les fichiers de restart ont été archivés sous :
`$IGCM_OUT/LMDZOR/TEST/amip/test-3days-SSP245-r35/`

- **Restarts journaliers produits (validés via `ncdump -h`) :**
  - `ATM/Restart/test-3days-SSP245-r35_2015010[1-3]_restart.nc` (414 Mo)
  - `ATM/Restart/test-3days-SSP245-r35_2015010[1-3]_restartphy.nc` (1,6 Go)
  - `SRF/Restart/test-3days-SSP245-r35_2015010[1-3]_sechiba_rest.nc` (565 Mo)
  - `SBG/Restart/test-3days-SSP245-r35_2015010[1-3]_stomate_rest.nc` (1,5 Go)
- **Fichiers historiques NetCDF :**
  - `ATM/Output/DA/test-3days-SSP245-r35_2015010[1-3]_2015010[1-3]_1D_histday.nc` (1,2 Go/jour)
  - `ATM/Output/MO/test-3days-SSP245-r35_2015010[1-3]_2015010[1-3]_1M_histmth.nc` (1,7 Mo/jour)
  - Sorties ORCHIDEE / STOMATE haute fréquence et mensuelles (`1M_stomate_history.nc`, `HF_sechiba_out_2.nc`, `1M_sechiba_history.nc`).
## 6. Distinction Explicite entre Validation Technique et Validation Scientifique
- **Validation Technique (Réussie) :**
  - Transition opérationnelle validée entre HIST (2014-12-31) et le scénario SSP2-4.5 sur 3 jours.
  - Utilisation des forçages nominaux SSP245 (`scenarios/ssp245` pour `PFTmap` et `woodharvest`).
  - Enchaînement séquentiel parfait des 3 périodes sous Slurm (`srun --mpi=pmix`) sans reliquat de processus orphelins.
  - Intégrité NetCDF des fichiers archivés sous `IGCM_OUT/`.
- **Validation Scientifique (Non couverte par un test de 3 jours) :**
  - Ce banc d'essai de 3 jours ne permet pas d'évaluer la trajectoire radiative, le réalisme de la réponse climatique ni l'impact de la végétation sous scénario SSP2-4.5.
  - Une validation scientifique nécessitera des simulations sur plusieurs décennies avec analyse des bilans d'énergie et de masse.
