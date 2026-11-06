# Test Couplé SSP5-8.5 3 Jours (`ssp585_3days`)

## 1. Objectif
Tester et valider sur le cluster HPC Toubkal (UM6P) la transition entre la simulation historique (terminée au 2014-12-31) et le scénario climatique futur de fort réchauffement SSP5-8.5 sur 3 jours consécutifs (2015-01-01 au 2015-01-03).

## 2. Configuration Scientifique et Technique
- **Période :** 2015-01-01 au 2015-01-03 (`PeriodLength = 1D`, `PeriodNb = 3`).
- **Restarts utilisés (date 2014-12-31, fin de HIST, via `simulations/lmdzor_restart_stage/sim35km-HIST-after2002`) :**
  - ATM : `sim35km-HIST-after2002_20141231_restart.nc` (433 Mo) et `restartphy.nc` (1,66 Go)
  - SRF : `sim35km-HIST-after2002_20141231_sechiba_rest.nc` (591 Mo)
  - SBG : `sim35km-HIST-after2002_20141231_stomate_rest.nc` (1,52 Go)
- **Forçages atmosphériques et surface SSP5-8.5 :**
  - GES : `CO2_CMIP6_2015_2500.txt`, `CH4`, `N2O`, `CFC` sous `ATM/GHG/CMIP6/ssp/ssp585/`
  - Aérosols : `aerosols2015_from_inca.nc` (4,6 Go) sous `ATM/AEROSOLS/CMIP6/Scenarios/ssp585/384x220/L79/`
  - Conditions aux limites : `ELI-384x220x79-ssp585-T01_2015_limit.nc` (1,9 Go) et `ELI-384x220x79-ssp585-T01_2015_climoz_LMDZ.nc` (300 Mo)
  - Solaire : `solarforcing_daily_3.2_6bands_2015.nc`
  - Végétation : `PFTmap_2015.nc` (60 Mo) et `woodharvest_2015.nc` (8,0 Mo) sous `SRF/PFTMAPS/CMIP6/ESA-LUH2v2/scenarios/ssp585/`
- **Ressources Slurm :**
  - 28 nœuds exclusifs (partition `compute`).
  - 111 rangs MPI : 110 rangs LMDZ/ORCHIDEE (14 threads OpenMP par rang, 4 rangs par nœud = 56 cœurs) + 1 rang serveur XIOS dédié.
  - Walltime : `01:15:00` (75 minutes).
  - Lanceur Slurm natif : `srun --mpi=pmix` avec distribution `multi-prog` et affinité `taskset`.
  - Environnement partagé : `common/env_toubkal.sh` (GCC 7.3.0, OpenMPI 3.1.1, NetCDF-C 4.6.1, NetCDF-Fortran 4.4.4, FFTW 3.3.8 - suite foss-2018b).

## 3. Résultats de Validation (Job `8235912`)
La simulation des 3 jours a été exécutée et validée avec succès le 2026-11-06 :
- **Job ID :** `8235912`
- **Statut Slurm :** `COMPLETED` (Code de sortie `0:0`)
- **Durée totale du job :** 36 minutes 03 secondes (début 18:30:20, fin 19:06:23)
- **Détail des 3 périodes au sein du même job :**
  - **Jour 1 (2015-01-01) :** Étape Slurm `8235912.0`, calcul modèle 15 min 09 s (`RealCpuTime = 908.45 s`), sortie `COMPLETED 0:0`.
  - **Jour 2 (2015-01-02) :** Étape Slurm `8235912.1`, calcul modèle 6 min 10 s (`RealCpuTime = 369.99 s`), sortie `COMPLETED 0:0`.
  - **Jour 3 (2015-01-03) :** Étape Slurm `8235912.2`, calcul modèle 6 min 12 s (`RealCpuTime = 371.97 s`), sortie `COMPLETED 0:0`.
- **État final `run.card` :** `CumulPeriod = 4`, `PeriodDateBegin = 2015-01-04`, `PeriodState = Completed`.

## 4. Sorties et Restarts Archivés (`IGCM_OUT/`)
Toutes les sorties scientifiques et les fichiers de restart ont été archivés sous :
`$IGCM_OUT/LMDZOR/TEST/amip/test-3days-SSP585-r35/`

- **Restarts journaliers produits (validés via `ncdump -h`) :**
  - `ATM/Restart/test-3days-SSP585-r35_2015010[1-3]_restart.nc` (414 Mo)
  - `ATM/Restart/test-3days-SSP585-r35_2015010[1-3]_restartphy.nc` (1,6 Go)
  - `SRF/Restart/test-3days-SSP585-r35_2015010[1-3]_sechiba_rest.nc` (565 Mo)
  - `SBG/Restart/test-3days-SSP585-r35_2015010[1-3]_stomate_rest.nc` (1,5 Go)
- **Fichiers historiques NetCDF :**
  - `ATM/Output/DA/test-3days-SSP585-r35_2015010[1-3]_2015010[1-3]_1D_histday.nc` (1,2 Go/jour)
  - `ATM/Output/MO/test-3days-SSP585-r35_2015010[1-3]_2015010[1-3]_1M_histmth.nc` (1,7 Mo/jour)
  - Sorties ORCHIDEE / STOMATE haute fréquence et mensuelles (`1M_stomate_history.nc`, `HF_sechiba_out_2.nc`, `1M_sechiba_history.nc`, `1M_sechiba_history_4dim.nc`).

## 5. Distinction Explicite entre Validation Technique et Validation Scientifique
- **Validation Technique (Réussie) :**
  - Intégration et chaîne de calcul validées sur l'infrastructure HPC Toubkal.
  - Résolution complète des forçages extrêmes SSP5-8.5 et des restarts HIST 2014-12-31.
  - Enchaînement séquentiel sans blocage de 3 périodes consécutives sous Slurm (`srun --mpi=pmix`).
  - Intégrité NetCDF des fichiers archivés (restarts ATM, SRF, SBG et sorties diagnostiques).
- **Validation Scientifique (Non couverte par un banc d'essai de 3 jours) :**
  - Ce test de 3 jours confirme uniquement la santé technique et la stabilité dynamique du code.
  - Il ne permet pas d'évaluer la qualité du climat simulé sous scénario SSP5-8.5, l'ajustement du cycle saisonnier, l'équilibre des flux de surface ni les réponses radiatives à long terme.
  - Toute analyse scientifique exigera une intégration pluriannuelle ou décennale confrontée aux données d'observation et aux protocoles CMIP6.
