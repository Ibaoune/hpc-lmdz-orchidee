# Test Couplé HIST 3 Jours (`hist_3days`)

## 1. Objectif
Valider sur le cluster HPC Toubkal (UM6P) l'exécution complète et la finalisation de la chaîne couplée LMDZ–ORCHIDEE–STOMATE–XIOS sur le domaine zoomé Maroc 35 km ($384 \times 220 \times 79$) sur trois journées consécutives.

## 2. Configuration Scientifique et Technique
- **Période :** Du 1980-01-01 au 1980-01-03 (`PeriodLength = 1D`, `PeriodNb = 3`).
- **Restarts initiaux (1979-12-31, via `simulations/lmdzor_restart_stage/sim35km-HIST`) :**
  - ATM : `sim35km-HIST_19791231_restart.nc` et `restartphy.nc`
  - SRF : `sim35km-HIST_19791231_sechiba_rest.nc`
  - SBG : `sim35km-HIST_19791231_stomate_rest.nc`
- **Forçages et conditions aux limites (1980) :** `limit.nc` (SST/glace AMIP), `climoz_LMDZ.nc` (ozone), `aerosols1980.nc` / `aerosols.nat.nc`, `PFTmap.nc`, `woodharvest.nc`, `solarforcing.nc`.
- **Ressources Slurm :**
  - 28 nœuds exclusifs (partition `compute`).
  - 111 rangs MPI : 110 rangs LMDZ/ORCHIDEE (14 threads OpenMP par rang, 4 rangs par nœud = 56 cœurs) + 1 rang serveur XIOS dédié.
  - Walltime : `01:15:00` (75 minutes).

## 3. Architecture d'Exécution et Finalisation MPI
Suite au diagnostic du Job `8235452` (blocage au shutdown d'OpenMPI 3.1.1 `mpirun`/`orted`), la chaîne de lancement dans `libIGCM_sys_toubkal.ksh` a été corrigée pour utiliser le lanceur Slurm natif :
- **Commande de lancement :**
  ```ksh
  typeset -r HOST_MPIRUN_COMMAND="/usr/bin/time srun --mpi=pmix"
  EXECUTION="${HOST_MPIRUN_COMMAND} --cpu-bind=none --distribution=arbitrary --multi-prog ./run_file"
  ```
- **Placement des cœurs :** Les wrappers MPMD (`prog_${ExeNameOut}.sh`) utilisent `exec taskset -c <slot>-<slot+omp-1> ./script_${ExeNameOut}.ksh` et les scripts de composants utilisent `exec ./${ExeNameOut}` pour garantir un placement strict sans laisser de sous-shells orphelins.
- **Environnement :** Chargement dynamique et autonome de `env_toubkal.sh` dans chaque script de composant.

## 4. Résultats de Validation (Job `8235576`)
La simulation des 3 jours a été exécutée et validée avec succès le 2026-11-06 :
- **Job ID :** `8235576`
- **Statut Slurm :** `COMPLETED` (Code de sortie `0:0`)
- **Durée totale du job :** 35 minutes 53 secondes (début 13:28:18, fin 14:04:11)
- **Détail des 3 périodes au sein du même job :**
  - **Jour 1 (1980-01-01) :** Étape Slurm `8235576.0`, calcul modèle 14 min 48 s (`RealCpuTime = 888.33 s`), sortie `COMPLETED 0:0`.
  - **Jour 2 (1980-01-02) :** Étape Slurm `8235576.1`, calcul modèle 6 min 17 s (`RealCpuTime = 377.33 s`), sortie `COMPLETED 0:0`.
  - **Jour 3 (1980-01-03) :** Étape Slurm `8235576.2`, calcul modèle 6 min 15 s (`RealCpuTime = 375.27 s`), sortie `COMPLETED 0:0`.
- **État final `run.card` :** `CumulPeriod = 4`, `PeriodDateBegin = 1980-01-04`, `PeriodState = Completed`.

## 5. Sorties et Restarts Archivés (`IGCM_OUT/`)
Toutes les sorties scientifiques et les fichiers de restart ont été archivés sous :
`$IGCM_OUT/LMDZOR/TEST/amip/test-3days-HIST-r35/`

- **Restarts journaliers produits (validés via `ncdump -h`) :**
  - `ATM/Restart/test-3days-HIST-r35_1980010[1-3]_restart.nc` (414 Mo)
  - `ATM/Restart/test-3days-HIST-r35_1980010[1-3]_restartphy.nc` (1,6 Go)
  - `SRF/Restart/test-3days-HIST-r35_1980010[1-3]_sechiba_rest.nc` (565 Mo)
  - `SBG/Restart/test-3days-HIST-r35_1980010[1-3]_stomate_rest.nc` (1,5 Go)
- **Fichiers historiques NetCDF :**
  - `ATM/Output/DA/test-3days-HIST-r35_1980010[1-3]_1980010[1-3]_1D_histday.nc` (1,2 Go/jour)
  - `ATM/Output/HF/test-3days-HIST-r35_1980010[1-3]_1980010[1-3]_HF_histhf.nc` (856 Mo/jour)
  - `ATM/Output/MO/test-3days-HIST-r35_1980010[1-3]_1980010[1-3]_1M_histmth.nc` (1,7 Mo/jour)
  - Sorties ORCHIDEE / STOMATE haute fréquence et mensuelles correspondantes.
## 6. Distinction Explicite entre Validation Technique et Validation Scientifique
- **Validation Technique (Réussie) :**
  - Stabilité et exécution complète de la chaîne couplée LMDZ–ORCHIDEE–STOMATE–XIOS sur 28 nœuds Toubkal (111 rangs MPI : 110 rangs LMDZ/ORCHIDEE avec 14 threads OpenMP et 1 rang XIOS).
  - Résolution des dépendances dynamiques via `env_toubkal.sh` et finalisation Slurm native propre via `srun --mpi=pmix`.
  - Lecture des restarts 1979-12-31 via l'espace de stage et production des restarts journaliers intègres pour les 3 jours.
  - Archivage libIGCM complet sous `IGCM_OUT/`.
- **Validation Scientifique (Non couverte par un test de 3 jours) :**
  - Un banc d'essai de 3 jours confirme uniquement la validité opérationnelle et numérique du modèle.
  - Il ne permet pas d'évaluer la représentativité climatologique (précipitations, températures, flux de surface) ni les équilibres des réservoirs de carbone et d'eau. Une validation scientifique nécessite des intégrations pluriannuelles ou décennales.
