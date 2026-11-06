# Bancs d'Essai Couplés 3 Jours (`tests/`)

Ce dossier contient les trois bancs d'essai de courte durée (3 jours simulés, `PeriodLength = 1D`, `PeriodNb = 3`) préparés pour valider le fonctionnement et la chaîne technique du modèle couplé LMDZ–ORCHIDEE 35 km sur le supercalculateur Toubkal (UM6P).

> [!IMPORTANT]
> **Validation Technique des Trois Tests Couplés (Toubkal) :**
> - **HIST (1980) :** Job `8235576` — **Validé techniquement** (3 jours consécutifs complétés, sorties NetCDF et restarts archivés).
> - **SSP2-4.5 (2015) :** Job `8235853` — **Validé techniquement** (3 jours consécutifs complétés, sorties NetCDF et restarts archivés).
> - **SSP5-8.5 (2015) :** Job `8235912` — **Validé techniquement** (3 jours consécutifs complétés, sorties NetCDF et restarts archivés).
>
> **Distinction Explicite entre Validation Technique et Validation Scientifique :**
> - **Validation Technique (Validée pour les 3 bancs d'essai) :** Vérifie la portabilité, la chaîne Slurm MPMD (`srun --mpi=pmix`), la stabilité numérique sur 111 rangs MPI (110 rangs LMDZ/ORCHIDEE avec 14 threads OpenMP par rang et 1 rang XIOS sur 28 nœuds Toubkal), l'enchaînement séquentiel sans blocage libIGCM, la transition correcte des restarts initiaux (1979-12-31 pour HIST, 2014-12-31 pour SSP) et la conformité des sorties NetCDF archivées sous `IGCM_OUT/`.
> - **Validation Scientifique (Non couverte par un banc d'essai de 3 jours) :** Un banc d'essai de 3 jours ne permet en aucun cas d'évaluer la représentativité climatologique, les cycles saisonniers, les bilans énergétiques ou hydrologiques, l'accumulation de carbone, ni la réponse climatique aux forçages anthropiques. La validation scientifique requiert des simulations pluriannuelles ou décennales analysées à l'aide de diagnostics climatiques dédiés (ERA5, flux FLUXNET, CMIP6).

---

## 1. Synthèse Comparative des Trois Bancs d'Essai

| Caractéristique | HIST (`hist_3days`) | SSP2-4.5 (`ssp245_3days`) | SSP5-8.5 (`ssp585_3days`) |
| :--- | :--- | :--- | :--- |
| **Objectif** | Validation couplage & dynamique 1980 | Transition HIST -> Scénario intermédiaire | Transition HIST -> Scénario extrême |
| **Période simulée** | 1980-01-01 au 1980-01-03 (3 jours) | 2015-01-01 au 2015-01-03 (3 jours) | 2015-01-01 au 2015-01-03 (3 jours) |
| **Job ID Slurm** | `8235576` | `8235853` | `8235912` |
| **Statut Slurm** | `COMPLETED` (ExitCode `0:0`) | `COMPLETED` (ExitCode `0:0`) | `COMPLETED` (ExitCode `0:0`) |
| **Durée totale** | 35 min 53 s | 36 min 03 s | 36 min 03 s |
| **Temps calcul J1 / J2 / J3** | 14 min 48 s / 6 min 17 s / 6 min 15 s | 14 min 43 s / 6 min 19 s / 6 min 17 s | 15 min 09 s / 6 min 10 s / 6 min 12 s |
| **Ressources HPC** | 28 nœuds (111 MPI, 14 OpenMP) | 28 nœuds (111 MPI, 14 OpenMP) | 28 nœuds (111 MPI, 14 OpenMP) |
| **Lanceur Slurm** | `srun --mpi=pmix --multi-prog` | `srun --mpi=pmix --multi-prog` | `srun --mpi=pmix --multi-prog` |
| **Restart initial** | 1979-12-31 (`sim35km-HIST`) | 2014-12-31 (`sim35km-HIST-after2002`) | 2014-12-31 (`sim35km-HIST-after2002`) |
| **Forçage Végétation** | Historique 1980 (`PFTmap`, `woodharvest`) | Nominal SSP245 (`scenarios/ssp245`) | Nominal SSP585 (`scenarios/ssp585`) |
| **Forçage GES / Aérosols** | CMIP6 HIST (1980) | CMIP6 SSP2-4.5 (2015) | CMIP6 SSP5-8.5 (2015) |
| **Statut Technique** | **VALIDÉ** | **VALIDÉ** | **VALIDÉ** |
| **Statut Scientifique** | Non validé (run long requis) | Non validé (run long requis) | Non validé (run long requis) |

---

## 2. Description Individuelle des Tests

- **[`hist_3days/`](hist_3days/)** : Test historique du **1980-01-01 au 1980-01-03** démarrant depuis les fichiers restart certifiés du 1979-12-31 (`sim35km-HIST`). Validé techniquement avec le Job `8235576`.
- **[`ssp245_3days/`](ssp245_3days/)** : Test scénario SSP2-4.5 nominal du **2015-01-01 au 2015-01-03** démarrant depuis le restart HIST du 2014-12-31 (`sim35km-HIST-after2002`) avec forçages réels SSP245. Validé techniquement avec le Job `8235853`.
- **[`ssp585_3days/`](ssp585_3days/)** : Test scénario SSP5-8.5 nominal du **2015-01-01 au 2015-01-03** démarrant depuis le restart HIST du 2014-12-31 (`sim35km-HIST-after2002`) avec forçages réels SSP585. Validé techniquement avec le Job `8235912`.

---

## 3. Environnement Logiciel & Versions Validées

- **Compilateur & MPI :** GCC 7.3.0, OpenMPI 3.1.1, NetCDF-C 4.6.1, NetCDF-Fortran 4.4.4, FFTW 3.3.8 (suite foss-2018b, chargés via `configs/lmdzor/r35/common/env_toubkal.sh`).
- **Noyaux modèles :**
  - LMDZ6 trunk révision 3855 (`IPSL-CM6A-MR`).
  - ORCHIDEE branch Quest révision 7086 (hydrologie CWRR 11 couches).
  - XIOS trunk révision 2466.
  - IOIPSL v2_2_5 révision 6058.
  - libIGCM révision 1572 adaptée Toubkal (`scripts/toubkal/libIGCM_sys_toubkal.ksh`).

---

## 4. Accès aux Restarts et Procédure après Clonage

libIGCM attend rigoureusement l'arborescence :
`${RestartPath}/${RestartJobName}/${Comp}/Restart/${RestartJobName}_${Date}_${File}`

Pour que les tests accèdent aux archives de restarts sans modifier les données originales ni dupliquer de gros volumes NetCDF, un espace de stage hors Git est utilisé sous `modipsl/simulations/lmdzor_restart_stage/` (exclu de Git par `.gitignore`).

### Procédure pour recréer les liens après un clonage du dépôt :
```bash
cd /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl

# 1. Créer le répertoire de stage
mkdir -p simulations/lmdzor_restart_stage

# 2. Créer les liens symboliques vers les archives certifiées
ln -sfn /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist simulations/lmdzor_restart_stage/sim35km-HIST
ln -sfn /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist simulations/lmdzor_restart_stage/sim35km-HIST-after2002
```

---

## 5. Emplacement des Sorties Archivées (`IGCM_OUT/`)

Les sorties et restarts des trois simulations de test sont archivés sous :
`/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/users/mohammad.elaabaribao/lmdz_simus/similar_dirs_jz_tbkl/gpfsscratch_tbkl/rech/NONE/mohammad.elaabaribao/IGCM_OUT/LMDZOR/TEST/amip/`

- **HIST 3 jours :** `test-3days-HIST-r35/`
- **SSP2-4.5 3 jours :** `test-3days-SSP245-r35/`
- **SSP5-8.5 3 jours :** `test-3days-SSP585-r35/`

Pour chaque simulation, les sous-répertoires contiennent :
- `ATM/Output/DA/` : fichiers journaliers atmosphériques `1D_histday.nc` (1,2 Go/jour).
- `ATM/Output/MO/` : fichiers mensuels atmosphériques `1M_histmth.nc`.
- `ATM/Restart/` : fichiers restart `restart.nc` (414 Mo) et `restartphy.nc` (1,6 Go) par jour.
- `SRF/Output/` & `SRF/Restart/` : sorties et restarts SECHIBA (`sechiba_rest.nc` 565 Mo/jour).
- `SBG/Output/` & `SBG/Restart/` : sorties et restarts STOMATE (`stomate_rest.nc` 1,5 Go/jour).

---

## 6. Procédure de Soumission Slurm
```bash
# Pour soumettre l'un des bancs d'essai sur Toubkal :
cd configs/lmdzor/r35/tests/<test_directory>
sbatch Job_test-3days-<EXPERIMENT>-r35
```
