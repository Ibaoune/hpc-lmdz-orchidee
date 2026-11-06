# Configurations LMDZ–ORCHIDEE 35 km sur Toubkal (LMDZOR_v6.2.2)

Ce répertoire rassemble les configurations scientifiques, paramètres et tests préparatoires pour le modèle couplé atmosphère–surface continentale **LMDZ–ORCHIDEE** sur la grappe **Toubkal** (HPC UM6P), à une résolution régionale zoomée de **35 km** ($384 	imes 220 	imes 79$) centrée sur le Maroc.

---

## 1. Organisation de l'arborescence

```
configs/lmdzor/r35/
├── README.md                      # Présentation générale et synthèse scientifique
├── common/                        # Paramètres physiques, dynamiques, XML et environnement partagés
│   ├── DRIVER/                    # Drivers libIGCM compatibles (lmdz, orchidee, stomate, xios)
│   ├── PARAM/                     # gcm.def, physiq.def, vert.def, orchidee.def, iodef.xml
│   ├── POST/                      # Fichiers de monitoring libIGCM
│   └── env_toubkal.sh             # Script d'environnement GCC 7.3.0 / OpenMPI 3.1.1 / NetCDF (foss-2018b) partagé
├── experiments/                   # Configurations de production pluriannuelles (1979–2100)
│   ├── etat0/                     # Création des états initiaux start.nc et limit.nc (ce0l.e)
│   ├── hist/                      # Simulation historique 1979–2014 (1M, forçages CMIP6)
│   ├── ssp245/                    # Scénario futur moyen 2015–2100 (1M, forçages SSP2-4.5)
│   └── ssp585/                    # Scénario futur extrême 2015–2100 (1M, forçages SSP5-8.5)
└── tests/                         # Bancs d'essai 3 jours (HIST, SSP245, SSP585 validés techniquement)
    ├── hist_3days/                # Test couplé 1980-01-01 -> 1980-01-03 (Job 8235576, COMPLETED)
    ├── ssp245_3days/              # Test couplé 2015-01-01 -> 2015-01-03 (Job 8235853, COMPLETED)
    └── ssp585_3days/              # Test couplé 2015-01-01 -> 2015-01-03 (Job 8235912, COMPLETED)
```

---

## 2. Provenance et Historique Git

- **Configurations de production (`experiments/`) :** importées depuis la branche historique `origin/master` (`Toubkal_Configs/`, commit Git `61124b2`), préservées dans leur état d'origine.
- **Paramètres physiques (`common/PARAM/`) :** importés depuis l'expérience de référence validée sur Toubkal `Sim35km-Maroc-T01/PARAM/`.
- **Drivers partagés (`common/DRIVER/`) :** importés depuis les sources de référence LMDZOR_v6.2.2.
- **Noyau modèle :** LMDZ6 trunk révision **3855** (`IPSL-CM6A-MR`), ORCHIDEE branch Quest révision **7086** (hydrologie CWRR 11 couches), XIOS trunk révision **2466**, IOIPSL v2_2_5 révision **6058**, libIGCM révision **1572** adaptée à Toubkal.

---

## 3. Points Scientifiques et Incohérences Documentées

### Incohérence n°1 : Cartes de végétation PFT dans SSP245
Dans la configuration originale `SSP245/COMP/orchidee.card`, les fichiers de conditions aux limites de surface :
```ini
List= (${R_IN}/SRF/PFTMAPS/CMIP6/ESA-LUH2v2/scenarios/ssp585/PFTmap_${year}.nc, PFTmap.nc)
```
pointaient vers le répertoire du scénario **`ssp585`** au lieu de `ssp245`. L'analyse des sorties d'exécution archivées sur Lustre (`Script_Output_sim35km-HIST-future-ssp245.*`) confirme que la simulation historique de SSP245 sur Toubkal a effectivement été exécutée avec les cartes de végétation de SSP585.
- La configuration de production dans `experiments/ssp245/` est conservée telle quelle pour l'intégrité historique.
- Le banc d'essai nominal `tests/ssp245_3days/` utilise les véritables forçages nominaux `/scenarios/ssp245/`.

### Incohérence n°2 : État réel des archives de scénarios sur Lustre
L'inventaire des fichiers de simulation sous `/srv/lustre01/.../ARCHIVE/Simulations/LMDZ/r35/` montre :
- **HIST (1979–2014) :** Intégralement archivé, 432 mois continus avec restarts mensuels du 1979-01-31 au 2014-12-31.
- **SSP245 (2015–2100) :** **Incomplet**. La série contient 974 mois. Elle présente un trou de 3 mois au début 2086 (`208601`, `208602`, `208603`), et s'arrête définitivement au **31 mai 2096** (55 mois manquants jusqu'à fin 2100).
- **SSP585 (2015–2100) :** Atteint le 2100-12-31 (1120 mois uniques), mais présente un trou d'un mois dans les restarts en août 2094 (`209408`), et un doublon à la date `20990131` marquant la transition entre segments de relance.

---

## 4. Exécutables et Environnement d'Exécution

Les simulations couplées LMDZ–ORCHIDEE 35 km requièrent :
- `modipsl/bin/gcm_384x220x79_prod.e` (binaire couplé MPI/OpenMP compilé avec `foss-2018b`).
- `modipsl/bin/xios_server_prod.exe` (serveur d'E/S parallèle XIOS 2).
- **Allocation Slurm :** 28 nœuds exclusifs (partition `compute`), 111 rangs MPI (110 rangs LMDZ/ORCHIDEE avec 14 threads OpenMP par rang + 1 rang serveur XIOS dédié, 4 rangs par nœud = 56 cœurs physiques).
- **Lanceur Slurm natif :** `srun --mpi=pmix --multi-prog` configuré dans `libIGCM_sys_toubkal.ksh`.
- **Environnement partagé :** `common/env_toubkal.sh` assurant le chargement de la suite foss-2018b (GCC 7.3.0, OpenMPI 3.1.1, NetCDF-C 4.6.1, NetCDF-Fortran 4.4.4, FFTW 3.3.8).

---

## 5. Synthèse des Validations des Bancs d'Essai (3 Jours)

Les trois tests de courte durée (3 périodes journalières) ont été validés techniquement avec succès sur Toubkal :
- **HIST (`hist_3days`) :** Job `8235576`, durée 35 min 53 s, code de sortie `0:0`.
- **SSP245 (`ssp245_3days`) :** Job `8235853`, durée 36 min 03 s, code de sortie `0:0`.
- **SSP585 (`ssp585_3days`) :** Job `8235912`, durée 36 min 03 s, code de sortie `0:0`.

Toutes les sorties et restarts sont archivés sous :
`$IGCM_OUT/LMDZOR/TEST/amip/test-3days-<EXPERIMENT>-r35/`

### Procédure pour recréer l'espace de liaison des restarts après clonage :
```bash
cd /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl
mkdir -p simulations/lmdzor_restart_stage
ln -sfn /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist simulations/lmdzor_restart_stage/sim35km-HIST
ln -sfn /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist simulations/lmdzor_restart_stage/sim35km-HIST-after2002
```

> [!NOTE]
> **Distinction Technique vs Scientifique :** Ces tests attestent de la robustesse technique et de la reproductibilité du couplage LMDZ–ORCHIDEE sur Toubkal. La validation scientifique des simulations climatiques nécessite des runs pluriannuels complets avec évaluation des biais moyens et des bilans d'énergie.
