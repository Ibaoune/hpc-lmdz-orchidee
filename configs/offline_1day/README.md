# Configuration Validée : Test ORCHIDEE Offline d'une journée (1er Février 1979)

Ce répertoire contient l'ensemble des fichiers de paramétrage (`run.def` et fichiers XML XIOS) configurés et validés pour exécuter avec succès une simulation ORCHIDEE Offline d'une journée sur le cluster Toubkal (Job Slurm validé : `8228332`).

Cette configuration réutilise les sorties atmosphériques haute fréquence et le fichier de restart issus des simulations historiques couplées LMDZ–ORCHIDEE à 35 km de résolution sur le Maroc et l'Afrique du Nord.

---

## 1. Provenance des données et conditions aux limites

### Forçage météorologique atmosphérique
- **Source :** Simulation couplée historique LMDZ–ORCHIDEE (`sim35km-HIST`)
- **Fichier d'origine :**  
  `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist/ATM/Output/HF/sim35km-HIST_19790201_19790228_HF_histhf.nc`
- **Fréquence :** Pas de 3 heures (8 pas sur 24 heures : 01h30, 04h30, 07h30, 10h30, 13h30, 16h30, 19h30, 22h30).
- **Format cible :** `forcing_19790201.nc` préparé via `scripts/prepare_forcing_1day.py` (variables `Tair`, `Qair`, `PSurf`, `Wind`, `Wind_E`, `Wind_N`, `Rainf`, `Snowf`, `SWdown`, `LWdown`, `contfrac`).

### Restart SECHIBA initial
- **Source :** Fin du mois de janvier 1979 de la simulation historique couplée
- **Fichier d'origine :**  
  `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist/SRF/Restart/sim35km-HIST_19790131_sechiba_rest.nc`  
- **Mode d'utilisation :** Strictement en lecture seule via le lien symbolique `sechiba_rest_in.nc`.
- **Compatibilité :** Grille $384 \times 221$, 15 PFTs, schéma hydrologique CWRR 11 couches.

### Données de surface IGCM (Toubkal)
Toutes les données de surface proviennent de la base de référence IGCM :
- Texture des sols : `igcmg/IGCM/SRF/SOIL/soils_param.nc`
- Température de référence : `igcmg/IGCM/SRF/reftemp.nc`
- Pentes topographiques : `igcmg/IGCM/SRF/cartepente2d_15min.nc`
- Albédo de fond MODIS : `igcmg/IGCM/SRF/albedo/alb_bg_modisopt_2D_ESA_v2.nc`
- Carte de végétation historique 1979 (15 PFTs) : `igcmg/IGCM/SRF/PFTMAPS/CMIP6/ESA-LUH2v2/historical/15PFT.v1/PFTmap_1979.nc`
- Réseau hydrographique : `igcmg/IGCM/SRF/ROUTING/routing.nc` et `floodplains.nc`

---

## 2. Paramétrage clé du modèle (`run.def`)

Le fichier `run.def` dérive de la configuration standard `OOL_SEC` (SECHIBA offline sans composante biochimique STOMATE) avec les réglages spécifiques suivants :

```ini
# Durée de la simulation : 1 jour (86 400 s)
TIME_LENGTH = 1D
TIME_SKIP = 0
DRIVER_reset_time = y
TIME_STEP = 1800

# Fichiers de forçage et restart
FORCING_FILE = forcing_19790201.nc
HEIGHT_LEV1 = 2.0
HEIGHT_LEVW = 10.0
STOMATE_OK_STOMATE = n
RESTART_FILEIN = NONE
SECHIBA_restart_in = sechiba_rest_in.nc
SECHIBA_rest_out = sechiba_rest_out.nc

# Diagnostics XIOS
XIOS_ORCHIDEE_OK = y
SECHIBA_HISTFILE2 = n
WRITE_STEP = 0
WRITE_STEP2 = 0
PRINTLEV = 2

# Physique et hydrologie
HYDROL_CWRR = y
RIVER_ROUTING = n
RIVER_DESC = n
SOILTYPE_CLASSIF = zobler
ROUGH_DYN = y
OK_FREEZE = y
OK_EXPLICITSNOW = y
READ_REFTEMP = y
DEPTH_MAX_T = 90
ALB_BG_MODIS = y
NVM = 15

# Désactivation du rééchantillonnage XIOS vers grille régulière
XIOS_REMAP_OUTPUT = n
```

---

## 3. Corrections critiques appliquées pour le succès du test

Pour garantir une exécution reproductible sans modifier le code source Fortran d'ORCHIDEE :

1. **Désactivation du rééchantillonnage XIOS (`XIOS_REMAP_OUTPUT = n`) :**  
   Par défaut, ORCHIDEE bascule les grilles non cartésiennes vers un domaine CMIP6 $144 \times 143$ (`domain_landpoints_regular`). Ce mécanisme échoue sur une grille régionale zoomée. En forçant `XIOS_REMAP_OUTPUT = n`, les sorties restent sur la grille native $384 \times 221$ à 35 km.
2. **Désactivation du routage fluvial externe (`RIVER_ROUTING = n`) :**  
   L'initialisation dynamique du routage fluvial (`routing_basins`) sans fichier précalculé `river_desc.nc` sur grille zoomée générait des index nuls hors bornes (`mask_coast_glo(0)`). En mode standalone SECHIBA, `RIVER_ROUTING = n` désactive le routage fluvial vers l'océan sans impacter la physique de surface (ruissellement, drainage, humidité du sol, bilan thermique restent intégralement calculés).
3. **Suppression de la compression NetCDF4 parallèle (`compression_level="0"`) :**  
   Dans `file_def_orchidee.xml`, la compression `deflate` avec accès collectif parallèle (`par_access="collective"`) est rejetée par NetCDF-4/HDF5 (`Impossible to use compression for variable 'lat' when using parallel mode`). Le niveau de compression a été fixé à `0`.
4. **Ajustement du niveau de diagnostic (`output_level="5"`) :**  
   Permet d'extraire toutes les variables physiques de surface (température, eau du sol, flux turbulents, ruissellement, drainage) tout en évitant de demander à XIOS des champs de routage inactifs.
5. **Lanceur OpenMPI (`mpirun` vs `srun`) :**  
   OpenMPI 3.1.1 sur Toubkal n'intègre pas le support PMI direct de Slurm (`srun` échoue avec `OPAL ERROR: Not initialized in file pmix2x_client.c`). L'appel doit utiliser `mpirun -np ${SLURM_NTASKS:-8} ./orchidee_ol`.
6. **Protection des résultats et gestion des relances :**  
   Le répertoire de référence partagé `simulations/offline_1day/` est quant à lui strictement protégé contre toute exécution ou purge.

---

## 4. Résultats et preuves de validation (Job 8228332)

- **Job Slurm :** `8228332` sur `slurm-compute-h21c8-u5-svn4`
- **Durée réelle :** 26 secondes (48 pas de temps, temps CPU modèle ~15,8 s)
- **État final :** `COMPLETED` (ExitCode 0:0)
- **Validation des variables physiques (37 206 points continentaux actifs) :**
  - Température de surface (`temp_sol`) : $-50,78\text{ }^\circ\text{C}$ à $+38,44\text{ }^\circ\text{C}$ (moyenne $+8,47\text{ }^\circ\text{C}$, 0 NaN).
  - Humidité totale du sol (`mrso`) : $319,69$ à $839,44\text{ kg/m}^2$ (moyenne $536,12\text{ kg/m}^2$, 0 NaN).
  - Humidité de surface (`mrsos`) : $7,14$ à $39,50\text{ kg/m}^2$ (moyenne $18,70\text{ kg/m}^2$, 0 NaN).
  - Évapotranspiration (`evap`) : $-1,19$ à $+6,66\text{ mm/j}$ (moyenne $+1,05\text{ mm/j}$, 0 NaN).
  - Ruissellement (`runoff`) : $0,00$ à $149,33\text{ mm/j}$ (moyenne $0,19\text{ mm/j}$, 0 NaN).
  - Drainage profond (`drainage`) : $0,00$ à $39,12\text{ mm/j}$ (moyenne $1,08\text{ mm/j}$, 0 NaN).
  - Restart produit (`sechiba_rest_out.nc`) : 589,4 Mo, horodaté à $t = 86400\text{ s}$ (`time_steps = 48`), prêt pour le jour suivant.

> [!NOTE]
> **Distinction historique vs reproduction :**
> - **Test historique initial :** Réalisé le 4 novembre 2026 sous le Job Slurm `8228332` (preuves conservées dans `simulations/offline_1day/` : `slurm-8228332.out` et forçage `forcing_19790201.nc`).
> - **Reproduction certifiée :** Exécutée le 5 novembre 2026 sous le Job Slurm `8234089` dans le répertoire dédié `simulations/offline_1day_reproduction/` (COMPLETED, 27 s, 48 pas de temps, 24h simulées). L'ensemble complet des sorties NetCDF (`sechiba_history.nc`, `sechiba_out_2.nc`, `sechiba_rest_out.nc`, `driver_rest_out.nc`, `used_run.def`) y est intégralement conservé avec son manifeste `SHA256SUMS` et sa documentation `REPRODUCTION_NOTE.md`.

---

## 5. Distinction entre configuration de référence et répertoire d'exécution

Pour assurer la reproductibilité et la pérennité du travail d'équipe, l'architecture sépare rigoureusement la configuration des espaces d'exécution :

1. **Configuration de référence (`configs/offline_1day/`) :**
   - Contient les fichiers modèles canoniques sous contrôle de version Git (`run.def`, `iodef.xml`, `context_orchidee.xml`, etc.).
   - Ce répertoire ne doit **jamais** servir de dossier d'exécution direct afin d'éviter la génération de fichiers temporaires ou de sorties dans l'espace de configuration.
   - Les modifications dans ce répertoire doivent être réservées aux évolutions pérennes de la configuration de référence commune.

2. **Répertoire d'exécution personnel (ex. `simulations/offline_1day/` ou `~/users/$USER/mon_test`) :**
   - Espace de calcul instancié automatiquement par la commande `bash scripts/setup_offline_1day.sh [TARGET_DIR]`.
   - Tous les fichiers de configuration et jeux de données y sont liés symboliquement par défaut.
   - Les fichiers de diagnostic (`sechiba_history.nc`, `sechiba_out_2.nc`), restarts (`sechiba_rest_out.nc`) et logs y sont générés sans impacter les collègues.
   - **Études de sensibilité :** Pour tester une modification de paramètre sans altérer la référence partagée, il suffit de supprimer le lien symbolique local `run.def` dans son propre répertoire et de le remplacer par une copie modifiable :
     ```bash
     cd ~/users/$USER/mon_test
     rm run.def
     cp /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/configs/offline_1day/run.def ./run.def
     # Éditer run.def puis relancer avec sbatch run_offline_1day.sh
     ```
