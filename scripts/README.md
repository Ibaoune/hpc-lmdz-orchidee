# Scripts d'automatisation pour ORCHIDEE Offline

Ce répertoire contient l'ensemble des outils nécessaires pour préparer, initialiser et exécuter la simulation ORCHIDEE Offline d'une journée (1er février 1979) sur le cluster Toubkal, soit dans l'espace de référence partagé, soit dans un répertoire d'exécution personnel.

---

## 1. Description des scripts

### `prepare_forcing_1day.py`
- **Rôle :** Extraction et formatage du forçage météorologique au format attendu par `orchidee_ol` à partir des sorties atmosphériques haute fréquence (HF) de la simulation couplée LMDZ–ORCHIDEE historique.
- **Fichier source par défaut :**  
  `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist/ATM/Output/HF/sim35km-HIST_19790201_19790228_HF_histhf.nc`
- **Variables traitées et converties :**
  - Température : `t2m` (K) $\rightarrow$ `Tair` (K)
  - Humidité spécifique : `q2m` (kg/kg) $\rightarrow$ `Qair` (kg/kg)
  - Pression au sol : `psol` (Pa) $\rightarrow$ `PSurf` (Pa)
  - Module du vent : `wind10m` (m/s) $\rightarrow$ `Wind` (m/s)
  - Composantes du vent : `u10m`, `v10m` (m/s) $\rightarrow$ `Wind_E`, `Wind_N` (m/s)
  - Précipitation liquide : `max(0, precip - snow)` (kg/m²/s) $\rightarrow$ `Rainf` (kg/m²/s)
  - Précipitation solide : `snow` (kg/m²/s) $\rightarrow$ `Snowf` (kg/m²/s)
  - Rayonnement solaire incident : `SWdnSFC` (W/m²) $\rightarrow$ `SWdown` (W/m²)
  - Rayonnement infrarouge incident : `LWdnSFC` (W/m²) $\rightarrow$ `LWdown` (W/m²)
  - Fraction continentale : `contfracOR` (-) $\rightarrow$ `contfrac` (-)
  - Coordonnées 2D : ajout de `nav_lon` et `nav_lat` pour la compatibilité IOIPSL/ORCHIDEE.
- **Protection contre l'écrasement :** Si le fichier de forçage cible existe déjà, le script refuse de le remplacer par défaut. L'option `--force` est nécessaire pour forcer sa régénération.
- **Usage :**
  ```bash
  python3 scripts/prepare_forcing_1day.py [--input FICHIER] [--output FICHIER] [--force]
  ```

---

### `setup_offline_1day.sh`
- **Rôle :** Configure de manière non destructive un répertoire d'exécution en créant les liens symboliques vers les configurations de référence, les jeux de données partagés, le restart initial et l'exécutable compilé.
- **Logique non destructive :**
  - Crée uniquement les liens manquants.
  - Conserve les liens déjà corrects.
  - Ne remplace et n'écrase **jamais** un fichier personnel existant (ex. un `run.def` modifié pour une étude de sensibilité).
  - Signale clairement tout conflit ou fichier préservé dans la console.
- **Usage :**
  ```bash
  # Configuration dans le répertoire par défaut partagé (simulations/offline_1day)
  bash scripts/setup_offline_1day.sh

  # Configuration dans un répertoire personnel de travail
  bash scripts/setup_offline_1day.sh ~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/users/$USER/mon_test_offline
  ```

---

### `run_offline_1day.sh`
- **Rôle :** Script de soumission Slurm pour Toubkal (8 tâches MPI, 1 nœud, partition `compute`, compte `CLIMAT-UM6P-ST-IWRI-7KSIFKVWKUY-DEFAULT-CPU`).
- **Emplacement des logs Slurm :**
  Les fichiers `slurm-<JOBID>.out` et `slurm-<JOBID>.err` sont générés par Slurm dans le répertoire depuis lequel la commande `sbatch` est invoquée. Pour que les logs soient créés directement dans votre répertoire personnel, lancez `cd <MON_REPERTOIRE> && sbatch run_offline_1day.sh`.
- **Gestion dynamique du répertoire d'exécution :**
  1. Argument en ligne de commande : `sbatch scripts/run_offline_1day.sh /chemin/vers/mon_repertoire`
  2. Variable d'environnement `SIM_DIR` : `SIM_DIR=/chemin sbatch scripts/run_offline_1day.sh`
  3. Répertoire courant `$PWD` si `run.def` et `orchidee_ol` y sont présents : `cd /mon/repertoire && sbatch run_offline_1day.sh`
- **Sécurité et protection des résultats :**
  - **Aucun nettoyage automatique :** Sans l'option explicite `--force`, aucun fichier n'est supprimé.
  - **Refus de lancement si résultats présents :** Si des sorties d'une exécution précédente existent (`sechiba_rest_out.nc`, `sechiba_history.nc`, etc.), le lancement est immédiatement refusé pour empêcher toute corruption ou écrasement.
  - **Protection absolue du répertoire de référence :** Toute exécution ou tentative de purge dans `simulations/offline_1day/` est strictement bloquée, avec ou sans l'option `--force`, afin de garantir l'intégrité des résultats certifiés (Job 8228332). Les simulations doivent utiliser des répertoires d'exécution personnels.
  - **Relance volontaire dans un répertoire personnel :** Dans un dossier personnel distinct de la référence, la relance est possible via l'option `--force`. La suppression est alors strictement limitée aux fichiers de résultats connus.
- **Usage :**
  ```bash
  # Lancement standard (refusé si des résultats existent déjà)
  sbatch scripts/run_offline_1day.sh [DIR]

  # Relance volontaire dans un répertoire personnel uniquement :
  sbatch scripts/run_offline_1day.sh [DIR] --force
  ```

---

## 2. Procédure de reproduction pour un collègue

### Option A : Dans un espace personnel indépendant (recommandé)

```bash
# 1. Charger l'environnement Toubkal
source modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env

# 2. Préparer le forçage météorologique (si non déjà généré)
python3 scripts/prepare_forcing_1day.py

# 3. Initialiser son répertoire personnel (non destructif)
MY_RUN_DIR="$HOME/lustre/climat-um6p-st-iwri-7ksifkvwkuy/users/$USER/offline_test_1day"
bash scripts/setup_offline_1day.sh "${MY_RUN_DIR}"

# 4. Lancer la simulation Slurm depuis son répertoire
cd "${MY_RUN_DIR}"
sbatch run_offline_1day.sh

# 5. Suivre l'exécution et vérifier les sorties
squeue -u $USER
tail -f slurm-*.out
ls -lh *.nc
```

### Option B : Consultation du test validé de référence

Le répertoire partagé `simulations/offline_1day/` conserve les résultats et preuves du test validé du 1er février 1979 (Job 8228332). Il est protégé contre toute nouvelle exécution ou écrasement.
