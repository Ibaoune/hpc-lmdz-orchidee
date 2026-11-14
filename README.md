# ORCHIDEE Offline on Toubkal HPC — Team User & Compilation Guide

> [!IMPORTANT]
> **Complete GCC 12 Reproduction Guide Available:**  
> For a fully self-contained, beginner-friendly reproduction workflow using **GCC 12.3.0 (`gompi-2023a`)** with all sources bundled directly in GitHub (no IPSL SVN checkouts or shared directory copies required), consult:  
> **[`docs/userguide/reproduction_guide.md`](docs/userguide/reproduction_guide.md)**

 M. El Aabaribaoune (Climate Assessment, Modeling, and Services, @UM6P)  

**Repository Location:** `~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl`  
**Shared Group:** `climat-um6p-st-iwri-7ksifkvwkuy-users`  

---

## Table of Contents
1. [Introduction & Purpose](#1-introduction--purpose)
2. [Quick Start: Validated 1-Day Offline Test (Feb 1, 1979)](#2-quick-start-validated-1-day-offline-test-feb-1-1979)
3. [Directory Structure & Components](#3-directory-structure--components)
4. [HPC Environment & Required Modules](#4-hpc-environment--required-modules)
5. [Toubkal Architecture Files (`arch-X64_TOUBKAL.*`)](#5-toubkal-architecture-files-arch-x64_toubkal)
6. [How to Compile or Recompile ANY ORCHIDEE Version](#6-how-to-compile-or-recompile-any-orchidee-version)
   - [Method 1: Automated Script (`compile_orchidee_toubkal.sh`)](#method-1-automated-script-compile_orchidee_toubkalsh)
   - [Method 2: Direct FCM Compilation (`makeorchidee_fcm`)](#method-2-direct-fcm-compilation-makeorchidee_fcm)
   - [How to switch to a new ORCHIDEE version (e.g. AgriVoltaics r4999)](#how-to-switch-to-a-new-orchidee-version-eg-agrivoltaics-r4999)
   - [Recompiling IOIPSL and XIOS](#recompiling-ioipsl-and-xios)
7. [Shared Datasets & Boundary Conditions on Toubkal](#7-shared-datasets--boundary-conditions-on-toubkal)
8. [How to Launch Simulations](#8-how-to-launch-simulations)
   - [Workflow A: Validated 1-Day Offline Test (Reproducible)](#workflow-a-validated-1-day-offline-test-reproducible)
   - [Workflow B: libIGCM Automated Workflow (Production / Spinup)](#workflow-b-libigcm-automated-workflow-production--spinup)
   - [Workflow C: Standalone Execution (Local Sites / Testing / AgriVoltaics)](#workflow-c-standalone-execution-local-sites--testing--agrivoltaics)
9. [Helper Scripts](#9-helper-scripts)
10. [Troubleshooting & Common Pitfalls](#10-troubleshooting--common-pitfalls)

---

## 1. Introduction & Purpose

This shared repository provides an operational, pre-configured, and fully validated build of the **IPSL ORCHIDEE Offline** land surface model and **XIOS 2** parallel I/O server for all researchers and colleagues in the **UM6P IWRI Climate Team**.

It is designed to:
1. Run standard offline simulations out of the box using existing verified binaries and configurations.
2. Reproduce the validated **1-day offline simulation (1979-02-01)** driven by historical LMDZ outputs in 5 simple commands.
3. Recompile the existing model with custom compiler flags or debug options.
4. Plug in and compile **any version** of ORCHIDEE (such as trunk, stable branches, or the AgriVoltaics `ORCHIDEE-AV` revision `r4999` from Lia) without fighting HPC library paths or linker errors.
5. Launch simulations either via **Standalone / Direct Slurm** (fast testing, regional zoom grids, AgriVoltaics) or via **libIGCM** (automated multi-year experiments).

---

## 2. Quick Start: Validated 1-Day Offline Test (Feb 1, 1979)

A reproducible reference test on the regional 35 km zoomed grid (384 × 221 points) over Morocco and North Africa has been executed and validated (Slurm Job **8228332**, 26 seconds runtime, exit code 0).

To reproduce this simulation from scratch:

```bash
# 1. Charger l'environnement Toubkal
source modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env

# 2. Extraire et formater le forçage LMDZ historique (si non présent)
python3 scripts/prepare_forcing_1day.py

# 3. Initialiser son répertoire personnel d'exécution (liens symboliques)
bash scripts/setup_offline_1day.sh /chemin/vers/mon_repertoire

# 4. Soumettre la simulation Slurm (8 cœurs MPI)
cd /chemin/vers/mon_repertoire
sbatch run_offline_1day.sh

# 5. Suivre l'exécution et vérifier les sorties NetCDF
squeue -u $USER
ls -lh /chemin/vers/mon_repertoire/*.nc
```

> [!TIP]
> Details on the parameters, datasets, and diagnostic checks are documented in [`configs/offline_1day/README.md`](configs/offline_1day/README.md) and [`scripts/README.md`](scripts/README.md).

---

## 3. Directory Structure & Components

The shared repository is organized as follows:

```
~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/
├── README.md                      # Guide général d'utilisation et de compilation IPSL
├── compile_orchidee_toubkal.sh    # Script de compilation automatisé
├── submit_standalone_template.sh  # Modèle générique de batch Slurm standalone
├── bin/                           # Exécutables opérationnels compilés :
│   ├── gcm_384x220x79_prod.e      # Binaire couplé LMDZ–ORCHIDEE 35 km (MPI/OpenMP)
│   ├── orchidee_ol                # Binaire ORCHIDEE Offline (XIOS 2 + OpenMPI)
│   ├── orchideedriver             # Driver standalone
│   ├── xios_server.exe            # Serveur d'E/S XIOS parallèle
│   └── xios_server_prod.exe       # Alias serveur d'E/S pour libIGCM
├── configs/
│   ├── offline_1day/              # Configuration validée pour le test ORCHIDEE Offline (1979-02-01)
│   │   ├── README.md              # Documentation détaillée du test validé
│   │   ├── run.def                # Paramétrage SECHIBA offline
│   │   └── *.xml                  # Définitions XML pour XIOS (iodef, context, field, file)
│   └── lmdzor/
│       └── r35/                   # Configurations couplées LMDZ–ORCHIDEE 35 km (Maroc)
│           ├── README.md          # Guide complet et synthèse scientifique
│           ├── common/            # Paramètres partagés (gcm.def, physiq.def, vert.def, CWRR, XIOS)
│           ├── experiments/       # Configurations de production (ETAT0, HIST, SSP245, SSP585)
│           └── tests/             # Bancs d'essai 3 jours préparés (HIST, SSP245, SSP585)
├── scripts/
│   ├── README.md                  # Documentation des outils et procédures
│   ├── prepare_forcing_1day.py    # Extraction et formatage du forçage LMDZ HF
│   ├── setup_offline_1day.sh      # Initialisation automatisée des répertoires d'exécution
│   └── run_offline_1day.sh        # Script Slurm de lancement opérationnel
├── simulations/
│   ├── offline_1day/              # Espace d'exécution historique (job 8228332, hors Git)
│   └── offline_1day_reproduction/ # Simulation de reproduction validée (job 8234089, hors Git)
├── modules/
│   └── foss-2018b/                # Modules EasyBuild démasqués pour Lmod Toubkal
├── lib/                           # Bibliothèques statiques (libioipsl.a, liborchidee.a)
├── modeles/
│   ├── IOIPSL/                    # Bibliothèque E/S IPSL (tag v2_2_5)
│   ├── XIOS/                      # Serveur d'E/S XML (XIOS 2, r2466)
│   └── ORCHIDEE/              # Code source du modèle de surface
│       ├── arch/                  # Fichiers d'architecture Toubkal (arch-X64_TOUBKAL.*)
│       ├── src_sechiba/       # Hydrologie et bilan hydrique/énergétique
│       ├── src_stomate/       # Dynamique de la végétation et cycle du carbone
│       └── makeorchidee_fcm   # Moteur de compilation FCM
├── config/
│   └── ORCHIDEE_OL/           # Configurations expérimentales OOL
└── libIGCM/                   # Gestionnaire de simulation multi-annuel IPSL
```

---

## 4. HPC Environment & Required Modules

The model stack on Toubkal is built upon the **foss-2018b** toolchain (GNU Compiler Collection 7.3.0, OpenMPI 3.1.1, NetCDF 4.6.1, NetCDF-Fortran 4.4.4, and HDF5 1.10.2).

> [!IMPORTANT]
> **Why custom modules directory?**  
> Toubkal cluster system administrators hide older EasyBuild modules by prepending a dot to their `.lua` filenames (e.g. `.7.3.0-2.30.lua`), causing standard `module load GCC/7.3.0-2.30` to fail with `Lmod: unknown module`.  
> An unhidden module symlink tree is maintained at:
> `~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modules/foss-2018b`

### How to Load the Environment
Before compiling any code or running simulations, run:
```bash
source ~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env
```

---

## 5. Toubkal Architecture Files (`arch-X64_TOUBKAL.*`)

The build configuration files are located in `modeles/ORCHIDEE/arch/`:

1. **`arch-X64_TOUBKAL.env`**:
   - Sets `MODULEPATH` to access unhidden modules.
   - Loads `foss/2018b`, `netCDF/4.6.1-foss-2018b`, and `netCDF-Fortran/4.4.4-foss-2018b`.
   - Exports `NETCDF_CONFIG`, `NETCDFF_CONFIG`, `HDF5_LIB_DIR`.
2. **`arch-X64_TOUBKAL.path`**:
   - Resolves include and library paths for NetCDF C, NetCDF Fortran, and HDF5.
3. **`arch-X64_TOUBKAL.fcm`**:
   - Compiler flags for GCC 7.3.0: `%COMPILER mpifort`, `-O3 -ffree-line-length-none`.
   - Linker options: explicitly links `%BASE_LD -lstdc++` (required for XIOS C++ runtime).

---

## 6. How to Compile or Recompile ANY ORCHIDEE Version

### Method 1: Automated Script (`compile_orchidee_toubkal.sh`)

```bash
cd ~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl
./compile_orchidee_toubkal.sh -clean -j 8
```

### Method 2: Direct FCM Compilation (`makeorchidee_fcm`)

```bash
source modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env
cd modeles/ORCHIDEE
./makeorchidee_fcm -arch X64_TOUBKAL -driver -parallel mpi -prod -j 8
```

The resulting executable is automatically copied to:
`bin/orchidee_ol`

---

## 7. Shared Datasets & Boundary Conditions on Toubkal

The shared IPSL input database is hosted directly on Lustre at:
```
/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/igcmg/IGCM/SRF/
```

### Essential Surface Datasets:
| Parameter / Field | Path / File on Toubkal |
| :--- | :--- |
| **Soil Texture & Hydraulics** | `/srv/lustre01/.../igcmg/IGCM/SRF/SOIL/soils_param.nc` |
| **PFT Vegetation Distribution** | `/srv/lustre01/.../igcmg/IGCM/SRF/PFTMAPS/CMIP6/ESA-LUH2v2/historical/15PFT.v1/` |
| **Surface Background Albedo** | `/srv/lustre01/.../igcmg/IGCM/SRF/albedo/alb_bg_modisopt_2D_ESA_v2.nc` |
| **Topography Slopes** | `/srv/lustre01/.../igcmg/IGCM/SRF/cartepente2d_15min.nc` |
| **Reference Temperature** | `/srv/lustre01/.../igcmg/IGCM/SRF/reftemp.nc` |
| **Hydrological River Routing** | `/srv/lustre01/.../igcmg/IGCM/SRF/ROUTING/routing.nc` |
| **Floodplains** | `/srv/lustre01/.../igcmg/IGCM/SRF/ROUTING/floodplains.nc` |

---

## 8. How to Launch Simulations

### Workflow A: Validated 1-Day Offline Test (Reproducible)
See [Section 2](#2-quick-start-validated-1-day-offline-test-feb-1-1979) and [`configs/offline_1day/README.md`](configs/offline_1day/README.md).

### Workflow B: libIGCM Automated Workflow (Production / Spinup)
For long-term multi-annual runs with automatic calendar cycling and archiving:
```bash
cd config/ORCHIDEE_OL/OOL_SEC
# Configure config.card and COMP/orchidee_ol.card
../../libIGCM/ins_job
sbatch Job_<EXPERIMENT>
```

### Workflow C: Standalone Execution (Local Sites / Testing / AgriVoltaics)
For localized runs, point tests, or custom sensitivity runs:
```bash
bash scripts/setup_offline_1day.sh /path/to/my_run_dir
# Edit run.def as needed
cd /path/to/my_run_dir
sbatch run_offline_1day.sh
# (ou sbatch scripts/run_offline_1day.sh /path/to/my_run_dir)
```

### Workflow D: Coupled Regional Simulations LMDZ–ORCHIDEE (35 km Zoom)
Pour exécuter ou tester les simulations couplées à l'échelle régionale sur le Maroc ($384 \times 220 \times 79$), utiliser les configurations structurées sous [`configs/lmdzor/r35/`](configs/lmdzor/r35/) :
- **Bancs d'essai validés (3 jours) :** voir [`configs/lmdzor/r35/tests/`](configs/lmdzor/r35/tests/) :
  - **HIST (1980-01-01 au 1980-01-03) :** Job `8235576` (`COMPLETED 0:0`, 35m53s).
  - **SSP2-4.5 (2015-01-01 au 2015-01-03) :** Job `8235853` (`COMPLETED 0:0`, 36m03s).
  - **SSP5-8.5 (2015-01-01 au 2015-01-03) :** Job `8235912` (`COMPLETED 0:0`, 36m03s).
- **Campagnes de production :** voir [`configs/lmdzor/r35/experiments/`](configs/lmdzor/r35/experiments/) (`etat0`, `hist`, `ssp245`, `ssp585`).
- **Recréation des liens de restarts (après clonage Git) :**
  ```bash
  mkdir -p simulations/lmdzor_restart_stage
  ln -sfn /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist simulations/lmdzor_restart_stage/sim35km-HIST
  ln -sfn /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist simulations/lmdzor_restart_stage/sim35km-HIST-after2002
  ```
Consultez le guide détaillé : [`configs/lmdzor/r35/README.md`](configs/lmdzor/r35/README.md).

---

## 9. Helper Scripts

1. **[`scripts/prepare_forcing_1day.py`](scripts/prepare_forcing_1day.py)**:
   - Extraction et conversion reproductibles du forçage atmosphérique HF LMDZ vers le format `orchidee_ol`.
2. **[`scripts/setup_offline_1day.sh`](scripts/setup_offline_1day.sh)**:
   - Création et vérification automatique des liens symboliques (exécutable, configurations, conditions aux limites, restart).
3. **[`scripts/run_offline_1day.sh`](scripts/run_offline_1day.sh)**:
   - Script de soumission Slurm optimisé pour Toubkal avec nettoyage préalable automatique et lanceur `mpirun`.
4. **[`compile_orchidee_toubkal.sh`](compile_orchidee_toubkal.sh)**:
   - Moteur de compilation automatisé pour Toubkal.
5. **[`submit_standalone_template.sh`](submit_standalone_template.sh)**:
   - Modèle générique de batch Slurm standalone.

---

## 10. Troubleshooting & Common Pitfalls

### Issue 1: "Lmod: module GCC/7.3.0-2.30 unknown"
- **Cause:** Toubkal system administrators hid legacy modules with leading dots.
- **Solution:** Always run `source .../arch-X64_TOUBKAL.env` which loads the unhidden module path `modules/foss-2018b`.

### Issue 2: "Impossible to use compression for variable in parallel mode"
- **Cause:** XIOS NetCDF-4 collective write mode (`par_access="collective"`) rejects `compression_level > 0`.
- **Solution:** Set `compression_level="0"` in `file_def_orchidee.xml`.

### Issue 3: "Domain domain_landpoints_out does not have overlapped regions"
- **Cause:** ORCHIDEE tries by default to interpolate unstructured/zoomed grids onto a regular CMIP6 144x143 grid (`XIOS_REMAP_OUTPUT = True`).
- **Solution:** Set `XIOS_REMAP_OUTPUT = n` in `run.def` to preserve the native regional zoom grid.

### Issue 4: Out-of-bounds crash in `routing_basins` during standalone runs
- **Cause:** On zoomed grids, on-the-fly basin calculation without a pre-computed `river_desc.nc` produces invalid indices.
- **Solution:** Set `RIVER_ROUTING = n` in `run.def` for standalone surface physics (hydrology, evaporation, temperature, and runoff remain fully active).

### Issue 5: OpenMPI abort with "OPAL ERROR: Not initialized in file pmix2x_client.c"
- **Cause:** OpenMPI 3.1.1 on Toubkal was built without Slurm PMI support, making direct `srun` fail.
- **Solution:** Use `mpirun -np ${SLURM_NTASKS:-8} ./orchidee_ol` in Slurm batch scripts.

### Issue 6: "The restart file already exists on the disc. IOIPSL will not overwrite it"
- **Cause:** IOIPSL `restopenout` refuses to overwrite existing restart files from an aborted run.
- **Solution:** In a personal directory, specify `--force` with `scripts/run_offline_1day.sh` to voluntarily purge known previous outputs before running. The shared reference directory `simulations/offline_1day` is strictly protected against purge.

---
