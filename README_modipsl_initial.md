# ORCHIDEE Offline on Toubkal HPC — Team User & Compilation Guide

**Author:** M. El Aabaribaoune (UM6P)  
**Affiliation:** Mohammed VI Polytechnic University (UM6P) — IWRI Climate Team  
**Contact:** `mohammad.elaabaribaoune@um6p.ma`  
**Host System:** Toubkal Supercomputer (UM6P HPC, Slurm / Lustre)  
**Repository Location:** `/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models`  
**Shared Group:** `climat-um6p-st-iwri-7ksifkvwkuy-users`  

---

## Table of Contents
1. [Introduction & Purpose](#1-introduction--purpose)
2. [Directory Structure & Components](#2-directory-structure--components)
3. [HPC Environment & Required Modules](#3-hpc-environment--required-modules)
4. [Toubkal Architecture Files (`arch-X64_TOUBKAL.*`)](#4-toubkal-architecture-files-arch-x64_toubkal)
5. [How to Compile or Recompile ANY ORCHIDEE Version](#5-how-to-compile-or-recompile-any-orchidee-version)
   - [Method 1: Automated Script (`compile_orchidee_toubkal.sh`)](#method-1-automated-script-compile_orchidee_toubkalsh)
   - [Method 2: Direct FCM Compilation (`makeorchidee_fcm`)](#method-2-direct-fcm-compilation-makeorchidee_fcm)
   - [How to switch to a new ORCHIDEE version (e.g. AgriVoltaics r4999)](#how-to-switch-to-a-new-orchidee-version-eg-agrivoltaics-r4999)
   - [Recompiling IOIPSL and XIOS](#recompiling-ioipsl-and-xios)
6. [Shared Datasets & Boundary Conditions on Toubkal](#6-shared-datasets--boundary-conditions-on-toubkal)
7. [How to Launch Simulations](#7-how-to-launch-simulations)
   - [Workflow A: libIGCM Automated Workflow (Production / Spinup)](#workflow-a-libigcm-automated-workflow-production--spinup)
   - [Workflow B: Standalone Execution (Local Sites / Testing / AgriVoltaics)](#workflow-b-standalone-execution-local-sites--testing--agrivoltaics)
8. [Helper Scripts Created by M. El Aabaribaoune](#8-helper-scripts-created-by-m-el-aabaribaoune)
9. [Troubleshooting & Common Pitfalls](#9-troubleshooting--common-pitfalls)
10. [Permissions and Team Collaboration](#10-permissions-and-team-collaboration)

---

## 1. Introduction & Purpose

This shared repository provides an operational, pre-configured, and fully validated build of the **IPSL ORCHIDEE Offline** land surface model and **XIOS 2** parallel I/O server for all researchers and colleagues in the **UM6P IWRI Climate Team**.

It is designed so that any colleague can:
1. Run standard offline simulations out of the box using existing verified binaries.
2. Recompile the existing model with custom compiler flags or debug options.
3. Plug in and compile **any version** of ORCHIDEE (such as trunk, stable branches, or the AgriVoltaics `ORCHIDEE-AV` revision `r4999` from Lia) without fighting HPC library paths or linker errors.
4. Launch simulations either via **libIGCM** (automated multi-year experiments) or **Standalone** (ideal for localized AgriVoltaics experiments and point simulations).

---

## 2. Directory Structure & Components

The shared repository is organized as follows:

```
/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/
├── README.md                      # This comprehensive guide (Author: M. El Aabaribaoune)
├── compile_orchidee_toubkal.sh    # One-click compilation script (Author: M. El Aabaribaoune)
├── submit_standalone_template.sh  # Slurm batch template for standalone runs (Author: M. El Aabaribaoune)
├── modules/
│   └── foss-2018b/                # Unhidden EasyBuild module definitions for Toubkal Lmod
└── modipsl/
    ├── bin/                       # Installed operational executables:
    │   ├── orchidee_ol            # Primary ORCHIDEE Offline binary (37 MB, XIOS 2 + OpenMPI)
    │   ├── orchideedriver         # Standalone driver binary
    │   └── xios_server.exe        # Asynchronous parallel XML I/O server (12 MB)
    ├── lib/                       # Compiled static libraries (libioipsl.a, liborchidee.a, etc.)
    ├── modeles/
    │   ├── IOIPSL/                # IPSL I/O library (tag v2_2_5)
    │   ├── XIOS/                  # XML I/O Server (XIOS 2, revision 2466)
    │   └── ORCHIDEE/              # Land surface model source code
    │       ├── arch/              # Toubkal architecture files (arch-X64_TOUBKAL.*)
    │       ├── src_sechiba/       # Energy and water balance hydrology
    │       ├── src_stomate/       # Carbon cycle and vegetation dynamics
    │       ├── src_xml/           # XIOS XML configuration templates
    │       └── makeorchidee_fcm   # FCM build command
    ├── config/
    │   └── ORCHIDEE_OL/           # Offline configurations and experiments
    │       └── OOL_SEC/           # SECHIBA offline simulation setup (config.card, PARAM/)
    ├── libIGCM/                   # IPSL simulation manager with Toubkal cluster drivers
    └── util/                      # Makefile generators (ins_make, AA_make.gdef)
```

---

## 3. HPC Environment & Required Modules

The model stack on Toubkal is built upon the **foss-2018b** toolchain (GNU Compiler Collection 7.3.0, OpenMPI 3.1.1, NetCDF 4.6.1, NetCDF-Fortran 4.4.4, and HDF5 1.10.2).

> [!IMPORTANT]
> **Why custom modules directory?**  
> Toubkal cluster system administrators hide older EasyBuild modules by prepending a dot to their `.lua` filenames (e.g. `.7.3.0-2.30.lua`), causing standard `module load GCC/7.3.0-2.30` to fail with `Lmod: unknown module`.  
> An unhidden module symlink tree is maintained at:
> `/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modules/foss-2018b`

### How to Load the Environment
Before compiling any code or running standalone simulations, simply run:
```bash
source /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env
```

This single command automatically executes:
```bash
module purge
module use /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modules/foss-2018b
module load Perl
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1-GCC-7.3.0-2.30
module load netCDF/4.6.1-foss-2018b
module load HDF5/1.10.2-foss-2018b
module load netCDF-Fortran/4.4.4-foss-2018b
```

---

## 4. Toubkal Architecture Files (`arch-X64_TOUBKAL.*`)

To ensure portability across any ORCHIDEE version on Toubkal, M. El Aabaribaoune created standardized architecture files located in:
- `modipsl/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env`
- `modipsl/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.path`
- `modipsl/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.fcm`
- `modipsl/util/AA_make.gdef` (section `#-Q- toubkal`)

### Key Compiler & Linker Flags:
1. **FCM Flags (`arch-X64_TOUBKAL.fcm`)**:
   - Compiler: `mpif90`
   - Base Fortran Flags: `-cpp -ffree-line-length-0 -fdefault-real-8 -DNC_DOUBLE`
   - Optimization: `-O3` (production)
   - Dynamic C++ Linking: `%BASE_LD -lstdc++` (required because XIOS 2 is written in C++ and linked into Fortran).
2. **Library Paths (`arch-X64_TOUBKAL.path`)**:
   - Explicitly links both NetCDF Fortran (`-lnetcdff`) and NetCDF C (`-lnetcdf`).
   - Links HDF5 runtime libraries (`-lhdf5_hl -lhdf5 -lz -lcurl`).

---

## 5. How to Compile or Recompile ANY ORCHIDEE Version

### Method 1: Automated Script (`compile_orchidee_toubkal.sh`)

M. El Aabaribaoune created a turnkey compilation script at the root of `shared/TEAM/models`:

```bash
cd /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models

# Compile using 8 parallel cores:
./compile_orchidee_toubkal.sh -j 8

# Or perform a complete clean rebuild:
./compile_orchidee_toubkal.sh -clean -j 16
```
This script handles loading the environment, checking prerequisite libraries (IOIPSL, XIOS), calling FCM, and installing `orchidee_ol` into `modipsl/bin/`.

---

### Method 2: Direct FCM Compilation (`makeorchidee_fcm`)

If you prefer to compile manually from within the `ORCHIDEE` source directory:

```bash
# 1. Load the Toubkal environment
source /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env

# 2. Navigate to ORCHIDEE source directory
cd /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles/ORCHIDEE

# 3. Launch compilation for offline driver with XIOS 2
./makeorchidee_fcm -xios2 -j 8 -parallel mpi -prod -arch X64_TOUBKAL -driver
```

**Compilation flags explained:**
- `-xios2` : Compile with XIOS 2 parallel I/O server client interface.
- `-j 8` : Compile in parallel using 8 tasks.
- `-parallel mpi` : Target distributed memory OpenMPI.
- `-prod` : High optimization level (`-O3`).
- `-arch X64_TOUBKAL` : Use Toubkal architecture definitions.
- `-driver` : Build the offline driver (`dim2_driver.f90` -> `orchidee_ol`).

The compiled binary will be placed in `modipsl/bin/orchidee_ol`.

---

### How to Switch to a New ORCHIDEE Version (e.g. AgriVoltaics r4999)

When receiving a new ORCHIDEE version (e.g. the AgriVoltaics AV package from Lia or a specific SVN revision):

1. **Backup current version:**
   ```bash
   cd /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles
   mv ORCHIDEE ORCHIDEE_r7086_quest_bak
   ```
2. **Place new ORCHIDEE version:**
   ```bash
   # Copy or extract new code into modeles/ORCHIDEE
   cp -r /path/to/ORCHIDEE_AV /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles/ORCHIDEE
   ```
3. **Copy the Toubkal Architecture Files into the new version:**
   ```bash
   cp ORCHIDEE_r7086_quest_bak/arch/arch-X64_TOUBKAL.* ORCHIDEE/arch/
   ```
4. **Compile:**
   ```bash
   /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/compile_orchidee_toubkal.sh -clean -j 8
   ```

---

### Recompiling IOIPSL and XIOS

In 99% of cases, `libioipsl.a` and `libxios.a` are already compiled and do not need to be touched. However, if starting from a fresh checkout:

- **To recompile IOIPSL:**
  ```bash
  source .../modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env
  cd /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles/IOIPSL/src
  gmake clean
  gmake -j 8
  ```
- **To recompile XIOS:**
  ```bash
  source .../modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env
  cd /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles/XIOS
  ./make_xios --arch X64_TOUBKAL --job 8
  ```

---

## 6. Shared Datasets & Boundary Conditions on Toubkal

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

### Atmospheric Meteorological Forcing:
Pre-existing forcing datasets (or your custom regional / site-level ERA5, WFDEI, or station data) should be linked in your run directory or `COMP/orchidee_ol.card`.

---

## 7. How to Launch Simulations

### Workflow A: libIGCM Automated Workflow (Production / Spinup)

`libIGCM` automatically manages calendar dates, file transfers, restarts, and Slurm job submissions.

#### Step 1: Go to Experiment Setup Directory
```bash
cd /home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/config/ORCHIDEE_OL/OOL_SEC
```

#### Step 2: Edit `config.card`
Configure your experiment name, simulation interval, and MPI core count:
```ini
[UserChoices]
JobName=SECHIBA01
TagName=OL2
SpaceName=PROD
ExperimentName=sechiba
CalendarType=noleap
DateBegin=1901-01-01
DateEnd=1901-12-31
PeriodLength=1Y

[Executable]
SRF= ("", "")
OOL= (orchidee_ol, orchidee_ol, 31MPI)
IOS= (xios_server.exe, xios.x, 1MPI)
```

#### Step 3: Check Forcing in `COMP/orchidee_ol.card`
Verify the line defining `forcing_file.nc` points to your meteorological NetCDF forcing file.

#### Step 4: Generate the Slurm Script
```bash
../../libIGCM/ins_job
```
This produces `Job_SECHIBA01`, pre-filled with the Toubkal Slurm partition and account:
```bash
#SBATCH --account=CLIMAT-UM6P-ST-IWRI-7KSIFKVWKUY-DEFAULT-CPU
#SBATCH --ntasks=32
#SBATCH --hint=nomultithread
```

#### Step 5: Submit & Monitor
```bash
# Submit
sbatch Job_SECHIBA01

# Monitor
squeue -u $USER
tail -f Script_Output_SECHIBA01.*
```

---

### Workflow B: Standalone Execution (Local Sites / Testing / AgriVoltaics)

For localized runs, point tests, or AgriVoltaics experiments where you want direct control over parameters and output files without `libIGCM`:

#### Step 1: Create a Personal Run Directory
```bash
RUN_DIR="/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/users/$USER/my_orchidee_run"
mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"
```

#### Step 2: Copy XML Control Files and `run.def`
```bash
MODELS="/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models"

# Copy XIOS XML definitions:
cp ${MODELS}/modipsl/modeles/ORCHIDEE/src_xml/iodef.xml .
cp ${MODELS}/modipsl/modeles/ORCHIDEE/src_xml/context_orchidee.xml .
cp ${MODELS}/modipsl/modeles/ORCHIDEE/src_xml/field_def_orchidee.xml .
cp ${MODELS}/modipsl/modeles/ORCHIDEE/src_xml/file_def_orchidee.xml .
cp ${MODELS}/modipsl/modeles/ORCHIDEE/src_xml/context_input_orchidee.xml .

# Copy parameter definition file:
cp ${MODELS}/modipsl/config/ORCHIDEE_OL/OOL_SEC/PARAM/run.def .

# Copy Slurm batch template:
cp ${MODELS}/submit_standalone_template.sh ./run_job.sh
```

#### Step 3: Adjust XML Server Mode
In `iodef.xml`:
- If running **without** a detached XIOS server (recommended for simple standalone runs): set `using_server="false"`.
- If running with detached `xios_server.exe`: set `using_server="true"`.

#### Step 4: Configure `run.def`
Edit `run.def` to specify:
- Meteorological forcing filename (`ATM_FILE = /path/to/meteo.nc`)
- Soil map (`SOIL_FILE = /srv/lustre01/.../soils_param.nc`)
- PFT map (`PFT_FILE = /srv/lustre01/.../PFTmap.nc`)
- `STOMATE_OK_STOMATE = y` (or `n` for Sechiba only)
- `TIME_STEP = 1800`

#### Step 5: Submit Job
```bash
sbatch run_job.sh
```

---

## 8. Helper Scripts Created by M. El Aabaribaoune

To make development and simulation seamless for team members, custom utilities are provided:

1. **[`compile_orchidee_toubkal.sh`](file:///home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/compile_orchidee_toubkal.sh)**:
   - Automated compilation engine for Toubkal.
   - Cleans, verifies dependencies, links NetCDF C/Fortran and C++ runtimes, builds `orchidee_ol`.
2. **[`submit_standalone_template.sh`](file:///home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/submit_standalone_template.sh)**:
   - Ready-to-use Slurm script pre-configured with account `CLIMAT-UM6P-ST-IWRI-7KSIFKVWKUY-DEFAULT-CPU`, proper module loading, and sanity checks.
3. **[`libIGCM_sys_toubkal.ksh`](file:///home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/libIGCM/libIGCM_sys/libIGCM_sys_toubkal.ksh)**:
   - Core system adapter for `libIGCM` on Toubkal, handling Lustre scratch storage and Slurm job queue management.

---

## 9. Troubleshooting & Common Pitfalls

### Issue 1: "Lmod: module GCC/7.3.0-2.30 unknown"
- **Cause:** Toubkal system administrators hid legacy modules with leading dots.
- **Solution:** Always run `source .../arch-X64_TOUBKAL.env` which loads the unhidden module path `/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modules/foss-2018b`.

### Issue 2: "Cannot find -lstdc++" or undefined XIOS symbols
- **Cause:** XIOS 2 is written in C++ and compiled with Boost. When linking the Fortran driver `dim2_driver.f90`, the GNU C++ runtime must be linked explicitly.
- **Solution:** Ensure `%BASE_LD -lstdc++` is present in `arch-X64_TOUBKAL.fcm` (already configured in this repository).

### Issue 3: "Fatal error: module format 14 vs 16 mismatch"
- **Cause:** Mixing Fortran `.mod` files compiled with GCC 7 (`format 14`) with newer GCC compilers (GCC 10/11 uses `format 16`).
- **Solution:** Keep using the `foss-2018b` toolchain as defined in `arch-X64_TOUBKAL.env`.

### Issue 4: Slurm Job Pending with "Invalid Account"
- **Cause:** Using default or expired Slurm project account.
- **Solution:** Always use `#SBATCH --account=CLIMAT-UM6P-ST-IWRI-7KSIFKVWKUY-DEFAULT-CPU`.

---

## 10. Permissions and Team Collaboration

- **POSIX Group:** All team members belong to `climat-um6p-st-iwri-7ksifkvwkuy-users`.
- **Inheritance (SGID):** All directories in `shared/TEAM/models` have the SGID bit enabled (`chmod g+s`) so new files automatically inherit the shared group.
- **Write Access:** Team members can create their own simulation branches or standalone runs. Please do not delete `modipsl/bin/` or shared reference libraries without consulting the team.
- **Questions & Maintenance:** Contact **M. El Aabaribaoune** (`mohammad.elaabaribaoune@um6p.ma`).
