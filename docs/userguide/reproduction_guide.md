# ORCHIDEE Offline on Toubkal HPC — Complete GCC 12 Reproduction Guide

**Author:** International Water Research Institute (IWRI), Mohammed VI Polytechnic University (UM6P)  
**Target Audience:** Researchers, climate modelers, and HPC users reproducing ORCHIDEE Offline on Toubkal HPC  
**Document Path:** `docs/userguide/reproduction_guide.md`  
**Toolchain:** GCC 12.3.0 / OpenMPI 4.1.5 (`gompi/2023a`), NetCDF-Fortran 4.6.1, HDF5 1.14.0  
**Benchmark Case:** Standalone 1-Day Simulation (February 1, 1979) on 35 km North Africa Zoom Grid  

---

## Executive Overview

This guide provides a self-contained, step-by-step walkthrough to reproduce the standalone (offline) land surface model **ORCHIDEE Offline** using modern **GCC 12.3.0** on the **Toubkal HPC supercomputer** at UM6P.

Unlike previous configurations that required external Subversion (SVN) checkouts or copying files from shared team directories, this repository is **completely self-contained**:
- All necessary source codes for **IOIPSL**, **XIOS 2**, and **ORCHIDEE Quest** are bundled directly inside the GitHub clone.
- All required build tools, including bundled **FCM distributions** (`FCM_OLD`, `FCM_V1.2`) and third-party C++ libraries (`Boost`, `Blitz`, `RapidXML`), are pre-packaged and ready to compile.
- Build automation scripts dynamically adapt paths to any user personal Lustre workspace without hardcoded paths.
- Input meteorological forcing and surface parameters on Toubkal HPC are referenced via lightweight symbolic links without duplicating large multi-gigabyte datasets.

---

## Step 1 — Clone the GitHub Repository

### 1. Objective
Establish an isolated personal workspace on the high-performance Lustre parallel filesystem and obtain the complete project source tree directly from GitHub.

### 2. Commands

```bash
# 1. Create and enter your personal directory on Lustre
export WORK_DIR="/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/users/$USER/orchidee_reproduction"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

# 2. Clone the dedicated Toubkal offline branch
git clone --branch orchidee-offline-toubkal --single-branch \
  https://github.com/Ibaoune/hpc-lmdz-orchidee.git modipsl

# 3. Enter the cloned repository and verify Git tracking
cd modipsl
git status -sb
```

> **Working Directory Policy on Toubkal HPC:**  
> All personal compilation and simulations must be conducted in your personal workspace under `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/users/$USER/`.  
> **Never modify or run simulations in the shared team directory** (`shared/TEAM/models/modipsl`).

### 3. Explanation
- **`--branch orchidee-offline-toubkal`:** Selects the branch pre-configured with GCC 12 architecture files and bundled model sources.
- **`--single-branch`:** Clones only the active branch, keeping the repository clone fast and lightweight.
- **`git status -sb`:** Confirms that the working tree cleanly tracks `origin/orchidee-offline-toubkal`.

### 4. Expected Result
The command completes with exit code 0 and reports:
```text
## orchidee-offline-toubkal...origin/orchidee-offline-toubkal
```

### 5. Status
- **Validated:** Clean repository clone verified on Toubkal HPC.

---

## Step 2 — Verify Bundled Model Sources

### 1. Objective
Confirm that all model components, build tools, and third-party libraries required for the GCC 12 build are already present inside the repository clone, without downloading from IPSL servers or copying from historical installations.

### 2. Commands

```bash
# Verify the presence of all three core components
ls -ld modeles/IOIPSL/src modeles/XIOS/src modeles/ORCHIDEE/src_driver

# Verify bundled FCM build tools
ls -l modeles/XIOS/tools/FCM_OLD/bin/fcm modeles/ORCHIDEE/tools/FCM_V1.2/bin/fcm

# Verify bundled third-party libraries in XIOS
ls -ld modeles/XIOS/extern/boost modeles/XIOS/extern/blitz modeles/XIOS/extern/rapidxml

# Verify that no legacy binary objects exist in the clean clone
find modeles/ -name "*.o" -o -name "*.mod" -o -name "*.a"
```

### 3. Explanation and Component Provenance

The repository bundles the exact scientific code revisions used for the certified Toubkal GCC 12 benchmark:

| Component | Provenance & Version | Role in the Modeling Chain | Key Bundled Directories |
| :--- | :--- | :--- | :--- |
| **IOIPSL** | `IOIPSL/tags/v2_2_5` (SVN r6058) | Input/Output library for NetCDF I/O, calendar handling, and spatial interpolation routines | `src/`, `tools/FCM_V1.2/`, `Makefile` |
| **XIOS** | `XIOS2/trunk` (SVN r2466) | High-performance asynchronous XML-IO-Server client/server parallel library | `src/`, `extern/` (Boost, Blitz, RapidXML), `tools/FCM_OLD/`, `make_xios` |
| **ORCHIDEE** | `ORCHIDEE_Quest` (SVN r7086) | Standalone land surface model (SECHIBA 11-layer soil hydrology, energy balance, STOMATE carbon) | `src_driver/`, `src_sechiba/`, `src_parallel/`, `tools/FCM_V1.2/`, `makeorchidee_fcm` |

- **No external downloads:** Users do not need access to upstream IPSL Subversion servers.
- **No file copying:** The obsolete utility `copy_model_sources.py` is no longer needed.
- **Relocatable:** All internal Makefiles and scripts dynamically discover paths relative to the clone root.

### 4. Expected Result
- The directory and file listings display existing files.
- The `find` command returns **no files** (confirming 0 legacy compiled objects).

### 5. Status
- **Validated:** Full source inventory verified in the cloned repository (~8,030 tracked source and configuration files, ~88 MB total volume).

---

## Step 3 — Prepare the GCC 12 Environment

### 1. Objective
Load the Toubkal HPC modern software stack (`gompi/2023a`), configuring GCC 12, OpenMPI, parallel HDF5, and NetCDF-Fortran.

### 2. Commands

```bash
# Load the pre-configured GCC 12 environment
source modeles/ORCHIDEE/arch/arch-X64_TOUBKAL_GCC12.env

# Inspect the loaded toolchain components
echo "Fortran Compiler : $(which mpif90)"
echo "GCC Version      : $(gcc --version | head -n 1)"
echo "GFortran Version : $(gfortran --version | head -n 1)"
echo "NetCDF-Fortran   : $(nf-config --version)"
echo "Parallel HDF5    : $(h5pcc -showconfig | grep -i "HDF5 Version" | awk "{print \$3}")"
```

### 3. Explanation
- **`module purge`:** Clears any pre-existing environment modules (e.g. legacy GCC 7 or Intel toolchains) to avoid compiler ABI mismatches.
- **`gompi/2023a`:** Provides GCC 12.3.0 and OpenMPI 4.1.5.
- **`netCDF-Fortran/4.6.1-gompi-2023a`:** Automatically pulls in NetCDF-C 4.9.2 and parallel HDF5 1.14.0.
- **`h5pcc` vs `h5cc`:** Using `h5pcc` verifies the parallel HDF5 build linked with OpenMPI, avoiding interference from any user-level Python/Anaconda environments.

### 4. Expected Result
```text
Fortran Compiler : /srv/software/easybuild/software/OpenMPI/4.1.5-GCC-12.3.0/bin/mpif90
GCC Version      : gcc (GCC) 12.3.0
GFortran Version : GNU Fortran (GCC) 12.3.0
NetCDF-Fortran   : netCDF-Fortran 4.6.1
Parallel HDF5    : 1.14.0
```

### 5. Status
- **Validated:** Environment tested and confirmed on Toubkal login and compute nodes.

---

## Step 4 — Compile ORCHIDEE Offline

### 1. Objective
Independently compile **IOIPSL**, **XIOS**, and **ORCHIDEE Offline** using the unified GCC 12 build driver `compile_orchidee_gcc12.sh`.

### 2. Commands

To compile using 8 parallel processes:

```bash
./compile_orchidee_gcc12.sh -clean -j 8
```

For submission as a Slurm compute batch job:

```bash
sbatch --job-name=build_orch --partition=compute --nodes=1 --ntasks=1 --cpus-per-task=8 --time=00:45:00 --wrap="./compile_orchidee_gcc12.sh -clean -j 8"
```

### 3. Explanation of Build Stages and Solved Challenges

The build workflow executes three sequential stages:

1. **IOIPSL Library (`lib/libioipsl.a`):**
   - **`MAKE_DIR` Dynamic Adaptation:** The script automatically updates `MAKE_DIR` in `modeles/IOIPSL/src/Makefile` to point to the current repository directory while preserving the Makefile directory safety verification.
   - **Serial Archive Generation:** IOIPSL compiles serially to prevent archive member race conditions (`ar -rs` collisions on parallel filesystems), ensuring all 11 object files (including `mathelp.o`) are reliably archived.
   - **GCC 12 Flag:** Uses `-fallow-argument-mismatch` to handle legacy Fortran argument signatures cleanly.

2. **XIOS Library (`modeles/XIOS/lib/libxios.a`):**
   - **Bundled FCM Execution:** Automatically invokes the pre-packaged `modeles/XIOS/tools/FCM_OLD/bin/fcm` tool.
   - **C++11 Standard:** Compiles XIOS C++ source files, bundled Boost, Blitz++, and RapidXML with `-std=c++11 -O3`.
   - **Parallel NetCDF Support:** Automatically creates `modeles/XIOS/.void_dir` to link against Toubkal external parallel NetCDF library.

3. **ORCHIDEE Standalone Driver (`bin/orchidee_ol_gcc12`):**
   - **FCM Compilation:** Calls `./makeorchidee_fcm -xios2 -j 8 -parallel mpi -prod -arch X64_TOUBKAL_GCC12 -driver`.
   - Compiles biophysical modules (`SECHIBA`), carbon dynamics (`STOMATE`), and the atmospheric offline driver (`dim2_driver`).
   - Links the final executable against `libioipsl.a`, `libxios.a`, OpenMPI, NetCDF-Fortran, NetCDF-C, and parallel HDF5.
   - Installs the final standalone binary to `bin/orchidee_ol_gcc12`.

### 4. Expected Result
The compilation script completes with:
```text
==========================================================
SUCCESS: Compilation completed successfully!
Executable : /.../modipsl/bin/orchidee_ol_gcc12
-rwxrwx--x+ 1 ... 39M ... /.../modipsl/bin/orchidee_ol_gcc12
==========================================================
```

### 5. Status
- **Validated:** Full clean build independently verified on Toubkal HPC (`libioipsl.a` 1.8 MB, `libxios.a` 49 MB, `orchidee_ol_gcc12` 39 MB).

---

## Step 5 — Prepare Benchmark Input Data

### 1. Objective
Establish a clean simulation run directory and create symbolic links to the required meteorological forcing, initial land surface restart, and boundary datasets stored on Toubkal HPC.

### 2. Commands

```bash
# Initialize the execution directory (simulations/offline_1day_gcc12)
bash scripts/setup_offline_1day_gcc12.sh simulations/offline_1day_gcc12

# Inspect the prepared simulation directory
cd simulations/offline_1day_gcc12
ls -l
```

### 3. Explanation of Input Files

ORCHIDEE Offline requires meteorological forcing and surface boundary conditions. These datasets are already hosted on Toubkal shared storage and linked into your personal workspace without copying:

| Local Link Name | Description | Source Path on Toubkal HPC |
| :--- | :--- | :--- |
| `forcing_19790201.nc` | Atmospheric forcing for Feb 1, 1979 (temperature, humidity, wind, rainfall, solar radiation) | `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/simulations/offline_1day/forcing_19790201.nc` |
| `sechiba_rest_in.nc` | Initial land surface restart state (soil moisture, temperature profiles on 1979-01-31) | `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist/SRF/Restart/sim35km-HIST_19790131_sechiba_rest.nc` |
| `PFTmap.nc` | 15 Plant Functional Types vegetation distribution map (ESA-CCI CMIP6) | `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/igcmg/IGCM/SRF/PFTMAPS/CMIP6/ESA-LUH2v2/historical/15PFT.v1/PFTmap_1979.nc` |
| `soils_param.nc` | Soil hydraulic properties, texture, and saturated conductivity (11 layers) | `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/igcmg/IGCM/SRF/SOIL/soils_param.nc` |
| `routing.nc` | River routing flow directions and basin catchments | `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/igcmg/IGCM/SRF/ROUTING/routing.nc` |
| `floodplains.nc` | Floodplain storage fractions | `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/igcmg/IGCM/SRF/ROUTING/floodplains.nc` |
| `cartepente2d_15min.nc`| Topographic slope distribution | `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/igcmg/IGCM/SRF/cartepente2d_15min.nc` |
| `alb_bg.nc` | Background surface albedo map | `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/igcmg/IGCM/SRF/albedo/alb_bg_modisopt_2D_ESA_v2.nc` |
| `reftemp.nc` | Reference surface temperature climatology | `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/igcmg/IGCM/SRF/reftemp.nc` |

In addition to datasets, `setup_offline_1day_gcc12.sh` links:
- The newly built executable: `orchidee_ol -> modipsl/bin/orchidee_ol_gcc12`
- The model namelist: `run.def`
- XIOS XML control files: `iodef.xml`, `context_orchidee.xml`, `context_input_orchidee.xml`, `field_def_orchidee.xml`, `file_def_orchidee.xml`
- The Slurm execution script: `run_offline_1day.sh`

### 4. Expected Result
The script confirms that all 17 target files and links are established:
```text
Succès : Répertoire simulations/offline_1day_gcc12 vérifié et prêt pour GCC 12.
Pour lancer la simulation :
  cd simulations/offline_1day_gcc12 && sbatch run_offline_1day.sh
```

### 5. Status
- **Validated:** All input links verified and tested.

---

## Step 6 — Run the One-Day Simulation

### 1. Objective
Submit and execute the 1-day benchmark simulation using Slurm on Toubkal HPC compute resources.

### 2. Commands

```bash
cd simulations/offline_1day_gcc12

# Submit the Slurm job (8 MPI tasks, 1 node)
sbatch run_offline_1day.sh

# Monitor job execution in the Slurm queue
squeue -u $USER
```

### 3. Explanation of Job Parameters
- **MPI Processes:** Runs across 8 MPI ranks (`--ntasks=8`) on 1 node (`--partition=compute`).
- **Memory & Threads:** `ulimit -s unlimited`, OpenMP disabled (`OMP_NUM_THREADS=1`) for pure MPI execution.
- **Runtime:** A 1-day simulation executes in approximately **25 to 30 seconds**.
- **Slurm Log:** Output is captured in `slurm-<JOBID>.out`, and any system errors in `slurm-<JOBID>.err`.

### 4. Expected Result
After ~30 seconds, `squeue -u $USER` shows the job has completed, and the end of `slurm-<JOBID>.out` displays:
```text
==========================================================
Fin normale de la simulation : Fri Nov 13 ... 2026
==========================================================
```
And `out_orchidee_0000` ends with:
```text
-------------------------------------------
------> CPU Time Global    14.8 s
------> Real Time Global   15.0 s
-------------------------------------------
END of dim2_driver
```

### 5. Status
- **Validated:** Executed successfully on Toubkal HPC (`slurm-8304516.out`, exit code 0, 0 bytes in `.err`).

---

## Step 7 — Validate the Results

### 1. Objective
Verify the generated NetCDF history and restart files, confirming that the simulation reproduced the certified scientific reference results within acceptable numerical tolerances.

### 2. Commands

```bash
cd simulations/offline_1day_gcc12

# 1. Verify that all 4 expected NetCDF simulation outputs were produced
ls -lh sechiba_history.nc sechiba_out_2.nc sechiba_rest_out.nc driver_rest_out.nc

# 2. Check the model namelist reproducibility
sha256sum used_run.def

# 3. Check simulation time bounds inside the output NetCDF
ncdump -h sechiba_history.nc | grep -E "time_counter =|time_origin"
```

### 3. Validation Metrics & Reference Comparison

The generated simulation outputs were rigorously verified against the certified GCC 12 reference simulation (`shared/TEAM/models/modipsl/simulations/offline_1day_gcc12`):

| Output File | Total Variables | Valid Physical Points Comparison | Scientific Validation Status |
| :--- | :---: | :---: | :--- |
| **`sechiba_rest_out.nc`** | 107 | 100% Identical | **BIT-IDENTICAL (100%)** |
| **`driver_rest_out.nc`** | 18 | 100% Identical | **BIT-IDENTICAL (100%)** |
| **`sechiba_out_2.nc`** | 21 | Max Abs Diff = `0.0` | **IDENTICAL (100% on valid data points)** |
| **`sechiba_history.nc`** | 120 | Max Abs Diff = `0.0` | **IDENTICAL (100% on valid data points)** |
| **`used_run.def`** | Namelist | SHA256: `c91727e43212c0e40dc2c669433661dfdbc7f722b9a87e8ea32f15e20887f627` | **BIT-IDENTICAL** |

- **Physical Variables:** Latent heat flux (`evap`), sensible heat flux (`qh`), net radiation (`fluxsens`), surface temperature, soil moisture per layer (`humrel`), drainage, and runoff match the certified GCC 12 reference simulation with **zero difference** on all valid grid points.
- **Tolerance:** `Max Relative Difference = 0.0` (exact numerical reproduction).

### 4. Expected Result
All NetCDF output files exist with matching file sizes:
- `sechiba_history.nc`: ~102 MB
- `sechiba_out_2.nc`: ~172 MB
- `sechiba_rest_out.nc`: ~589 MB
- `driver_rest_out.nc`: ~9.5 MB

### 5. Status
- **Validated:** Exact scientific reproducibility confirmed on November 13, 2026.

---

## Troubleshooting Guide

This section explains root causes and permanent solutions for issues previously encountered during development. The workflow in this guide already prevents these issues, but they are documented here for complete technical transparency.

### 1. IOIPSL Error: `You must invoke this Makefile from its directory`
- **Symptom:** IOIPSL compilation fails immediately with an error from GNU Make.
- **Root Cause:** Historical IPSL `AA_make` templates write an absolute path into `MAKE_DIR` during installation. When the code is cloned into a new workspace, `$(shell /bin/pwd)` does not match `$(MAKE_DIR)`.
- **Permanent Solution:** In `compile_orchidee_gcc12.sh`, `MAKE_DIR` is dynamically updated before compiling:
  ```bash
  sed -i "s|^MAKE_DIR = .*|MAKE_DIR = ${ROOT_MODIPSL}/modeles/IOIPSL/src|" "${IOIPSL_MAKEFILE}"
  ```
  This retains the directory safety check while allowing relocation to any workspace.

### 2. IOIPSL Undefined References (`__mathelp_MOD_mathop_r11`) during Link
- **Symptom:** Linking ORCHIDEE fails with undefined symbols belonging to `mathelp.f90` inside `libioipsl.a`.
- **Root Cause:** In GNU Make, archive member rules (`$(MODEL_LIB)(%.o): %.f90`) trigger multiple simultaneous `ar -rs` executions when run with `-j N`. On parallel filesystems (Lustre), concurrent `ar` calls overwrite the archive and drop object files.
- **Permanent Solution:** IOIPSL contains only 11 files and is now compiled with serial `gmake` (taking ~15-20 seconds total), ensuring every object file is cleanly added to `libioipsl.a`.

### 3. XIOS Error: `./make_xios: line 296: fcm: command not found`
- **Symptom:** XIOS compilation stops because the FCM build command is missing from `$PATH`.
- **Root Cause:** Previous copy utilities omitted `tools/FCM_OLD/bin` because the directory was named `bin/`.
- **Permanent Solution:** `modeles/XIOS/tools/FCM_OLD` is fully tracked in Git, and `compile_orchidee_gcc12.sh` verifies its presence and sets executable permissions before building.

### 4. HDF5 Version Reporting 1.10.4 instead of 1.14.0
- **Symptom:** Environment diagnostics report an older HDF5 version (1.10.4) despite loading `gompi/2023a`.
- **Root Cause:** `h5cc` resolved to a serial Anaconda/Miniconda Python path in `$PATH`.
- **Permanent Solution:** Use `h5pcc -showconfig`, which explicitly invokes the OpenMPI-parallel HDF5 1.14.0 provided by Toubkal EasyBuild module stack.

### 5. Stale Compilation Lock: `fcm.bld.lock`
- **Symptom:** ORCHIDEE compilation aborts with `ERROR: Active compilation lock detected`.
- **Root Cause:** A previous compilation was killed unexpectedly, leaving a lockfile behind.
- **Permanent Solution:** Run `./compile_orchidee_gcc12.sh -clean -j 8`. The `-clean` option verifies if the locking process is still alive and safely removes stale lockfiles.

---

## Summary of Reproduction Workflow

```bash
# 1. Clone repository
git clone --branch orchidee-offline-toubkal --single-branch \
  https://github.com/Ibaoune/hpc-lmdz-orchidee.git modipsl
cd modipsl

# 2. Compile all components with GCC 12
./compile_orchidee_gcc12.sh -clean -j 8

# 3. Setup benchmark run directory
bash scripts/setup_offline_1day_gcc12.sh simulations/offline_1day_gcc12

# 4. Submit simulation
cd simulations/offline_1day_gcc12
sbatch run_offline_1day.sh

# 5. Verify outputs
ls -lh sechiba_history.nc sechiba_out_2.nc sechiba_rest_out.nc driver_rest_out.nc
```
