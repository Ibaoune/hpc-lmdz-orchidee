# ORCHIDEE Offline on Toubkal HPC — Step-by-Step Reproduction Guide (GCC 12)

**Author:** Mohammad El Aabaribaoune  
**Affiliation:** International Water Research Institute (IWRI), Mohammed VI Polytechnic University (UM6P), Ben Guerir, Morocco  
**Document:** `docs/userguide/reproduction_guide.md`  
**Target Audience:** New HPC users and researchers reproducing ORCHIDEE Offline on Toubkal HPC  

---

## Overview

This guide provides a progressive, step-by-step walkthrough to reproduce the standalone (offline) land surface model **ORCHIDEE** using the modern toolchain **GCC 12.3.0 / OpenMPI 4.1.5 (`gompi-2023a`)** on the **Toubkal HPC supercomputer** at UM6P.

The instructions are designed to be pedagogical, modular, and clear for beginners.

---

## Step 1 — Retrieving the Project from GitHub

### 1. Objective
Obtain a clean copy of the project scripts, Toubkal architecture definitions, and experiment configurations from the dedicated repository branch, while establishing an isolated personal workspace.

### 2. Commands

First, navigate to your personal workspace on the high-performance **Lustre parallel filesystem**:

```bash
# Set up and enter your personal directory on Lustre
export WORK_DIR="/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/users/$USER/orchidee_reproduction_gcc12_manual"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

# Clone the Toubkal offline branch into a local modipsl directory
git clone --branch orchidee-offline-toubkal --single-branch \
https://github.com/Ibaoune/hpc-lmdz-orchidee.git modipsl

# Enter the cloned repository and verify branch tracking
cd modipsl
git status -sb
```

> **Important Working Directory Rule:**  
> Every user must perform their work within their own personal directory (on Lustre under `users/$USER/` or alternatively in `$HOME`).  
> **Never work directly inside the shared team installation** (`shared/TEAM/models/modipsl`). Working in an isolated space protects certified reference simulations from unintended modifications and prevents conflicts between colleagues.

### 3. Explanation

- **`--branch orchidee-offline-toubkal`:** Selects the specific development branch designated for Toubkal HPC offline runs.
- **`--single-branch`:** Clones only the specified branch, avoiding downloading unnecessary commit history from other branches and keeping the download fast and lightweight.
- **`git status -sb`:** Confirms that your local working copy is on branch `orchidee-offline-toubkal` and cleanly tracks `origin/orchidee-offline-toubkal`.

#### What the Cloned Workspace Provides
The workspace provides the complete structural and operational framework for the reproduction:

> **Repository Maintenance vs. User Reproduction:**  
> The 1-day benchmark configurations in `configs/offline_1day/` are already versioned on GitHub. The GCC 12 build automation scripts, architecture files, staging tools, and reproduction documentation have been developed and validated locally, and are staged for publication to branch `orchidee-offline-toubkal`.  
> Users reproducing the experiment (such as Khadija) follow the operational steps in this guide (`git clone`, staging, compile, run) and do not need to run repository maintenance commands (`git add`, `git commit`).
- **Benchmark Experiment Configurations (versioned in Git):** 1-day offline benchmark XML descriptors and namelists in `configs/offline_1day/`.
- **GCC 12 Toolchain & Scripts (integrated for this reproduction):**
  - Build automation driver (`compile_orchidee_gcc12.sh`).
  - IOIPSL GCC 12 patch (`scripts/patches/ioipsl_makefile_toubkal.patch`).
  - Toubkal HPC GCC 12 architecture files (`modeles/ORCHIDEE/arch/` and `modeles/XIOS/arch/`).
  - Source code staging utility (`scripts/copy_model_sources.py`).
  - Slurm simulation drivers (`scripts/setup_offline_1day_gcc12.sh`, `scripts/run_offline_1day_gcc12.sh`).
  - Step-by-step user documentation (`docs/userguide/reproduction_guide.md`).

#### What is NOT Included in the Git Clone
To keep the Git repository lightweight and responsive, several components are intentionally excluded via `.gitignore`:
1. **Model Source Code Trees:** Full Fortran/C++ sources for ORCHIDEE, XIOS, and IOIPSL are not stored in Git. In the IPSL ecosystem, these large codebases are maintained in Subversion (SVN) repositories.
2. **NetCDF Datasets:** Git does not store large binary data. Atmospheric forcing files (`forcing_19790201.nc`), initial land surface restarts (`sechiba_rest_in.nc`), and multi-layer surface parameter maps (soil types, PFT vegetation maps, albedo) are not part of the repository.

### 4. Expected Result
- A `modipsl/` directory is created inside your personal workspace:
  `/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/users/$USER/orchidee_reproduction_gcc12_manual/modipsl`
- Running `git status -sb` returns:
  ```text
  ## orchidee-offline-toubkal...origin/orchidee-offline-toubkal
  ```

### 5. Status
- **Validated:** The GitHub project clone was successfully performed and verified on Toubkal HPC on **November 12, 2026**.

---

## Step 2 — Retrieving Model Source Trees (IOIPSL, XIOS, ORCHIDEE)

### 1. Objective
Safely populate the clean source code trees for **IOIPSL**, **XIOS**, and **ORCHIDEE** from Toubkal's shared team repository into your personal clone, while excluding all pre-existing binary compilation artifacts and strictly preserving existing GCC 12 architecture files.

---

### 2. Strategy: Local Sources vs. Shared NetCDF Datasets

An examination of the IPSL modeling tools confirmed that **IOIPSL, XIOS, and ORCHIDEE perform "in-source" builds**:
- The build tools generate object files (`*.o`), compiled module files (`*.mod`), configuration files, and temporary symlinks directly inside their own source directories.
- Compiling directly against a shared read-only directory is therefore not supported.

To balance user independence with storage efficiency, we make a clear distinction:

1. **Why Source Code is Copied Locally (~78 MB):**  
   Each user needs full write permissions in their own `modeles/` tree to compile with GCC 12 without interfering with the shared installation or other users. The total source code is very lightweight (~78 MB), making a local staging fast, clean, and completely independent.
2. **Why NetCDF Datasets Use Symbolic Links (Zero Data Duplication):**  
   Unlike source code, simulation inputs (meteorological forcing files, soil/vegetation parameters, and initial restart files) are static, read-only binary NetCDF files totaling many gigabytes. These files are never modified by the model during execution. Referencing them via lightweight symbolic links to the shared storage saves substantial disk space and avoids unnecessary file transfers.
3. **Why Pre-Existing Binaries Are Excluded:**  
   The shared installation contains legacy `.o`, `.mod`, and `.a` files compiled with GCC 7.3.0. Because GFortran module ABIs are strictly incompatible across major GCC versions, excluding all legacy binaries ensures a 100% clean baseline for GCC 12.

---

### 3. Why Model Sources Are Staged Separately from Git
After running `git clone`, the local directories `modeles/IOIPSL`, `modeles/XIOS`, and `modeles/ORCHIDEE` contain only the versioned architecture files (`arch/`), but no Fortran or C++ source files.

Sources are retrieved separately for three practical reasons:
1. **Upstream Source Size & SVN Distribution:** Upstream IPSL models comprise hundreds of thousands of lines of Fortran and C++ code and are officially distributed via Subversion (SVN) servers at IPSL, not Git. Storing these large external codebases in Git would cause unnecessary repository bloat.
2. **Network Isolation on Toubkal HPC:** Compute and login nodes on Toubkal operate in a restricted network environment without direct external internet access to IPSL SVN servers.
3. **Certified Shared Reference:** A complete, certified reference version of the source trees (ORCHIDEE Quest r7086, XIOS 2 r2466, IOIPSL v2_2_5) is already available on Toubkal's shared parallel storage at:
   ```text
   /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles/
   ```
   Users can stage clean copies of these source trees directly into their personal workspaces without needing external downloads.

---

### 4. The Three Retrieved Components

| Component | Role in the Modeling Chain | Key Subdirectories Copied |
| :--- | :--- | :--- |
| **IOIPSL** | Input/Output library for IPSL climate models. Handles NetCDF file operations, calendar management, and spatial interpolation routines. | `src/`, `tools/`, `Makefile` |
| **XIOS** (v2) | XML-IO-Server. High-performance asynchronous parallel I/O server and client library designed for massively parallel climate simulations. | `src/`, `extern/`, `inputs/`, `make_xios` |
| **ORCHIDEE** (Quest) | Standalone land surface model. Simulates vegetation dynamics, surface energy balance, soil hydrology, and carbon cycles. | `src_parallel/`, `src_driver/`, `src_sechiba/`, `src_parameters/`, `arch/` |

---

### 5. Role of the Staging Script: `scripts/copy_model_sources.py`

To eliminate the risk of manual copy errors or accidental overwrites, the repository provides an automated, non-destructive staging utility:
[`scripts/copy_model_sources.py`](../../scripts/copy_model_sources.py).

#### Key Design and Safety Mechanisms:
1. **Targeted Traversal:** Scans only the three required components (`IOIPSL`, `XIOS`, `ORCHIDEE`).
2. **Selective Filtering:** Filters out compiled objects (`*.o`, `*.mod`, `*.a`), build directories (`build/`, `obj/`, `ppsrc/`, `lib/`, `inc/`, `bin/`), version control metadata (`.svn/`), and NetCDF datasets (`*.nc`).
3. **Local Architecture Preservation:** If a destination file already exists in your local workspace (notably your Toubkal GCC 12 architecture files in `modeles/*/arch/`), the script **strictly preserves it** and skips copying from the shared space.
4. **Strictly Read-Only:** Accesses the shared team directory strictly in read-only mode, guaranteeing that reference models cannot be altered.
5. **Linear and Readable:** Written in standard Python with clear variable names and inline explanations, making it easily maintainable by researchers.

---

### 6. Difference Between `--dry-run` and Real Execution

| Mode | Command | What It Does | Modifies Disk? |
| :--- | :--- | :--- | :--- |
| **Preview (`--dry-run`)** | `python3 scripts/copy_model_sources.py --dry-run` | Traverses directories, evaluates filters, checks for existing files, and prints a detailed summary table. | **No** (0 bytes written) |
| **Real Execution** | `python3 scripts/copy_model_sources.py` | Performs the actual file copy of the identified clean source files into your local `modeles/` directory tree. | **Yes** (staged into personal clone) |

---

### 7. Execution and Actual Results

#### A. Preview (Dry-Run Mode)
The preview was first executed to inspect all planned operations without modifying any files:

```bash
python3 scripts/copy_model_sources.py --dry-run
```

#### B. Real Source Copy Execution — Completed Successfully
The actual source staging was executed from the root of your personal `modipsl/` clone:

```bash
python3 scripts/copy_model_sources.py
```

**Actual Execution Summary:**
```text
==============================================================================
 SUMMARY OF ACTIONS (EXECUTE (ACTUAL COPY))
==============================================================================
Component    | Files to Copy  | Volume (MB)  | Preserved  | Excluded Files
------------------------------------------------------------------------------
IOIPSL       |             61 |         0.76 |          3 |              0
XIOS         |           7308 |        71.55 |          6 |             20
ORCHIDEE     |            371 |         5.77 |          6 |              2
------------------------------------------------------------------------------
TOTAL        |           7740 |        78.09 |         15 |             22

Excluded Directories Count : 51
Shared Directory State     : UNTOUCHED (Strictly Read-Only)

--> COPY COMPLETED SUCCESSFULLY.
    Model source trees are now staged and ready for GCC 12 compilation.
==============================================================================
```

---

### 8. Verification of Staged Sources

Following execution, the staging was confirmed with the following verification checks:

#### 1. Confirm Source Directory Presence
```bash
ls -ld modeles/IOIPSL/src modeles/XIOS/src modeles/ORCHIDEE/src_parallel
```
*Result:* All three core directories are populated with clean Fortran/C++ sources.

#### 2. Confirm Complete Absence of Legacy Compilation Artifacts
```bash
find modeles/ -name "*.o" -o -name "*.mod" -o -name "*.a"
```
*Result:* **Empty** (0 legacy binary files present).

#### 3. Confirm Preservation of GCC 12 Architecture Files
```bash
ls -l modeles/ORCHIDEE/arch/arch-X64_TOUBKAL_GCC12.* modeles/XIOS/arch/arch-X64_TOUBKAL_GCC12.*
```
*Result:* All 6 GCC 12 architecture files (`.env`, `.fcm`, `.path`) are present intact with their original timestamps.

---

### 9. Validation Status

- **Step 1 (Git Clone):** **VALIDATED** (November 12, 2026).
- **Step 2 (Model Source Retrieval & Staging):** **VALIDATED** (Successfully executed on Toubkal HPC on November 13, 2026: 7,740 source files staged, 78.09 MB volume, 15 local architecture files preserved, 51 directories excluded, shared reference repository untouched).
