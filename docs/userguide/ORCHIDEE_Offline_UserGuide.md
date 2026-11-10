# ORCHIDEE Offline on Toubkal HPC — User Guide
**Practical Handbook for Atmospheric Forcing Preparation, Standalone Execution, Sensitivity Studies, and AgriPV Modeling**

 M. El Aabaribaoune (Climate Assessment, Modeling, and Services, @UM6P)  

**Document Version:** `v0.1 — November 2026`  
**Repository Path:** `~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl`  
**Benchmark Reference:** Jobs `8228332` / `8234089` (Validated 1979-02-01 Simulation)

---

## Executive Summary & Target Audience

This user guide provides researchers, hydrologists, climate scientists, and agronomists with an operational, step-by-step technical manual for configuring, executing, and analyzing **ORCHIDEE Offline** simulations on the **Toubkal HPC supercomputer** at UM6P.

Unlike computationally intensive coupled ocean-atmosphere-land simulations, ORCHIDEE in standalone (offline) mode runs uncoupled from an atmospheric general circulation model, driven entirely by prescribed meteorological boundary conditions. This architecture makes it exceptionally fast (simulating one complete day in **~26 seconds** on 8 MPI cores), making it the tool of choice for:
- Evaluating alternative atmospheric forcing datasets (e.g., ERA5, AgMERRA, dynamically downscaled WRF).
- Assessing hydrological impacts (soil moisture, evapotranspiration, deep drainage, and surface runoff).
- Investigating Agrivoltaic (**AgriPV**) microclimates where solar panels redistribute precipitation and shade the underlying crop canopy.
- Performing parameter sensitivity experiments without modifying the shared reference repository.

---

## Table of Contents

- [Chapter 1 — Introduction to ORCHIDEE Offline](#chapter-1--introduction-to-orchidee-offline)
- [Chapter 2 — Toubkal HPC Environment & Repository Architecture](#chapter-2--toubkal-hpc-environment--repository-architecture)
- [Chapter 3 — Reproducing the Certified Reference Simulation](#chapter-3--reproducing-the-certified-reference-simulation)
- [Chapter 4 — Running ORCHIDEE With Custom Meteorological Forcing](#chapter-4--running-orchidee-with-custom-meteorological-forcing)
  - [4.1 Use Case A: Comparing Alternative Forcing Datasets](#41-use-case-a-comparing-alternative-forcing-datasets)
  - [4.2 Use Case B: Multi-Year & Multi-File Forcing Sequences](#42-use-case-b-multi-year--multi-file-forcing-sequences)
  - [4.3 Technical NetCDF Forcing Specifications & Constraints](#43-technical-netcdf-forcing-specifications--constraints)
  - [4.4 Adapting Meteorological Forcing for Agrivoltaic (AgriPV) Studies](#44-adapting-meteorological-forcing-for-agrivoltaic-agripv-studies)
- [Chapter 5 — Configuring and Submitting an Experiment](#chapter-5--configuring-and-submitting-an-experiment)
- [Chapter 6 — Model Outputs and Basic Diagnostics](#chapter-6--model-outputs-and-basic-diagnostics)
- [Chapter 7 — Comprehensive Troubleshooting Guide](#chapter-7--comprehensive-troubleshooting-guide)
- [Chapter 8 — Quick Start Checklist & Pre-Flight Reference](#chapter-8--quick-start-checklist--pre-flight-reference)

---

## Chapter 1 — Introduction to ORCHIDEE Offline

### 1.1 What is ORCHIDEE?
**ORCHIDEE** (*Organising Carbon and Hydrology In Dynamic Ecosystems*) is the state-of-the-art land surface model (LSM) developed by the Institut Pierre-Simon Laplace (IPSL). It simulates energy, water, carbon, and momentum exchanges between the terrestrial biosphere and the atmosphere. 

ORCHIDEE is structured into two core modular components:
1. **SECHIBA** (*Schématisation des Échanges Hydriques à l'Interface Biosphère-Atmosphère*): Solves the biophysical processes at short time steps (30 minutes), including the 11-layer **CWRR** (*Center for Water Resources Research*) Richards-based unsaturated soil hydrology, canopy water interception, explicit 3-layer snow thermodynamics, sensible and latent heat flux partitioning, and surface temperature prognostic equations.
2. **STOMATE** (*Saclay-Toulouse-Orsay Model for the Analysis of Terrestrial Ecosystems*): Solves carbon allocation, plant phenology, litter decomposition, soil organic carbon dynamics, and vegetation demographic transitions across 15 Plant Functional Types (**PFTs**) at daily to multi-annual time steps.

### 1.2 What is ORCHIDEE Offline?
In coupled Earth System simulations (such as LMDZ–ORCHIDEE), the atmospheric model (LMDZ) computes atmospheric dynamics and physics, passing near-surface temperature, humidity, wind, and precipitation to ORCHIDEE, which in turn feeds back surface albedo, temperature, and turbulent heat fluxes.

In **ORCHIDEE Offline** mode (`orchidee_ol`), the atmospheric model is replaced by an atmospheric driver (`orchideedriver.f90`). The model reads pre-computed or observed meteorological variables directly from external NetCDF files. Because no atmospheric Navier-Stokes equations are solved, offline simulations run several orders of magnitude faster than coupled simulations while preserving the complete biophysical and hydrological fidelity of the land surface scheme.

| Simulation Mode | Atmosphere Interaction | Computational Cost | Typical Focus |
| :--- | :--- | :--- | :--- |
| **Coupled (LMDZ–ORCHIDEE)** | Interactive 2-way dynamic feedback | Heavy (e.g. 28 Toubkal nodes, 111 MPI ranks, 36 min/3 days) | Global/regional climate projections, land-atmosphere feedbacks |
| **Offline (`orchidee_ol`)** | Prescribed 1-way atmospheric driving | Ultra-light (1 node, 8 MPI tasks, ~26 seconds/day) | Hydrological impacts, AgriPV microclimates, agronomy, parameter calibration |

### 1.3 Typical Applications of Standalone Simulations
- **Regional Hydrology:** Simulating soil moisture profiles down to 2 meters depth, evaporation, transpiration, deep percolation, and surface runoff across Morocco's agricultural basins.
- **Microclimate & Agrivoltaics (AgriPV):** Evaluating the agronomic impact of elevated photovoltaic panels on crop water stress, soil shading, and rainfall redistribution.
- **Dataset Benchmarking:** Quantifying surface water balance differences when forcing the same soil-vegetation system with ERA5 vs. high-resolution local weather stations.

### 1.4 Conceptual Workflow Diagram

```
+-------------------------------------------------------------------------+
|                  Prescribed Meteorological Forcing                      |
| (Tair, Qair, PSurf, Wind, Rainf, Snowf, SWdown, LWdown @ 3h or 30-min)  |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
|                ORCHIDEE Standalone Engine (orchidee_ol)                 |
|  +--------------------------------+  +-------------------------------+  |
|  |     SECHIBA Biophysics         |  |      STOMATE Biochemistry     |  |
|  | - 11-layer CWRR Hydrology      |  | - 15 Plant Functional Types   |  |
|  | - Canopy Interception & Energy |  | - Carbon Pools & Phenology    |  |
|  | - Explicit Multi-Layer Snow    |  | - Dynamic Allocation & Litter |  |
|  +--------------------------------+  +-------------------------------+  |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
|               Land Surface Diagnostics & State Outputs                  |
|  - Evapotranspiration (evap), Surface Temperature (temp_sol)            |
|  - Soil Moisture (mrso, mrsos), Total Runoff (mrro), Water Storage (TWS)|
|  - Net Radiation (netrad), Restart NetCDF Files (sechiba_rest_out.nc)   |
+-------------------------------------------------------------------------+
```

---

## Chapter 2 — Toubkal HPC Environment & Repository Architecture

### 2.1 Toubkal HPC Architecture
The **Toubkal** supercomputer at UM6P comprises Dell EMC PowerEdge compute nodes equipped with dual-socket Intel Xeon Platinum (Cascade Lake) processors delivering 56 physical cores per node, high-speed Mellanox HDR InfiniBand interconnects, and a high-performance Lustre parallel filesystem (`/srv/lustre01`).

### 2.2 Shared Repository Organization
All validated models, executables, scripts, and configurations are centrally maintained under:
```
~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl
```

Key subdirectories:
```
modipsl/
├── bin/
│   └── orchidee_ol             # Compiled standalone ORCHIDEE executable
├── configs/
│   └── offline_1day/           # Reference configuration files (run.def, XMLs)
├── scripts/
│   ├── prepare_forcing_1day.py # Atmospheric forcing extraction tool
│   ├── setup_offline_1day.sh   # Non-destructive directory setup utility
│   └── run_offline_1day.sh     # Slurm batch submission script
├── simulations/
│   ├── offline_1day/           # Reference benchmark workspace (Read-only)
│   └── offline_1day_reproduction/ # Certified reproduction workspace
├── modules/
│   └── foss-2018b/             # Unhidden module tree on Toubkal
└── modeles/
    └── ORCHIDEE/arch/          # Architecture and environment definitions
```

### 2.3 Loading the HPC Environment
The model is compiled with the **foss-2018b** toolchain:
- **Compilers:** GCC 7.3.0
- **MPI:** OpenMPI 3.1.1
- **I/O & Math:** NetCDF-C 4.6.1, NetCDF-Fortran 4.4.4, HDF5 1.10.2

> [!IMPORTANT]
> **Why custom module path?**  
> Toubkal system administrators hide older EasyBuild modules by prepending a dot to their `.lua` filenames (e.g. `.7.3.0-2.30.lua`), causing standard `module load GCC/7.3.0-2.30` to fail with `Lmod: unknown module`.  
> An unhidden module symlink tree is maintained in `modules/foss-2018b/`.

To load the environment correctly in any shell session:
```bash
source ~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env
```

### 2.4 Separation of Configuration, Run, and Output Spaces
To ensure strict reproducibility and protect reference data:
1. **Reference Configuration (`configs/offline_1day/`):** Contains the official XML and `run.def` templates tracked in Git. It must **never** be used as an active execution folder.
2. **Execution Directory (`simulations/` or personal scratch):** Created dynamically. All configuration and dataset files are linked symbollically. Outputs and logs are written locally without polluting the repository.
3. **Reference Results (`simulations/offline_1day/`):** Strictly write-protected against purge or rerun.

---

## Chapter 3 — Reproducing the Certified Reference Simulation

### 3.1 Overview of the Reference Benchmark
The reference experiment validates SECHIBA biophysics for **February 1, 1979** across the regional Morocco 35 km grid ($384 	imes 221$, 37,206 active land points) driven by atmospheric outputs from the coupled historical run.
- **Job ID:** `8228332` (initial validation) / `8234089` (certified reproduction).
- **Execution Time:** ~26 seconds walltime on 8 MPI cores (1 node).
- **Time Stepping:** 48 timesteps (dt = 1800 s = 30 minutes).

### 3.2 Step-by-Step Reproduction Procedure

#### Step 1: Navigate to Project Root
```bash
cd ~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl
```

#### Step 2: Source the Toubkal Architecture Environment
```bash
source modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env
```

#### Step 3: Instantiate a Personal Run Directory
Run the automated non-destructive setup script, supplying the path to your desired workspace (e.g., in your scratch directory):
```bash
MY_RUN_DIR="$HOME/lustre/climat-um6p-st-iwri-7ksifkvwkuy/users/$USER/test_offline_reprod"
bash scripts/setup_offline_1day.sh "${MY_RUN_DIR}"
```
This script automatically creates symlinks to:
- The compiled executable (`orchidee_ol`).
- Surface datasets (`soils_param.nc`, `PFTmap.nc`, `reftemp.nc`, `cartepente2d_15min.nc`, `alb_bg.nc`).
- Initial SECHIBA restart (`sechiba_rest_in.nc`).
- Atmospheric forcing (`forcing_19790201.nc`).
- XML definitions and `run.def`.

#### Step 4: Submit the Slurm Job
Navigate into the run directory and launch the Slurm batch script:
```bash
cd "${MY_RUN_DIR}"
sbatch run_offline_1day.sh
```

#### Step 5: Monitor Execution
```bash
squeue -u $USER
tail -f slurm-*.out
```

#### Step 6: Verify Benchmark Success
Check that all output NetCDF files were generated:
```bash
ls -lh sechiba_history.nc sechiba_out_2.nc sechiba_rest_out.nc driver_rest_out.nc
```
Expected execution metrics:
- Slurm job state: `COMPLETED` with exit code `0:0`.
- Execution walltime: 25 to 28 seconds.
- Physical domain sanity check: 0 NaNs across all 37,206 active continental grid points.
  - Mean surface temperature (`temp_sol`): $+8.47^\circ	ext{C}$ (range: $-50.8^\circ	ext{C}$ to $+38.4^\circ	ext{C}$).
  - Mean total soil moisture (`mrso`): $536.1	ext{ kg/m}^2$ (range: $319.7$ to $839.4	ext{ kg/m}^2$).
  - Mean evapotranspiration (`evap`): $1.05	ext{ mm/day}$ (range: $-1.2$ to $+6.7	ext{ mm/day}$).

---

## Chapter 4 — Running ORCHIDEE With Custom Meteorological Forcing

*(Priority Chapter)*

This chapter details how researchers can replace the default atmospheric forcing with custom meteorological datasets, including downscaled reanalyses, weather station observations, or Agrivoltaic (**AgriPV**) modified microclimates.

### 4.1 Use Case A: Comparing Alternative Forcing Datasets
When comparing multiple forcing products (e.g., Simulation A = Baseline ERA5, Simulation B = Custom AgriPV, Simulation C = Bias-Corrected WRF), **always maintain strictly isolated execution directories**.

```bash
# 1. Initialize three independent run directories
bash scripts/setup_offline_1day.sh ~/scratch/sim_control
bash scripts/setup_offline_1day.sh ~/scratch/sim_agripv_low
bash scripts/setup_offline_1day.sh ~/scratch/sim_agripv_high

# 2. In each directory, replace the forcing symlink with the target custom file:
cd ~/scratch/sim_agripv_low
rm forcing_19790201.nc
ln -s /path/to/forcing_agripv_low.nc forcing_19790201.nc
sbatch run_offline_1day.sh

cd ~/scratch/sim_agripv_high
rm forcing_19790201.nc
ln -s /path/to/forcing_agripv_high.nc forcing_19790201.nc
sbatch run_offline_1day.sh
```

### 4.2 Use Case B: Multi-Year & Multi-File Forcing Sequences
ORCHIDEE's standalone driver (`orchideedriver.f90`) natively supports reading multiple sequential NetCDF files covering multi-month or multi-year simulations.

In `run.def`, specify an array of filenames enclosed in parentheses:
```ini
# Multi-file annual sequence (up to 100 files supported natively)
FORCING_FILE = (forcing_1979.nc, forcing_1980.nc, forcing_1981.nc)

# Set simulation duration accordingly:
TIME_LENGTH = 3Y
```

#### Multi-file Rules Enforced by ORCHIDEE Driver:
1. **Time Axis Monotonicity:** The `time_counter` coordinate across `file1.nc`, `file2.nc`, etc. must be strictly monotonic and continuous.
2. **Matching Calendars:** All files must share the exact same calendar attribute (`gregorian`, `standard`, or `noleap`).
3. **No Overlapping Timesteps:** The end timestamp of file $N$ must immediately precede the start timestamp of file $N+1$ without duplicate records.

### 4.3 Technical NetCDF Forcing Specifications & Constraints

> [!WARNING]
> **NetCDF Compatibility Check:**  
> A custom meteorological NetCDF file is **NOT** plug-and-play with ORCHIDEE unless it adheres strictly to the CF-compliant variables, units, and spatial coordinates described below.

#### Required NetCDF Variables Dictionary

| Variable Name | Standard Name | Units | Dimension / Shape | Description & Typical Range |
| :--- | :--- | :--- | :--- | :--- |
| **`Tair`** | `air_temperature` | **K** | `(time_counter, lat, lon)` | 2m near-surface air temperature [220 K to 325 K] |
| **`Qair`** | `specific_humidity` | **kg/kg** | `(time_counter, lat, lon)` | 2m specific humidity [0.0001 to 0.030 kg/kg] |
| **`PSurf`** | `surface_air_pressure` | **Pa** | `(time_counter, lat, lon)` | Surface barometric pressure [50,000 to 105,000 Pa] |
| **`Wind`** | `wind_speed` | **m/s** | `(time_counter, lat, lon)` | 10m scalar wind speed [0.0 to 45.0 m/s] |
| **`Wind_E`** | `eastward_wind` | **m/s** | `(time_counter, lat, lon)` | Eastward (u) 10m wind vector component |
| **`Wind_N`** | `northward_wind` | **m/s** | `(time_counter, lat, lon)` | Northward (v) 10m wind vector component |
| **`Rainf`** | `rainfall_flux` | **kg m-2 s-1** | `(time_counter, lat, lon)` | Liquid precipitation rate (1 kg/m²/s = 1 mm/s) |
| **`Snowf`** | `snowfall_flux` | **kg m-2 s-1** | `(time_counter, lat, lon)` | Solid snowfall rate (set to 0.0 if no snow) |
| **`SWdown`** | `surface_downwelling_shortwave_flux_in_air` | **W m-2** | `(time_counter, lat, lon)` | Downwelling shortwave solar radiation [0 to 1200 W/m²] |
| **`LWdown`** | `surface_downwelling_longwave_flux_in_air` | **W m-2** | `(time_counter, lat, lon)` | Downwelling longwave infrared radiation [150 to 500 W/m²] |
| **`contfrac`** | `land_area_fraction` | **-** (0 to 1) | `(lat, lon)` or 3D | Continental land fraction (1 = 100% land, 0 = ocean) |

#### Coordinates and Dimensions
- **1D Coordinates:** `lon(lon)` and `lat(lat)`. On the regional Morocco zoom, dimensions must be `lon = 384`, `lat = 221`.
- **2D Auxiliary Coordinates:** `nav_lon(lat, lon)` and `nav_lat(lat, lon)` are **mandatory**. They allow IOIPSL and XIOS to locate curvilinear grid points.
- **Time Dimension:** `time_counter` (unlimited or fixed dimension), with attribute `units = "seconds since YYYY-MM-DD 00:00:00"`.
- **Temporal Frequency:** Supported native frequencies are 30 minutes, 1 hour, or 3 hours. ORCHIDEE automatically performs linear interpolation down to the model physics timestep (`TIME_STEP = 1800` s).

#### Established Constraints vs. Items to Confirm

| Parameter | Established Code Constraint | User Item to Confirm |
| :--- | :--- | :--- |
| **Precipitation Units** | **Must be $	ext{kg}\cdot	ext{m}^{-2}\cdot	ext{s}^{-1}$** | Convert if dataset is in mm/day or mm/hour ($\div 86400$ or $\div 3600$). |
| **Precipitation Partition** | `Rainf` and `Snowf` are separate fields | If dataset has total `precip`, compute `Rainf = max(0, precip - snow)` and `Snowf = snow`. |
| **Temperature Units** | **Must be Kelvin (K)** | Convert if dataset is in Celsius ($T_K = T_C + 273.15$). |
| **Radiation Units** | **Must be $	ext{W}\cdot	ext{m}^{-2}$** | Convert if dataset stores accumulated Joules ($J/	ext{m}^2 \div \Delta t$). |
| **Missing Values** | `_FillValue = 1.0e+20` or NaN on ocean | **Zero NaNs allowed on active land points.** |
| **Measurement Heights** | Specified in `run.def` | `HEIGHT_LEV1 = 2.0` (m) for T/Q, `HEIGHT_LEVW = 10.0` (m) for Wind. |

### 4.4 Adapting Meteorological Forcing for Agrivoltaic (AgriPV) Studies
Agrivoltaic systems alter the boundary layer microclimate through two primary mechanisms:
1. **Precipitation Interception & Concentration (`Rainf`):** Photovoltaic panels act as impermeable barriers. Part of the rainfall is intercepted and redistributed as concentrated drip lines along the bottom edges of the panel arrays:
   $$P_{	ext{under}} = (1 - f_{	ext{panel}} \cdot \eta_{	ext{capture}}) P_{	ext{ambient}}$$
   $$P_{	ext{dripline}} = P_{	ext{ambient}} + \left(rac{A_{	ext{panel}}}{A_{	ext{dripline}}}
ight) \cdot \eta_{	ext{capture}} \cdot P_{	ext{ambient}}$$
   When creating an AgriPV forcing file, modify the `Rainf` field accordingly to assess whether drip-line infiltration enhances deep recharge or causes surface runoff.
2. **Solar Radiation Shading (`SWdown`):** Fixed or tracking solar panels reduce the downwelling shortwave radiation reaching the underlying crop canopy:
   $$SW_{	ext{under}} = (1 - 	ext{ShadingFraction}(t)) \cdot SW_{	ext{ambient}}$$
   Reducing `SWdown` in the forcing NetCDF allows ORCHIDEE to calculate the resulting reduction in crop transpiration and soil evaporation.

---

## Chapter 5 — Configuring and Submitting an Experiment

### 5.1 Decoupling Configuration for Sensitivity Runs
By default, `setup_offline_1day.sh` creates a symlink to `configs/offline_1day/run.def`. To modify parameters for your experiment without affecting colleagues:
```bash
cd "${MY_RUN_DIR}"
rm run.def
cp ~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/configs/offline_1day/run.def ./run.def
```

### 5.2 Key Parameters in `run.def`

```ini
# Duration: e.g., 1D (1 day), 1M (1 month), 1Y (1 year), or 30D
TIME_LENGTH = 1D

# Physics time step in seconds (1800s = 30min)
TIME_STEP = 1800

# Name or array of atmospheric forcing files
FORCING_FILE = my_custom_forcing.nc

# Observation heights for temperature/humidity (2m) and wind (10m)
HEIGHT_LEV1 = 2.0
HEIGHT_LEVW = 10.0

# Initial Restart Configuration:
# Warm start from certified restart:
RESTART_FILEIN = NONE
SECHIBA_restart_in = sechiba_rest_in.nc
SECHIBA_rest_out = sechiba_rest_out.nc

# Cold start (spinup initialization without restart):
# SECHIBA_restart_in = NONE
```

### 5.3 Critical XML Rules for Stability on Toubkal
1. **Disable XIOS Remapping (`XIOS_REMAP_OUTPUT = n`):** On zoomed non-regular grids ($384 	imes 221$), setting this to `y` causes an MPI crash trying to interpolate to a global $144 	imes 143$ grid. Keep it set to `n`.
2. **Disable River Routing (`RIVER_ROUTING = n`):** Without precomputed coastal routing topology maps, routing calculates out-of-bound indices. `n` disables fluvial transport while preserving 100% of soil column infiltration, drainage, and runoff physics.
3. **Set Compression Level to 0 (`compression_level = "0"` in `file_def_orchidee.xml`):** NetCDF-4 parallel collective writing rejects deflate compression.

### 5.4 Slurm Batch Script Settings (`run_offline_1day.sh`)
```bash
#SBATCH --job-name=orch_ol_custom
#SBATCH --account=CLIMAT-UM6P-ST-IWRI-7KSIFKVWKUY-DEFAULT-CPU
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:30:00
#SBATCH --hint=nomultithread
```
Submit with:
```bash
sbatch run_offline_1day.sh
```

To rerun in the same directory, use the `--force` flag to clear previous outputs:
```bash
sbatch run_offline_1day.sh --force
```

---

## Chapter 6 — Model Outputs and Basic Diagnostics

### 6.1 Output Files Generated
1. **`sechiba_history.nc`:** Daily aggregated physical diagnostics (1 snapshot/day).
2. **`sechiba_out_2.nc`:** High-frequency diagnostics (3-hourly snapshots, 8 steps/day).
3. **`sechiba_rest_out.nc`:** Full model state restart at $t = 	ext{end}$, used as `sechiba_rest_in.nc` for subsequent runs.
4. **`driver_rest_out.nc`:** Atmospheric driver state restart.
5. **`used_run.def`:** Exhaustive record of all parameters used during the simulation.

### 6.2 Key Physical Variables Dictionary

| Output Variable | Source File | Standard Units | Physical Meaning | Typical Range (Morocco Zoom) |
| :--- | :--- | :--- | :--- | :--- |
| **`evap`** | `sechiba_history.nc` | **mm/day** | Total evapotranspiration rate | $-1.5$ to $+8.0	ext{ mm/d}$ |
| **`temp_sol`** | `sechiba_history.nc` | **°C** | Average surface skin temperature | $-50.0$ to $+50.0	ext{ }^\circ	ext{C}$ |
| **`rain`** | `sechiba_history.nc` | **mm/day** | Prescribed liquid precipitation received | $0.0$ to $200.0	ext{ mm/d}$ |
| **`netrad`** | `sechiba_history.nc` | **W m-2** | Net surface radiation absorbed | $-100.0$ to $+600.0	ext{ W/m}^2$ |
| **`mrsos`** | `sechiba_out_2.nc` | **kg m-2** | Topsoil moisture (0 to 10 cm depth) | $5.0$ to $45.0	ext{ kg/m}^2$ |
| **`mrro`** | `sechiba_out_2.nc` | **kg m-2 s-1** | Surface runoff rate | $0.0$ to $0.005	ext{ kg/m}^2/	ext{s}$ |
| **`TWS`** | `sechiba_out_2.nc` | **kg m-2** | Total terrestrial water storage in column | $300.0$ to $900.0	ext{ kg/m}^2$ |
| **`ptn`** | `sechiba_out_2.nc` | **K** | Multi-layer deep ground temperature (18 levels) | $260.0$ to $320.0	ext{ K}$ |
| **`lai`** | `sechiba_history.nc` | **m2/m2** | Leaf Area Index per Plant Functional Type | $0.0$ to $6.0	ext{ m}^2/	ext{m}^2$ |

### 6.3 Command-Line Inspection
Quickly inspect variable dimensions and metadata:
```bash
ncdump -h sechiba_history.nc | grep -E "temp_sol|evap|mrsow|rain"
```

### 6.4 Python Validation and Plotting Script
Researchers can use the following complete Python script to verify physical integrity and absence of NaNs:

```python
import netCDF4 as nc
import numpy as np

ds = nc.Dataset("sechiba_history.nc", "r")
evap = ds.variables["evap"][:]
temp_sol = ds.variables["temp_sol"][:]
contfrac = ds.variables["Contfrac"][:]

# Mask ocean points
land_mask = contfrac > 0.1
evap_land = evap[0, land_mask]
temp_land = temp_sol[0, land_mask]

print(f"Total active land points: {np.sum(land_mask)}")
print(f"NaN count (evap): {np.sum(np.isnan(evap_land))}")
print(f"Evaporation: Mean = {np.mean(evap_land):.2f} mm/d, Min = {np.min(evap_land):.2f}, Max = {np.max(evap_land):.2f}")
print(f"Surface Temp: Mean = {np.mean(temp_land):.2f} °C, Min = {np.min(temp_land):.2f}, Max = {np.max(temp_land):.2f}")

ds.close()
```

---

## Chapter 7 — Comprehensive Troubleshooting Guide

| Issue / Error Symptom | Root Cause | Diagnosis Command | Verified Resolution |
| :--- | :--- | :--- | :--- |
| **`Lmod: unknown module GCC/7.3.0`** | Toubkal hides legacy modules behind dot-prefixes. | `module avail GCC` | Source `modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env` which exposes `modules/foss-2018b`. |
| **`FORCING_FILE: No such file`** | The file specified in `run.def` is missing or symlink is broken. | `ls -l $(readlink -f FORCING_FILE)` | Ensure the forcing file exists or update `FORCING_FILE` in `run.def`. |
| **`iim_loc != iim dimension mismatch`** | NetCDF spatial grid does not match model domain ($384 	imes 221$). | `ncdump -h forcing.nc \| head -n 10` | Regrid your forcing to exactly $384 	imes 221$ coordinates matching the Morocco zoom. |
| **`Missing variable: precip or Rainf`** | Forcing NetCDF uses non-standard naming. | `ncdump -h forcing.nc` | Rename variables: liquid rain to `Rainf`, snow to `Snowf`, air temperature to `Tair`. |
| **`Negative values in Rainf`** | Calculation of `precip - snow` produced negative residuals. | `cdo fldmin -select,name=Rainf forcing.nc` | Ensure $Rainf = \max(0, precip - snow)$ during forcing preprocessing. |
| **`NetCDF: Impossible to use compression`** | XIOS collective parallel E/S rejects NetCDF-4 deflate compression. | Check `compression_level` in `file_def_orchidee.xml` | Set `compression_level="0"` in `file_def_orchidee.xml`. |
| **`domain_landpoints_out overlap crash`** | ORCHIDEE tries to re-interpolate regional zoom to global 144x143 grid. | Check `used_run.def` | Set `XIOS_REMAP_OUTPUT = n` in `run.def`. |
| **`mask_coast_glo(0) index error`** | River routing lacks precomputed coastal basins for zoom domain. | Check `out_orchidee_0000` | Set `RIVER_ROUTING = n` and `RIVER_DESC = n` in `run.def`. |
| **`Refusing to overwrite previous results`** | Safety check in `run_offline_1day.sh` prevents corrupting outputs. | `ls -l sechiba_rest_out.nc` | Submit with `--force` flag: `sbatch run_offline_1day.sh --force`. |
| **`OPAL ERROR in pmix2x_client.c`** | Direct `srun` fails with OpenMPI 3.1.1 on Toubkal. | Review batch script command | Always use `mpirun -np ${SLURM_NTASKS:-8} ./orchidee_ol`. |

---

## Chapter 8 — Quick Start Checklist & Pre-Flight Reference

### 8.1 Pre-Flight Verification Checklist
Before submitting any custom offline experiment, verify the following 5 checkpoints:

- [ ] **1. HPC Environment:** `source modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env` executed without error.
- [ ] **2. Forcing NetCDF Variables:** `Tair` (K), `Qair` (kg/kg), `PSurf` (Pa), `Wind` (m/s), `Rainf` (kg/m²/s), `Snowf` (kg/m²/s), `SWdown` (W/m²), `LWdown` (W/m²), and `contfrac` are present with 0 NaNs on land.
- [ ] **3. Spatial Grid:** Dimensions are exactly `lon = 384`, `lat = 221` with `nav_lon` and `nav_lat` 2D coordinates.
- [ ] **4. Configuration Safety:** `run.def` decoupled in your personal directory (`XIOS_REMAP_OUTPUT = n`, `RIVER_ROUTING = n`, `compression_level = "0"`).
- [ ] **5. Workspace Isolation:** Job launched from a personal directory (never from `simulations/offline_1day`).

### 8.2 Parameter Quick Reference Table

| Target Configuration | Parameter in `run.def` | Recommended Setting | Purpose |
| :--- | :--- | :--- | :--- |
| **Standard 1-Day Run** | `TIME_LENGTH` | `1D` | Simulates 24 hours (86,400 seconds). |
| **Monthly Run** | `TIME_LENGTH` | `1M` | Simulates 1 calendar month. |
| **Multi-Year Run** | `TIME_LENGTH` | `3Y` (or `365D`) | Extended multi-annual simulation. |
| **Physics Timestep** | `TIME_STEP` | `1800` | 30-minute numerical integration. |
| **Forcing File** | `FORCING_FILE` | `my_forcing.nc` | Direct path or array of forcing files. |
| **Initial Condition** | `SECHIBA_restart_in` | `sechiba_rest_in.nc` | Point to valid restart, or `NONE` for cold start. |
| **Diagnostic Level** | `output_level` in XML | `5` | Exports full water and energy balance fields. |
| **MPI Parallel Tasks** | `#SBATCH --ntasks` | `8` | 8 MPI ranks provides optimal performance on 1 node. |

### 8.3 3-Command Execution Cheat Sheet

```bash
# 1. Setup execution directory:
bash ~/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/scripts/setup_offline_1day.sh ~/scratch/my_new_run

# 2. Link custom forcing:
cd ~/scratch/my_new_run && rm forcing_19790201.nc && ln -s /path/to/my_forcing.nc forcing_19790201.nc

# 3. Launch Slurm job:
sbatch run_offline_1day.sh
```

---
*ORCHIDEE Offline on Toubkal HPC User Guide — v0.1 (November 2026)*  
