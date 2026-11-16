# ORCHIDEE Offline on Toubkal HPC — User Guide for Customizing Simulations

**Author:** International Water Research Institute (IWRI), Mohammed VI Polytechnic University (UM6P)  
**Target Audience:** Researchers, climate modelers, agronomists, and graduate students using standalone ORCHIDEE Offline on Toubkal HPC  
**Document Path:** `docs/userguide/customize_simulation_guide.md`  
**Prerequisite:** Verified reproduction of the reference benchmark (`docs/userguide/reproduction_guide.md`)  
**Toolchain & Model:** GCC 12.3.0 / OpenMPI 4.1.5 (`gompi/2023a`), ORCHIDEE Quest (SVN r7086), XIOS 2 (SVN r2466), IOIPSL v2_2_5  
**Reference Case:** 1-Day Standalone Simulation (February 1, 1979) on 35 km North Africa Zoom Grid  

---

## Executive Overview

You have successfully completed the 7-step reproduction guide and verified that ORCHIDEE Offline runs and reproduces certified numerical outputs on Toubkal HPC. You now have a working executable (`bin/orchidee_ol_gcc12`) and a validated benchmark directory (`simulations/offline_1day_gcc12/`).

**This guide is the second step in your modeling journey.** Its purpose is to teach you how to transition from simply *reproducing* an existing simulation to *customizing* your own scientific experiments:
- Running **climate sensitivity experiments** (e.g., assessing the impact of a +20% increase in rainfall, drought stress, heatwaves, or solar radiation changes).
- Extending the **simulation duration** from a single day to multi-day or seasonal episodes.
- Understanding the architectural prerequisites for **changing the spatial domain** (e.g., regional catchments or Agri-Photovoltaics [AgriPV] study sites).
- Modifying **soil physics, surface albedo, and vegetation parameters** in the model namelist.
- Controlling **diagnostic output frequencies and variables** via XIOS XML definitions.

> [!IMPORTANT]
> **Cardinal Golden Rule of Experimentation on Toubkal HPC:**  
> **Never modify the shared reference installation or configuration files.**  
> Always conduct your experiments in isolated, dedicated directories under your personal Lustre workspace (`/srv/lustre01/project/.../users/$USER/`). The reference configuration (`configs/offline_1day/`) must remain strictly pristine.

---

## Part I — Understanding an ORCHIDEE Offline Simulation

A standalone (offline) land surface simulation does not run in a vacuum. It is a coupled ecosystem of software libraries, configuration namelists, XML descriptors, meteorological boundary conditions, and initial land surface states.

```
+-----------------------------------------------------------------------------------+
|                            ORCHIDEE Standalone Driver                             |
|                                (orchidee_ol_gcc12)                                |
+-----------------------------------------+-----------------------------------------+
                                          |
        +---------------------------------+---------------------------------+
        |                                                                   |
        v                                                                   v
+-------------------------------+                         +-------------------------------+
|            IOIPSL             |                         |            XIOS 2             |
|   (Input & Time Management)   |                         |     (Output Diagnostics)      |
+---------------+---------------+                         +---------------+---------------+
                |                                                         |
  +-------------+-------------+                             +-------------+-------------+
  |                           |                             |                           |
  v                           v                             v                           v
run.def            forcing_19790201.nc               iodef.xml           file_def_orchidee.xml
(Namelist switches  (Atmospheric 3-hourly            (Server parameters   (Output file names,
 & physical params)  T, Q, Wind, Rain, SW, LW)        & context links)     levels & frequencies)
```

### 1. The Core Components
1. **ORCHIDEE (`modeles/ORCHIDEE`):**  
   The biophysical and biogeochemical engine. In our offline setup (`OOL_SEC`), it calculates the 11-layer soil hydrology (CWRR diffusion scheme), energy balance, latent and sensible heat fluxes, surface runoff, deep drainage, soil temperature profiles, and snow dynamics.
2. **IOIPSL (`modeles/IOIPSL`):**  
   The Input/Output library developed by IPSL. It reads the model namelist (`run.def`), handles atmospheric forcing reading and calendar time-stepping, and manages spatial interpolation routines.
3. **XIOS (`modeles/XIOS`):**  
   The XML-IO-Server client/server library. It handles all output diagnostics asynchronously. Instead of hardcoding output fields in Fortran, XIOS dynamically controls which variables are written, at what temporal frequency, and into which NetCDF files based on XML configuration files.

### 2. The Key Control Files
- **`run.def`:**  
  The primary Fortran namelist. It specifies runtime parameters such as simulation duration (`TIME_LENGTH`), time step (`TIME_STEP`), the atmospheric forcing file name (`FORCING_FILE`), hydrology options (`HYDROL_CWRR`), and initial restart file names (`SECHIBA_restart_in`).
- **XIOS XML Files:**
  - `iodef.xml`: Top-level entry point defining the XIOS server context (client mode vs server mode).
  - `context_orchidee.xml`: Defines spatial grids and operational axes for ORCHIDEE.
  - `file_def_orchidee.xml`: Defines which output NetCDF files are created (e.g. `sechiba_history.nc`, `sechiba_out_2.nc`), their output frequency (`output_freq`), and diagnostic depth (`output_level`).
  - `field_def_orchidee.xml`: Master catalog of all available physical variables and their scientific units.
- **Atmospheric Forcing NetCDF (`forcing_19790201.nc`):**  
  The meteorological boundary condition at 2 m and 10 m heights driving the land surface. It contains 3-hourly time series of near-surface air temperature (`Tair`), specific humidity (`Qair`), surface pressure (`PSurf`), wind speed (`Wind`), precipitation fluxes (`Rainf`, `Snowf`), and downwelling shortwave/longwave radiation (`SWdown`, `LWdown`).
- **Initial State / Restart File (`sechiba_rest_in.nc`):**  
  Stores the 3D state of the land surface at the beginning of the run (soil moisture in 11 layers, soil temperatures, canopy water interception, snow pack).
- **Slurm Execution Script (`run_offline_1day.sh`):**  
  Allocates compute resources (8 MPI tasks on 1 compute node), loads the `gompi/2023a` environment, sets stack limits (`ulimit -s unlimited`), and launches the model via `mpirun`.

### 3. Reference Configuration vs. Execution Directory
To ensure team reproducibility, this repository strictly separates *code/config templates* from *active execution directories*:
- **Reference Configuration (`configs/offline_1day/`):**  
  Version-controlled templates (`run.def`, XML files). Never launch simulations directly inside this directory.
- **Execution Directory (e.g., `simulations/offline_1day_gcc12/`):**  
  The working directory where the model binary is executed, containing symbolic links to the binary, configs, and input data. Diagnostic outputs (`*.nc`) and execution logs (`slurm-*.out`) are generated here.

### 4. Customization Reference Table

| What do you want to change? | Which file to modify? | What must be checked? |
| :--- | :--- | :--- |
| **Simulation duration** | `run.def` (`TIME_LENGTH`) | Forcing temporal coverage, start date, calendar consistency |
| **Atmospheric forcing** | `run.def` (`FORCING_FILE`), NetCDF | Variable names, units (`kg/m^2/s`, `K`), time axis, 2D coordinates |
| **Precipitation flux** | Custom NetCDF (`Rainf`, `Snowf`) | Units (`kg/m^2/s`), strictly non-negative values, mask consistency |
| **Air temperature / Radiation** | Custom NetCDF (`Tair`, `SWdown`) | Physical realism ($T > 0\text{ K}$, $SW \ge 0\text{ W/m}^2$), units |
| **Output diagnostic frequency** | `file_def_orchidee.xml` | `output_freq` attribute (`1h`, `3h`, `1d`, `1mo`), file sizes |
| **Output variable selection** | `file_def_orchidee.xml` | Field references must exist in `field_def_orchidee.xml`, `output_level` |
| **Soil hydrology physics** | `run.def` (`HYDROL_CWRR`) | `y` for 11-layer CWRR diffusion, `n` for 2-layer bucket scheme |
| **Vegetation parameters** | `run.def` (`VCMAX25__*`, `LAI_MAX__*`) | Target PFT index (1 to 15), physical parameter bounds |
| **Spatial domain / Grid** | Forcing, Restart, PFT, Soil NetCDF | Resolution ($N_{lon} \times N_{lat}$), `nav_lon`/`nav_lat`, `contfrac` mask |

---

## Part II — Creating an Independent Experiment

To run any new simulation or sensitivity test, you must create a dedicated execution directory.

### Step-by-Step Setup Workflow

```bash
# 1. Navigate to your modipsl repository root
cd /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/users/$USER/.../modipsl

# 2. Initialize a new experiment directory using the automated setup script
bash scripts/setup_offline_1day_gcc12.sh simulations/my_new_experiment

# 3. Enter your new experiment workspace
cd simulations/my_new_experiment

# 4. Inspect the directory contents
ls -l
```

### How the Setup Script Works
`setup_offline_1day_gcc12.sh` automatically creates 17 lightweight symbolic links:
- `orchidee_ol -> .../modipsl/bin/orchidee_ol_gcc12`
- `run.def -> .../modipsl/configs/offline_1day/run.def`
- XML files: `iodef.xml`, `context_orchidee.xml`, `file_def_orchidee.xml`, etc.
- Boundary conditions: `soils_param.nc`, `PFTmap.nc`, `alb_bg.nc`, `routing.nc`, etc.
- Restart & Forcing: `sechiba_rest_in.nc`, `forcing_19790201.nc`
- Slurm runner: `run_offline_1day.sh -> .../scripts/run_offline_1day_gcc12.sh`

### Decoupling Symbolic Links for Customization
Because the links point to read-only or shared files, you customize an experiment by **breaking specific links and replacing them with local, editable files**:

```bash
# Example: Customizing the namelist run.def
# 1. Remove the symbolic link
rm run.def

# 2. Copy the reference file to create a local editable copy
cp ../../configs/offline_1day/run.def ./run.def

# 3. Edit your local run.def safely
nano run.def   # or vi, emacs, etc.
```

> [!TIP]
> `scripts/setup_offline_1day_gcc12.sh` features a built-in conflict detector. If a local file already exists (such as your modified `run.def`), running the setup script again will **preserve your personal file** without overwriting it.

---

## Part III — Exercise 1: Changing an Atmospheric Forcing Variable (Precipitation Sensitivity: +20% Rain)

### 1. Objective
Simulate a climate sensitivity experiment where rainfall is increased by **+20%** across the entire domain on February 1, 1979. We will compare the hydrological response (evapotranspiration, runoff, drainage, soil moisture) against the reference simulation.

> **Status:** Fully Validated Methodology. Tested on Toubkal HPC.

### 2. Files to Inspect and Modify
- **Source Forcing File:** `forcing_19790201.nc` (38 MB NetCDF4 file, 8 time steps of 3 hours).
- **Target Modified File:** `forcing_rain_plus20.nc`.
- **Target Variable:** `Rainf` (Rainfall flux, units: `kg/m^2/s`).

### 3. Step-by-Step Instructions

#### Step 1: Create an Independent Experiment Workspace
```bash
cd modipsl
bash scripts/setup_offline_1day_gcc12.sh simulations/exp_precip_plus20
cd simulations/exp_precip_plus20
```

#### Step 2: Generate the Modified Forcing File in Python
ORCHIDEE requires that all metadata, dimensions, coordinate arrays (`lon`, `lat`, `nav_lon`, `nav_lat`), and auxiliary variables remain perfectly intact. We use a Python script with `netCDF4` to safely scale `Rainf`:

```bash
cat << 'PYEOF' > modify_rain.py
import netCDF4 as nc
import numpy as np

src_path = "forcing_19790201.nc"
dst_path = "forcing_rain_plus20.nc"

print(f"Reading reference forcing from: {src_path}")
with nc.Dataset(src_path, "r") as src, nc.Dataset(dst_path, "w", format="NETCDF4") as dst:
    # 1. Copy global attributes
    dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
    dst.title = "ORCHIDEE Offline Forcing - Precipitation Sensitivity (+20% Rainf)"

    # 2. Copy dimensions
    for name, dim in src.dimensions.items():
        dst.createDimension(name, (len(dim) if not dim.isunlimited() else None))

    # 3. Copy all variables, scaling Rainf by 1.20
    for name, var in src.variables.items():
        fill_val = getattr(var, "_FillValue", None)
        out_var = dst.createVariable(
            name, var.dtype, var.dimensions,
            zlib=True, complevel=4, fill_value=fill_val
        )
        out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs() if k != "_FillValue"})

        data = var[:]
        if name == "Rainf":
            print(f"Applying +20% scaling to '{name}'...")
            out_var[:] = (data * 1.20).astype(var.dtype)
        else:
            out_var[:] = data

print(f"Successfully generated: {dst_path}")
PYEOF
python3 modify_rain.py
```

#### Step 3: Verify the Modified Forcing File
```bash
# Check file size (~38 MB)
ls -lh forcing_rain_plus20.nc

# Verify that Rainf increased by exactly 20% while Snowf and Tair remained unchanged
python3 -c "
import netCDF4 as nc
import numpy as np

ref = nc.Dataset('forcing_19790201.nc')
mod = nc.Dataset('forcing_rain_plus20.nc')

r_ref = ref.variables['Rainf'][:]
r_mod = mod.variables['Rainf'][:]
diff_ratio = np.mean(r_mod[r_ref > 0]) / np.mean(r_ref[r_ref > 0])

print(f'Mean Rainf Reference : {np.mean(r_ref):.6e} kg/m^2/s')
print(f'Mean Rainf Modified  : {np.mean(r_mod):.6e} kg/m^2/s')
print(f'Scaling ratio        : {diff_ratio:.4f} (Expected: 1.2000)')

s_ref = ref.variables['Snowf'][:]
s_mod = mod.variables['Snowf'][:]
print(f'Snowf Max Abs Diff   : {np.max(np.abs(s_mod - s_ref))}')
ref.close()
mod.close()
"
```

#### Step 4: Link the Modified Forcing into the Experiment
Replace the original `forcing_19790201.nc` symlink with your new forcing file:
```bash
rm forcing_19790201.nc
ln -s forcing_rain_plus20.nc forcing_19790201.nc
```

#### Step 5: Submit the Simulation via Slurm
```bash
sbatch run_offline_1day.sh
```

#### Step 6: Monitor and Verify Normal Completion
```bash
# Check queue
squeue -u $USER

# Once finished (~25-30 seconds), inspect the end of the Slurm log
tail -n 10 slurm-*.out
```
*Expected log output:*
```text
==========================================================
Fin normale de la simulation : ... 2026
==========================================================
```

#### Step 7: Compare the Scientific Results
Compare the outputs of your new experiment against the reference benchmark (`../offline_1day_gcc12`):

```bash
python3 - << 'PYEOF'
import netCDF4 as nc
import numpy as np

ref_path = "../offline_1day_gcc12/sechiba_history.nc"
exp_path = "./sechiba_history.nc"

with nc.Dataset(ref_path) as ds_ref, nc.Dataset(exp_path) as ds_exp:
    print("="*75)
    print("HYDROLOGICAL IMPACT OF +20% PRECIPITATION (sechiba_history.nc)")
    print("="*75)
    
    vars_to_check = [
        ("evap", "Evapotranspiration", "mm/d"),
        ("runoff", "Surface Runoff", "mm/d"),
        ("drainage", "Deep Soil Drainage", "mm/d"),
        ("mrso", "Total Soil Moisture", "kg/m^2"),
        ("mrsos", "Top-soil Moisture (10cm)", "kg/m^2")
    ]
    
    for vname, desc, unit in vars_to_check:
        v_ref = ds_ref.variables[vname][0, :, :]
        v_exp = ds_exp.variables[vname][0, :, :]
        
        mask = np.isnan(v_ref) | (v_ref == 0.0)
        mean_ref = np.mean(v_ref[~mask])
        mean_exp = np.mean(v_exp[~mask])
        delta = mean_exp - mean_ref
        pct = (delta / mean_ref) * 100.0 if mean_ref != 0 else 0.0
        
        print(f"{desc:<26} ({vname}):")
        print(f"  Reference : {mean_ref:10.4f} {unit}")
        print(f"  +20% Rain : {mean_exp:10.4f} {unit}")
        print(f"  Change    : {delta:+10.4f} {unit} ({pct:+.2f}%)")
        print("-" * 75)
PYEOF
```

### 4. Expected Results & Scientific Interpretation
- **Rainfall Increase:** Precipitation is directly partitioned between soil infiltration and direct surface runoff.
- **Runoff (`runoff`):** Shows immediate, sensitive response over saturated or low-permeability soils (expected increase > 20% in mountainous or high-slope regions).
- **Drainage (`drainage`):** Gradual response as water percolates down the 11 CWRR layers.
- **Evapotranspiration (`evap`):** Increases in moisture-limited regions (semi-arid areas of North Africa) where evaporation is water-limited rather than energy-limited.

---

## Part IV — Exercise 2: Extending the Simulation Period (From 1 Day to Multi-Day)

> [!WARNING]
> **Status: PROPOSED ARCHITECTURAL DESIGN PATTERN — NOT YET TESTED ON TOUBKAL**  
> Chaining multi-day runs requires custom multi-day forcing extraction and restart management. This section explains the exact technical mechanism for future experiments.

### 1. Objective
Understand how ORCHIDEE handles time, how to extend a simulation from 1 day (48 time steps) to 3 days (144 time steps), and what constraints govern the time axis.

### 2. Critical Caveat: Why Changing `TIME_LENGTH` Alone is Not Enough
A common beginner mistake is simply changing `TIME_LENGTH = 1D` to `TIME_LENGTH = 3D` in `run.def`.  
**This will immediately crash the simulation with an end-of-file (EOF) error from IOIPSL.**  
ORCHIDEE reads atmospheric forcing at every time step. If your forcing file only contains 8 time steps (24 hours at 3-hour resolution), the driver will abort when attempting to read step 9.

### 3. The Prerequisites for Multi-Day Simulations
To extend a simulation, four components must be aligned:
1. **The Forcing File:** Must contain data for the entire target duration. For 3 days at 3-hourly frequency, the file must contain $3 \times 8 = 24$ records.
2. **The Namelist Parameter (`TIME_LENGTH`):**  
   In `run.def`, change:
   ```ini
   TIME_LENGTH = 3D
   ```
3. **The Time Step (`TIME_STEP`):**  
   Default is `1800` seconds (30 minutes). 1 day = 48 steps. 3 days = 144 steps.
4. **XIOS Diagnostics Frequency:**  
   In `file_def_orchidee.xml`, daily diagnostics (`output_freq="1d"`) will now record 3 distinct temporal slices in `sechiba_history.nc`.

### 4. The Time Axis and `DRIVER_reset_time`
In standalone ORCHIDEE, the parameter `DRIVER_reset_time` controls how the model calendar synchronizes with the forcing:
```ini
DRIVER_reset_time = y
```
- When set to `y`, the offline driver resets its internal calendar clock to `1979-01-01 00:00:00` for step 1, regardless of the calendar date stored in the forcing file.
- The output files reflect this origin (`time_centered:time_origin = "1979-01-01 00:00:00"`).
- When chaining successive simulations (e.g. Day 1 restart into Day 2), `SECHIBA_restart_in` must point to the restart produced at the end of the previous run (`sechiba_rest_out.nc`).

### 5. Step-by-Step Procedure for a 3-Day Experiment (Proposed)

#### Step 1: Extract 3 Days of Forcing from the Master Archive
The master historical archive file on Toubkal (`sim35km-HIST_19790201_19790228_HF_histhf.nc`) contains 224 time steps (the full month of February 1979). Using `scripts/prepare_forcing_1day.py`:

```bash
python3 scripts/prepare_forcing_1day.py \
  --input /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist/ATM/Output/HF/sim35km-HIST_19790201_19790228_HF_histhf.nc \
  --output simulations/exp_3days/forcing_19790201_3days.nc \
  --start 0 \
  --steps 24 \
  --force
```

#### Step 2: Configure the 3-Day Run Directory
```bash
# 1. Initialize directory
bash scripts/setup_offline_1day_gcc12.sh simulations/exp_3days
cd simulations/exp_3days

# 2. Decouple run.def
rm run.def
cp ../../configs/offline_1day/run.def ./run.def

# 3. Update run.def parameters
sed -i "s|^TIME_LENGTH = .*|TIME_LENGTH = 3D|" run.def
sed -i "s|^FORCING_FILE = .*|FORCING_FILE = forcing_19790201_3days.nc|" run.def

# 4. Link the new forcing file
ln -sf forcing_19790201_3days.nc forcing_19790201_3days.nc
```

#### Step 3: Adjust Slurm Walltime and Launch
A 1-day simulation runs in ~25 seconds. A 3-day simulation takes approximately 75–90 seconds. The default Slurm script walltime (`#SBATCH --time=00:30:00`) is more than sufficient.
```bash
sbatch run_offline_1day.sh
```

#### Step 4: Verification of 3-Day Outputs
```bash
# Check that the time dimension has 3 steps in daily history
ncdump -h sechiba_history.nc | grep "time_counter = "
# Expected: time_counter = UNLIMITED ; // (3 currently)

# Check that the 3-hourly history has 24 steps
ncdump -h sechiba_out_2.nc | grep "time_counter = "
# Expected: time_counter = UNLIMITED ; // (24 currently)
```

---

## Part V — Exercise 3: Changing the Simulation Domain

> [!WARNING]
> **Status: CONCEPTUAL ARCHITECTURAL WORKFLOW — REQUIRES NEW ANCILLARY DATASETS**  
> Changing the spatial domain cannot be done by editing a configuration switch. It requires preparing a complete, mutually consistent family of spatial input files.

### 1. Objective
Understand the data architecture required to relocate ORCHIDEE Offline to a different geographic domain (e.g. an agricultural catchment, an AgriPV site in southern Morocco, or a Mediterranean sub-basin).

### 2. Why Domain Change Requires New Datasets
In ORCHIDEE, the computational grid is defined by the input NetCDF files. Unlike atmospheric models that can generate an arbitrary analytical grid at runtime, ORCHIDEE extracts its grid dimensions ($N_{lon} \times N_{lat}$) and coordinate arrays directly from the forcing and boundary files.

All spatial input files **must have the exact same dimensions and matching land/sea masks**:

```
+-----------------------------------------------------------------------------------+
|                        MANDATORY DOMAIN-CONSISTENT FAMILY                         |
+-----------------------------------------------------------------------------------+
| 1. Meteorological Forcing | forcing.nc        | Tair, Rainf, Wind, nav_lon, etc.  |
| 2. Land-Sea Mask          | contfrac          | Continental fraction (0 to 1)     |
| 3. Initial Restart        | sechiba_rest_in.nc| Soil moisture & temp (11 layers)  |
| 4. Vegetation Map         | PFTmap.nc         | 15 Plant Functional Types (PFTs)  |
| 5. Soil Hydraulics        | soils_param.nc    | Soil texture, conductivity        |
| 6. Surface Albedo         | alb_bg.nc         | Background albedo from MODIS      |
| 7. Reference Climatology  | reftemp.nc        | Climatological reference temp     |
| 8. Topographic Slopes     | cartepente2d.nc   | Surface slope distributions       |
+-----------------------------------------------------------------------------------+
```

### 3. Step-by-Step Data Preparation Pipeline

```
  [Target Domain Bounds]
 (lon_min, lon_max, lat_min, lat_max)
           |
           +---> 1. Atmospheric Forcing (ERA5-Land, CORDEX, or WRF)
           |        -> Interpolate to target grid
           |        -> Calculate nav_lon, nav_lat, Rainf, Snowf
           |
           +---> 2. Surface Parameters (IGCM Database)
           |        -> PFTmap.nc (Cut/Remap CMIP6 15-PFT map)
           |        -> soils_param.nc (USDA or Zobler soil textures)
           |        -> alb_bg.nc (MODIS background albedo)
           |        -> cartepente2d.nc (Topographic slopes)
           |
           +---> 3. Initial State / Spin-up
                    -> Option A: Spatial interpolation of existing restart
                    -> Option B: Cold start (uniform wet soil) + multi-year spinup
```

#### Detailed Steps:
1. **Define the Grid:** Choose your coordinate system (regular lat/lon or rotated grid). Generate 2D coordinate matrices `nav_lon(lat, lon)` and `nav_lat(lat, lon)`.
2. **Prepare Meteorological Forcing:** Extract atmospheric variables from high-resolution reanalyses (e.g. ERA5-Land at 9 km) or regional models. Ensure variable units match ORCHIDEE standards (`Rainf` in `kg/m^2/s`, `Tair` in `K`, `PSurf` in `Pa`).
3. **Clip/Interpolate IGCM Surface Files:** Use CDO (Climate Data Operators) or Python scripts to remap the global Toubkal IGCM boundary datasets (`/srv/lustre01/project/.../igcmg/IGCM/SRF/`) to your target grid.
4. **Generate the Restart (`sechiba_rest_in.nc`):**  
   - If migrating to a new domain where no coupled restart exists, you must initialize from a "cold start" (setting `RESTART_FILEIN = NONE` and `SECHIBA_restart_in = NONE` in `run.def`).
   - Run a multi-year spinup (3 to 10 years) using repeated meteorological forcing to allow root-zone soil moisture and deep drainage to reach equilibrium.
5. **Disable XIOS Output Remapping:**  
   Always ensure `XIOS_REMAP_OUTPUT = n` in `run.def` when using regional or zoom grids.

---

## Part VI — Other Useful Customizations

### 1. Modifying Model Physics Switches in `run.def`

All switches can be edited directly in your local `run.def`:

- **Soil Hydrology Scheme (`HYDROL_CWRR`):**
  ```ini
  # y = 11-layer multilayer diffusion scheme (recommended for realistic soil hydrology)
  # n = Legacy 2-layer bucket scheme
  HYDROL_CWRR = y
  ```
- **Soil Freezing (`OK_FREEZE`):**
  ```ini
  # Activates latent heat release/absorption during soil water freezing and thawing
  OK_FREEZE = y
  ```
- **Explicit Multi-Layer Snow Physics (`OK_EXPLICITSNOW`):**
  ```ini
  # Multi-layer snow pack scheme with variable snow density and albedo aging
  OK_EXPLICITSNOW = y
  ```
- **Surface Roughness Length (`ROUGH_DYN`):**
  ```ini
  # Calculates dynamic aerodynamic roughness length based on vegetation height
  ROUGH_DYN = y
  ```

### 2. Customizing Vegetation & Canopy Parameters
You can adjust photosynthetic capacity ($V_{cmax}$) and Maximum Leaf Area Index ($LAI_{max}$) for specific Plant Functional Types (PFTs). For example, for **PFT 10 (Temperate C3 Grass)**:

```ini
# Maximum carboxylation rate at 25°C (micromol CO2 / m^2 / s)
VCMAX25__10 = 50.0

# Maximum seasonal Leaf Area Index (m^2 leaf / m^2 ground)
LAI_MAX__10 = 2.5
```

### 3. Customizing Output Frequencies in `file_def_orchidee.xml`
Decouple `file_def_orchidee.xml` in your experiment directory (`rm file_def_orchidee.xml && cp ../../configs/offline_1day/file_def_orchidee.xml .`).

To change the output frequency of daily history (`sechiba_history`) to hourly or 3-hourly:
```xml
<!-- Default: daily (1d) -->
<file id="sechiba1" name="sechiba_history" output_level="5" output_freq="1d" enabled=".TRUE.">

<!-- Customized: Hourly output -->
<file id="sechiba1" name="sechiba_history" output_level="5" output_freq="3600s" enabled=".TRUE.">

<!-- Customized: Monthly output -->
<file id="sechiba1" name="sechiba_history" output_level="5" output_freq="1mo" enabled=".TRUE.">
```

### 4. Selecting Diagnostic Output Variables in XIOS
To add a specific variable to `sechiba_history.nc`:
1. Check that the variable exists in `field_def_orchidee.xml` (e.g. `field_ref="temp_sol"`).
2. Add the variable to the `<field_group>` inside `file_def_orchidee.xml`:
   ```xml
   <field field_ref="temp_sol" level="1" />
   <field field_ref="evap"     level="1" />
   <field field_ref="runoff"   level="1" />
   <field field_ref="drainage" level="1" />
   ```
3. Set `output_level="5"` to ensure level <= 5 fields are written.

---

## Part VII — Checking and Comparing Results

### 1. Monitoring Slurm Execution
```bash
# Monitor active jobs in queue
squeue -u $USER

# Inspect resource consumption and elapsed time after completion
sacct -j <JOBID> --format=JobID,JobName,Partition,AllocCPUS,Elapsed,State,ExitCode
```

### 2. Inspecting Log Files
Always verify that the simulation completed without silent errors:
```bash
# Check standard output for normal termination
grep "Fin normale de la simulation" slurm-*.out

# Check model internal driver summary
grep -E "CPU Time Global|Real Time Global|END of dim2_driver" out_orchidee_0000

# Check that standard error is completely empty
cat slurm-*.err
```

### 3. Inspecting NetCDF Output Files
```bash
# Check dimensions and variables of daily output
ncdump -h sechiba_history.nc | head -n 35

# Check simulation time coordinates
ncdump -v time_centered sechiba_history.nc
```

### 4. Python Quick-Look Analysis Script
Use this script to compute spatial statistics on your simulation outputs:

```bash
python3 - << 'PYEOF'
import netCDF4 as nc
import numpy as np

with nc.Dataset("sechiba_history.nc") as ds:
    print(f"Simulation Output: {ds.title if hasattr(ds, 'title') else 'ORCHIDEE Diagnostics'}")
    
    for vname in ["temp_sol", "evap", "runoff", "drainage", "mrso"]:
        if vname in ds.variables:
            var = ds.variables[vname]
            data = var[0, :, :]
            mask = np.isnan(data) | (data == 0.0)
            valid = data[~mask]
            
            unit = getattr(var, "units", "N/A")
            print(f"Variable: {vname:<10} | Min: {np.min(valid):10.2f} | Max: {np.max(valid):10.2f} | Mean: {np.mean(valid):10.2f} {unit}")
PYEOF
```

### 5. Final Experimenter Pre-Flight Checklist

Before submitting `sbatch run_offline_1day.sh`, verify each item:

- [ ] I am working in an independent directory under `simulations/` (not in `configs/` or `shared/`).
- [ ] My `run.def` is a local editable file if I modified namelist parameters.
- [ ] My atmospheric forcing file (`FORCING_FILE`) exists and covers the requested `TIME_LENGTH`.
- [ ] Precipitation and radiation fluxes in custom forcings have strictly non-negative values.
- [ ] `XIOS_REMAP_OUTPUT = n` is set in `run.def`.
- [ ] Old simulation outputs (`sechiba_*.nc`, `out_orchidee_*`) from previous runs in the same folder are removed (or `--force` flag is used).
- [ ] `orchidee_ol` points to a valid, executable binary.

---

## Part VIII — Troubleshooting and Next Steps

### 1. Common Pitfalls & Diagnostic Solutions

#### Symptom 1: Job crashes immediately with `ERROR: Previous output files exist`
- **Cause:** `scripts/run_offline_1day_gcc12.sh` includes safety protection against accidental overwrites.
- **Solution:** Either remove existing outputs (`rm -f *.nc out_orchidee_*`) or force execution:
  ```bash
  sbatch run_offline_1day.sh --force
  ```

#### Symptom 2: Crash with `NetCDF: Index exceeds dimension bound` or `EOF on forcing`
- **Cause:** Your simulation duration (`TIME_LENGTH`) is longer than the number of time steps available in your atmospheric forcing file.
- **Solution:** Check time step count in forcing: `ncdump -h forcing.nc | grep "time_counter = "`. Ensure $(N_{steps} \times \text{forcing\_interval}) \ge \text{TIME\_LENGTH}$.

#### Symptom 3: Crash with `segmentation fault` or `stack overflow`
- **Cause:** Linux shell stack size limit was not set to unlimited.
- **Solution:** Ensure `ulimit -s unlimited` is active in your submission script (this is already pre-configured in `run_offline_1day_gcc12.sh`).

#### Symptom 4: XIOS XML Parse Error (`XML error: attribute '...' not found`)
- **Cause:** Syntax error or unclosed tag in your customized `file_def_orchidee.xml`.
- **Solution:** Validate XML syntax using `xmllint --noout file_def_orchidee.xml`.

#### Symptom 5: All land points show NaN or physically impossible values
- **Cause:** Unit mismatch in forcing. For example, supplying rainfall in `mm/day` instead of `kg/m^2/s` ($1\text{ mm/day} \approx 1.1574 \times 10^{-5}\text{ kg/m}^2/\text{s}$).

---

### 2. Next Steps: Multi-Year Simulations and libIGCM

For multi-year, decadal, or coupled climate simulations, manually managing daily restarts and forcing files becomes cumbersome.

In the standard IPSL modeling environment, long-term simulations are automated using **libIGCM** (the IPSL model workflow and job manager). libIGCM automatically slices multi-year simulations into monthly or annual segments, archives outputs to mass storage, manages calendar rollovers, and generates automated diagnostic web pages.

> [!NOTE]
> **libIGCM Status on Toubkal HPC:**  
> **libIGCM is NOT yet validated for this standalone GCC 12 Offline configuration on Toubkal.**  
> Our current offline setup uses lightweight, standalone bash drivers (`scripts/setup_offline_1day_gcc12.sh`, `scripts/run_offline_1day_gcc12.sh`), providing maximum transparency and direct control for short-term and sensitivity experiments. Future work by the IWRI climate team will address libIGCM adaptation to Toubkal's modern software stack.

---

## Summary of Useful Commands

```bash
# 1. Initialize any new experiment
bash scripts/setup_offline_1day_gcc12.sh simulations/<my_experiment>
cd simulations/<my_experiment>

# 2. Decouple config to customize
rm run.def && cp ../../configs/offline_1day/run.def ./run.def

# 3. Submit simulation
sbatch run_offline_1day.sh

# 4. Monitor
squeue -u $USER

# 5. Check termination
tail -n 10 slurm-*.out
```
