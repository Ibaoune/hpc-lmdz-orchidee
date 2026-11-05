#!/bin/bash
# ==============================================================================
# Script: submit_standalone_template.sh
# Description: Slurm batch submission template for ORCHIDEE Offline on Toubkal HPC
# Author: M. El Aabaribaoune (UM6P)
# Affiliation: Mohammed VI Polytechnic University - IWRI Climate Team
# Contact: mohammad.elaabaribaoune@um6p.ma
# Date: 2026
# ==============================================================================
#
#SBATCH --job-name=orchidee_ol
#SBATCH --account=CLIMAT-UM6P-ST-IWRI-7KSIFKVWKUY-DEFAULT-CPU
#SBATCH --ntasks=8
#SBATCH --time=02:00:00
#SBATCH --hint=nomultithread
#SBATCH --output=run_orch_%j.out
#SBATCH --error=run_orch_%j.err

set -e

# 1. Base Paths
SHARED_MODELS="/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models"
MODIPSL="${SHARED_MODELS}/modipsl"
BIN="${MODIPSL}/bin/orchidee_ol"

echo "Job started on $(hostname) at $(date)"
echo "Working directory: $(pwd)"

# 2. Load the Toubkal HPC environment (foss-2018b)
source "${MODIPSL}/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env"

# 3. Sanity check: Ensure required input/control files exist in execution directory
for file in run.def iodef.xml context_orchidee.xml field_def_orchidee.xml file_def_orchidee.xml; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Missing required control file: $file"
        echo "Please copy XML files from ${MODIPSL}/modeles/ORCHIDEE/src_xml/ and run.def from ${MODIPSL}/config/ORCHIDEE_OL/OOL_SEC/PARAM/"
        exit 1
    fi
done

# 4. Launch simulation with Slurm srun
echo "Launching orchidee_ol with ${SLURM_NTASKS} MPI tasks..."
srun "${BIN}"

echo "Job completed successfully at $(date)"
