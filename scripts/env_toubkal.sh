#!/bin/bash
# ==============================================================================
# env_toubkal.sh
# Environnement d'exécution HPC Toubkal (UM6P) pour LMDZ - ORCHIDEE - XIOS
# Compatible ksh, bash, sh
# ==============================================================================

# 1. Détermination de la racine modipsl
if [ -n "${MODIPSL:-}" ] && [ -d "${MODIPSL}" ]; then
    _MODIPSL="${MODIPSL}"
elif [ -d "/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl" ]; then
    _MODIPSL="/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl"
else
    _MODIPSL="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../../.." && pwd)"
fi

# 2. Initialisation Lmod si nécessaire (support des sous-shells non-interactifs / ksh)
if [ -z "${LMOD_CMD:-}" ] && [ -f /srv/software/lmod/lmod/libexec/lmod ]; then
    export LMOD_CMD=/srv/software/lmod/lmod/libexec/lmod
    module() { eval "$($LMOD_CMD sh "$@")"; }
fi

# 3. Purge et ajout de l'arborescence des modules foss-2018b
module purge > /dev/null 2>&1 || true
module use "${_MODIPSL}/modules/foss-2018b"

# 4. Chargement des modules compatibles (compilateur, MPI, I/O, math)
module load Perl/5.42.0-GCCcore-15.2.0 > /dev/null 2>&1 || module load Perl > /dev/null 2>&1 || true
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1-GCC-7.3.0-2.30
module load HDF5/1.10.2-foss-2018b
module load netCDF/4.6.1-foss-2018b
module load netCDF-Fortran/4.4.4-foss-2018b
module load FFTW/3.3.8-gompi-2018b

# 5. Export explicite pour propagation MPI
export LD_LIBRARY_PATH
export PATH
