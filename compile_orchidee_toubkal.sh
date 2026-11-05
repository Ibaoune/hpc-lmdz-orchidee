#!/bin/bash
# ==============================================================================
# Script: compile_orchidee_toubkal.sh
# Description: Helper script to compile ORCHIDEE Offline on Toubkal HPC (UM6P)
# Author: M. El Aabaribaoune (UM6P)
# Affiliation: Mohammed VI Polytechnic University - IWRI Climate Team
# Contact: mohammad.elaabaribaoune@um6p.ma
# Date: 2026
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -d "${SCRIPT_DIR}/modeles" ]; then
    ROOT_MODIPSL="${SCRIPT_DIR}"
elif [ -d "${SCRIPT_DIR}/modipsl" ]; then
    ROOT_MODIPSL="${SCRIPT_DIR}/modipsl"
elif [ -d "${SCRIPT_DIR}/../modeles" ]; then
    ROOT_MODIPSL="$(cd "${SCRIPT_DIR}/.." && pwd)"
else
    ROOT_MODIPSL="/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl"
fi

NJOBS=8
CLEAN_BUILD=0
TARGET_DRIVER=1

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -j N         Number of parallel compile processes (default: 8)"
    echo "  -clean       Perform full clean before compiling"
    echo "  -h, --help   Display this help message"
    echo ""
    echo "Example:"
    echo "  $0 -j 16"
    echo "  $0 -clean -j 8"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j)
            NJOBS="$2"
            shift 2
            ;;
        -clean)
            CLEAN_BUILD=1
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

echo "=========================================================="
echo " ORCHIDEE Offline Compilation for Toubkal HPC (UM6P)"
echo " Author: M. El Aabaribaoune (UM6P)"
echo " modipsl root: ${ROOT_MODIPSL}"
echo " Parallel jobs: ${NJOBS}"
echo "=========================================================="

# 1. Environment Loading
ENV_SCRIPT="${ROOT_MODIPSL}/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env"
if [ ! -f "${ENV_SCRIPT}" ]; then
    echo "ERROR: Environment file not found: ${ENV_SCRIPT}"
    exit 0
fi
echo "==> Loading Toubkal HPC environment (foss-2018b)..."
source "${ENV_SCRIPT}"

# 2. Check/Build IOIPSL
IOIPSL_LIB="${ROOT_MODIPSL}/lib/libioipsl.a"
if [ ! -f "${IOIPSL_LIB}" ] || [ ${CLEAN_BUILD} -eq 1 ]; then
    echo "==> Building IOIPSL library..."
    cd "${ROOT_MODIPSL}/modeles/IOIPSL/src"
    if [ ${CLEAN_BUILD} -eq 1 ]; then
        gmake clean || true
    fi
    gmake -j "${NJOBS}"
    cd "${ROOT_MODIPSL}"
fi

# 3. Check XIOS Library
XIOS_LIB="${ROOT_MODIPSL}/modeles/XIOS/lib/libxios.a"
if [ ! -f "${XIOS_LIB}" ]; then
    echo "ERROR: XIOS library not found at: ${XIOS_LIB}"
    echo "Please ensure XIOS is built or contact M. El Aabaribaoune (UM6P)."
    exit 0
fi

# 4. Compile ORCHIDEE Offline Driver
echo "==> Compiling ORCHIDEE Offline..."
cd "${ROOT_MODIPSL}/modeles/ORCHIDEE"

FCM_ARGS="-xios2 -j ${NJOBS} -parallel mpi -prod -arch X64_TOUBKAL -driver"

if [ ${CLEAN_BUILD} -eq 1 ]; then
    echo "==> Cleaning previous ORCHIDEE build..."
    ./makeorchidee_fcm -clean -arch X64_TOUBKAL || true
fi

echo "==> Running: ./makeorchidee_fcm ${FCM_ARGS}"
./makeorchidee_fcm ${FCM_ARGS}

# 5. Verify Executables
mkdir -p "${ROOT_MODIPSL}/bin"
if [ -f "${ROOT_MODIPSL}/bin/orchidee_ol" ]; then
    echo "=========================================================="
    echo "SUCCESS: orchidee_ol successfully generated!"
    echo "Binary: ${ROOT_MODIPSL}/bin/orchidee_ol"
    ls -lh "${ROOT_MODIPSL}/bin/orchidee_ol"
    echo "=========================================================="
else
    echo "ERROR: orchidee_ol executable was not found in ${ROOT_MODIPSL}/bin/"
    exit 0
fi
