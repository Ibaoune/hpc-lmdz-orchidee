#!/bin/bash
# ==============================================================================
# Script: compile_orchidee_gcc12.sh
# Description: Fully independent build of ORCHIDEE Offline with modern GCC 12
#              suite on Toubkal HPC (IOIPSL, XIOS, ORCHIDEE).
# Toolchain: GCC 12.3.0 / OpenMPI 4.1.5 / NetCDF 4.9.2 / NetCDF-Fortran 4.6.1 / HDF5 1.14.0
# Author: M. El Aabaribaoune (IWRI / UM6P)
# Affiliation: Mohammed VI Polytechnic University - IWRI Climate Team
# Date: November 2026
# ==============================================================================

set -euo pipefail

# 1. Dynamic discovery of modipsl root (completely independent of user path)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -d "${SCRIPT_DIR}/modeles" ]; then
    ROOT_MODIPSL="${SCRIPT_DIR}"
elif [ -d "${SCRIPT_DIR}/../modeles" ]; then
    ROOT_MODIPSL="$(cd "${SCRIPT_DIR}/.." && pwd)"
elif [ -d "${PWD}/modeles" ]; then
    ROOT_MODIPSL="${PWD}"
else
    echo "ERROR: Unable to locate modipsl root containing modeles/ directory." >&2
    echo "Please execute this script from within the modipsl repository root." >&2
    exit 1
fi

NJOBS=8
CLEAN_BUILD=0

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -j N         Number of parallel compile processes (default: 8)"
    echo "  -clean       Perform full clean before compiling (rebuilds IOIPSL & XIOS)"
    echo "  -h, --help   Display this help message"
    echo ""
    echo "Example:"
    echo "  $0 -j 8"
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
echo " ORCHIDEE Offline Compilation (GCC 12 / gompi-2023a)"
echo " Toubkal HPC — Mohammed VI Polytechnic University (UM6P)"
echo " Workspace : ${ROOT_MODIPSL}"
echo " Jobs      : ${NJOBS}"
echo "=========================================================="

# 2. Environment Loading
ENV_SCRIPT="${ROOT_MODIPSL}/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL_GCC12.env"
if [ ! -f "${ENV_SCRIPT}" ]; then
    echo "ERROR: Environment file not found: ${ENV_SCRIPT}" >&2
    exit 1
fi

echo "==> Loading Toubkal HPC modern environment (gompi-2023a)..."
# Disable set -u temporarily for Lmod
set +u
source "${ENV_SCRIPT}"
set -u

echo "Compiler       : $(which mpif90)"
echo "GCC Version    : $(gcc --version | head -n 1)"
echo "NetCDF-Fortran : $(nf-config --version 2>/dev/null || echo 'loaded')"
echo "HDF5 Version   : $(h5pcc -showconfig 2>/dev/null | grep -i "HDF5 Version" | head -n 1 | awk '{print $3}' || echo 'loaded')"

# 3. Pre-flight verification of build tools and bundled FCM distributions
FCM_XIOS="${ROOT_MODIPSL}/modeles/XIOS/tools/FCM_OLD/bin/fcm"
FCM_ORCH="${ROOT_MODIPSL}/modeles/ORCHIDEE/tools/FCM_V1.2/bin/fcm"

if [ ! -f "${FCM_XIOS}" ]; then
    echo "ERROR: XIOS FCM build tool not found: ${FCM_XIOS}" >&2
    echo "Ensure modeles/XIOS/tools/FCM_OLD is bundled in the repository." >&2
    exit 1
fi

if [ ! -f "${FCM_ORCH}" ]; then
    echo "ERROR: ORCHIDEE FCM build tool not found: ${FCM_ORCH}" >&2
    echo "Ensure modeles/ORCHIDEE/tools/FCM_V1.2 is bundled in the repository." >&2
    exit 1
fi

chmod +x "${FCM_XIOS}" "${FCM_ORCH}"
chmod +x "${ROOT_MODIPSL}/modeles/XIOS/make_xios"
chmod +x "${ROOT_MODIPSL}/modeles/ORCHIDEE/makeorchidee_fcm"
mkdir -p "${ROOT_MODIPSL}/modeles/XIOS/.void_dir"

# 4. Check and adapt IOIPSL Makefile for local directory and GCC 12
IOIPSL_MAKEFILE="${ROOT_MODIPSL}/modeles/IOIPSL/src/Makefile"
if [ ! -f "${IOIPSL_MAKEFILE}" ]; then
    echo "ERROR: IOIPSL Makefile not found: ${IOIPSL_MAKEFILE}" >&2
    exit 1
fi

# Dynamically adapt MAKE_DIR to local modipsl directory, preserving the Makefile safety check
echo "==> Adapting IOIPSL Makefile MAKE_DIR to local workspace..."
sed -i "s|^MAKE_DIR = .*|MAKE_DIR = ${ROOT_MODIPSL}/modeles/IOIPSL/src|" "${IOIPSL_MAKEFILE}"

if grep -q "4.4.4-foss-2018b" "${IOIPSL_MAKEFILE}" 2>/dev/null || ! grep -q -- "-fallow-argument-mismatch" "${IOIPSL_MAKEFILE}" 2>/dev/null; then
    echo "==> Adapting IOIPSL Makefile compiler flags for GCC 12 (gompi-2023a)..."
    PATCH_FILE="${ROOT_MODIPSL}/scripts/patches/ioipsl_makefile_toubkal.patch"
    if [ -f "${PATCH_FILE}" ]; then
        patch -p1 -d "${ROOT_MODIPSL}" < "${PATCH_FILE}" || {
            echo "WARNING: Patch application returned non-zero. Verifying Makefile directly..."
        }
    fi
fi

# Explicit check for GCC 12 compatibility in IOIPSL Makefile
if grep -q "4.4.4-foss-2018b" "${IOIPSL_MAKEFILE}" 2>/dev/null || ! grep -q -- "-fallow-argument-mismatch" "${IOIPSL_MAKEFILE}" 2>/dev/null; then
    echo "ERROR: IOIPSL Makefile is incompatible with GCC 12 (-fallow-argument-mismatch missing or legacy NetCDF path present)." >&2
    exit 1
fi

# 5. Check/Build IOIPSL
IOIPSL_LIB="${ROOT_MODIPSL}/lib/libioipsl.a"
if [ ! -f "${IOIPSL_LIB}" ] || [ ${CLEAN_BUILD} -eq 1 ]; then
    echo "==> Building IOIPSL library with GCC 12..."
    mkdir -p "${ROOT_MODIPSL}/lib"
    cd "${ROOT_MODIPSL}/modeles/IOIPSL/src"
    if [ ${CLEAN_BUILD} -eq 1 ]; then
        gmake clean || true
    fi
    gmake || {
        echo "ERROR: IOIPSL compilation failed." >&2
        exit 1
    }
    cd "${ROOT_MODIPSL}"
else
    echo "==> IOIPSL library already present: ${IOIPSL_LIB}"
fi

# 6. Check/Build XIOS Library
XIOS_LIB="${ROOT_MODIPSL}/modeles/XIOS/lib/libxios.a"
if [ ! -f "${XIOS_LIB}" ] || [ ${CLEAN_BUILD} -eq 1 ]; then
    echo "==> Building XIOS library with GCC 12..."
    cd "${ROOT_MODIPSL}/modeles/XIOS"
    if [ ${CLEAN_BUILD} -eq 1 ]; then
        ./make_xios --arch X64_TOUBKAL_GCC12 --clean || true
    fi
    ./make_xios --arch X64_TOUBKAL_GCC12 --job "${NJOBS}" || {
        echo "ERROR: XIOS compilation failed." >&2
        exit 1
    }
    cd "${ROOT_MODIPSL}"
else
    echo "==> XIOS library already present: ${XIOS_LIB}"
fi

# 7. Check Active Locks in ORCHIDEE
LOCK_FILE="${ROOT_MODIPSL}/modeles/ORCHIDEE/build/fcm.bld.lock"
if [ -f "${LOCK_FILE}" ]; then
    LOCK_PID="$(grep -o '[0-9]\+' "${LOCK_FILE}" 2>/dev/null | head -n 1 || echo "")"
    if [ -n "${LOCK_PID}" ] && kill -0 "${LOCK_PID}" 2>/dev/null; then
        echo "ERROR: Active compilation lock detected in ${LOCK_FILE} (PID ${LOCK_PID} is running)." >&2
        echo "Please wait for the current compilation to finish or terminate PID ${LOCK_PID}." >&2
        exit 1
    elif [ ${CLEAN_BUILD} -eq 1 ]; then
        echo "==> Removing stale lock (${LOCK_FILE}) as requested by -clean..."
        rm -f "${LOCK_FILE}"
    else
        echo "WARNING: Stale lockfile detected (${LOCK_FILE})."
        echo "If no compilation is running, use -clean to remove stale artifacts or remove manually:"
        echo "  rm -f ${LOCK_FILE}"
        exit 1
    fi
fi

# 8. Compile ORCHIDEE Offline Driver
echo "==> Compiling ORCHIDEE Offline (dim2_driver / orchidee_ol)..."
cd "${ROOT_MODIPSL}/modeles/ORCHIDEE"

FCM_ARGS="-xios2 -j ${NJOBS} -parallel mpi -prod -arch X64_TOUBKAL_GCC12 -driver"

if [ ${CLEAN_BUILD} -eq 1 ]; then
    echo "==> Cleaning previous ORCHIDEE build..."
    ./makeorchidee_fcm -clean -arch X64_TOUBKAL_GCC12 || true
fi

echo "==> Running: ./makeorchidee_fcm ${FCM_ARGS}"
./makeorchidee_fcm ${FCM_ARGS} || {
    echo "ERROR: ORCHIDEE FCM compilation failed." >&2
    exit 1
}

# 9. Install GCC 12 executable in local bin/
mkdir -p "${ROOT_MODIPSL}/bin"

COMPILED_EXE="${ROOT_MODIPSL}/modeles/ORCHIDEE/build/bin/dim2_driver.exe"
if [ -f "${COMPILED_EXE}" ]; then
    echo "==> Installing GCC 12 executable into bin/orchidee_ol_gcc12..."
    cp -f "${COMPILED_EXE}" "${ROOT_MODIPSL}/bin/orchidee_ol_gcc12"
    chmod +x "${ROOT_MODIPSL}/bin/orchidee_ol_gcc12"
    
    echo "=========================================================="
    echo "SUCCESS: Compilation completed successfully!"
    echo "Executable : ${ROOT_MODIPSL}/bin/orchidee_ol_gcc12"
    ls -lh "${ROOT_MODIPSL}/bin/orchidee_ol_gcc12"
    echo "=========================================================="
else
    echo "ERROR: Compiled executable not found at: ${COMPILED_EXE}" >&2
    exit 1
fi
