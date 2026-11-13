#!/bin/bash
#SBATCH --job-name=orch_ol_gcc12
#SBATCH --account=CLIMAT-UM6P-ST-IWRI-7KSIFKVWKUY-DEFAULT-CPU
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:30:00
#SBATCH --hint=nomultithread
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# ==============================================================================
# Script de soumission Slurm pour le test ORCHIDEE Offline 1 jour (1979-02-01)
# Compilateur : GCC 12.3.0 / OpenMPI 4.1.5 (gompi-2023a)
# Machine     : Toubkal (HPC UM6P)
# Exécutable  : modipsl/bin/orchidee_ol_gcc12 (Mode XIOS intégré, 8 rangs MPI)
# ==============================================================================

set -euo pipefail

# 1. Découverte dynamique de la racine modipsl
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -d "${SCRIPT_DIR}/../modeles" ]]; then
    MODIPSL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
elif [[ -d "${SCRIPT_DIR}/modeles" ]]; then
    MODIPSL_DIR="${SCRIPT_DIR}"
elif [[ -d "${PWD}/../../modeles" ]]; then
    MODIPSL_DIR="$(cd "${PWD}/../.." && pwd)"
elif [[ -d "${PWD}/modeles" ]]; then
    MODIPSL_DIR="${PWD}"
else
    echo "ERREUR : Impossible de localiser la racine modipsl contenant le dossier modeles/." >&2
    exit 1
fi

# 2. Détermination du répertoire d'exécution (RUN_DIR)
FORCE_RERUN="${FORCE_RERUN:-0}"
TARGET_ARG=""

for arg in "$@"; do
    if [[ "$arg" == "--force" || "$arg" == "-f" ]]; then
        FORCE_RERUN=1
    elif [[ -z "${TARGET_ARG}" && "${arg}" != --* ]]; then
        TARGET_ARG="${arg}"
    fi
done

if [[ -n "${TARGET_ARG}" ]]; then
    RUN_DIR="${TARGET_ARG}"
elif [[ -n "${SIM_DIR:-}" ]]; then
    RUN_DIR="${SIM_DIR}"
elif [[ -f "./run.def" && -f "./orchidee_ol" ]]; then
    RUN_DIR="${PWD}"
else
    RUN_DIR="${MODIPSL_DIR}/simulations/offline_1day_gcc12"
fi

if [[ ! -d "${RUN_DIR}" ]]; then
    echo "ERREUR : Le répertoire d'exécution n'existe pas : ${RUN_DIR}" >&2
    echo "Initialisez-le d'abord avec : bash ${MODIPSL_DIR}/scripts/setup_offline_1day_gcc12.sh ${RUN_DIR}" >&2
    exit 1
fi

cd "${RUN_DIR}"

REF_DIR_PHYSICAL="$(cd "${MODIPSL_DIR}/simulations/offline_1day" 2>/dev/null && pwd -P || echo "/nonexistent_ref")"
SHARED_REF1="/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/simulations/offline_1day"
SHARED_REF2="/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/simulations/offline_1day_reproduction"
RUN_DIR_PHYSICAL="$(pwd -P)"

# 3. Protection absolue des répertoires de référence partagés
if [[ "${RUN_DIR_PHYSICAL}" == "${REF_DIR_PHYSICAL}" || \
      "${RUN_DIR_PHYSICAL}" == "$(cd "${SHARED_REF1}" 2>/dev/null && pwd -P || echo "")" || \
      "${RUN_DIR_PHYSICAL}" == "$(cd "${SHARED_REF2}" 2>/dev/null && pwd -P || echo "")" ]]; then
    echo "=========================================================="
    echo "ERREUR DE SÉCURITÉ : Répertoire de référence partagé protégé !"
    echo "Le répertoire ${RUN_DIR} est réservé à la référence certifiée."
    echo "Toute nouvelle exécution y est strictement interdite."
    echo "=========================================================="
    exit 1
fi

# 4. Vérification de la présence de résultats antérieurs dans le répertoire
KNOWN_OUTPUTS=(
    "sechiba_history.nc"
    "sechiba_out_2.nc"
    "sechiba_rest_out.nc"
    "driver_rest_out.nc"
    "used_run.def"
    "Load_balance_orchidee.dat"
    "out_orchidee_0000"
)

FOUND_OUTPUTS=()
for out_name in "${KNOWN_OUTPUTS[@]}"; do
    if [[ -f "${out_name}" ]]; then
        FOUND_OUTPUTS+=("${out_name}")
    fi
done

if [[ ${#FOUND_OUTPUTS[@]} -gt 0 ]]; then
    if [[ "${FORCE_RERUN}" != "1" ]]; then
        echo "=========================================================="
        echo "ERREUR : Des fichiers de résultats précédents existent dans :"
        echo "  ${RUN_DIR}"
        echo "Fichiers détectés : ${FOUND_OUTPUTS[*]}"
        echo ""
        echo "Le lancement est refusé pour préserver ces résultats."
        echo "Pour forcer la relance dans ce répertoire : sbatch ${0} --force"
        echo "=========================================================="
        exit 1
    fi

    echo "Relance autorisée (--force). Nettoyage des sorties précédentes..."
    for f in "${KNOWN_OUTPUTS[@]}"; do
        rm -f "${f}"
    done
    rm -f out_orchidee_*
fi

# 5. Chargement de l'environnement GCC 12 / gompi-2023a
ENV_FILE="${MODIPSL_DIR}/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL_GCC12.env"
if [[ ! -f "${ENV_FILE}" ]]; then
    echo "ERREUR : Fichier d'environnement introuvable : ${ENV_FILE}" >&2
    exit 1
fi

set +u
source "${ENV_FILE}"
set -u

# 6. Paramétrage des limites système et threads
ulimit -s unlimited
export OMP_NUM_THREADS=1

# 7. Vérification de l'exécutable local
if [[ ! -x "./orchidee_ol" ]]; then
    echo "ERREUR : Exécutable ./orchidee_ol manquant ou non exécutable dans ${RUN_DIR}." >&2
    exit 1
fi

echo "=========================================================="
echo "Démarrage de la simulation ORCHIDEE Offline GCC 12 (1 jour)"
echo "Date                   : $(date)"
echo "Host                   : $(hostname)"
echo "Répertoire d'exécution : ${PWD}"
echo "Exécutable             : $(readlink -f ./orchidee_ol)"
echo "Nombre de rangs MPI    : ${SLURM_NTASKS:-8}"
echo "Toolchain              : GCC 12.3.0 / OpenMPI 4.1.5 (gompi-2023a)"
echo "=========================================================="

# 8. Lancement avec OpenMPI
mpirun -np ${SLURM_NTASKS:-8} ./orchidee_ol

echo "=========================================================="
echo "Fin normale de la simulation : $(date)"
echo "=========================================================="
