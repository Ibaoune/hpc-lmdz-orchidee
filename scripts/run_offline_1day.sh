#!/bin/bash
#SBATCH --job-name=orch_ol_1day
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
# Machine    : Toubkal (HPC UM6P)
# Exécutable : modipsl/bin/orchidee_ol (Mode XIOS intégré, 8 coeurs MPI)
#
# Emplacement des logs Slurm :
#   Slurm écrit slurm-<JOBID>.out et slurm-<JOBID>.err dans le répertoire
#   depuis lequel sbatch est exécuté (ou celui défini via --chdir).
#   Recommandation : lancez 'cd <RUN_DIR> && sbatch run_offline_1day.sh' pour
#   que les logs soient directement créés dans votre dossier d'exécution.
#
# Règles de sécurité :
#   - Aucun nettoyage automatique n'est exécuté sans l'option explicite --force.
#   - Le lancement est refusé si des fichiers de résultats précédents existent.
#   - Le répertoire de référence simulations/offline_1day/ est strictement
#     protégé contre toute purge, même avec l'option --force.
#   - Dans un répertoire personnel, la relance volontaire (--force) supprime
#     exclusivement les sorties de simulation connues.
# ==============================================================================

set -euo pipefail

# 1. Découverte de la racine du projet modipsl
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/../bin/orchidee_ol" ]]; then
    MODIPSL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
elif [[ -f "${SCRIPT_DIR}/../modipsl/bin/orchidee_ol" ]]; then
    MODIPSL_DIR="$(cd "${SCRIPT_DIR}/../modipsl" && pwd)"
else
    MODIPSL_DIR="/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl"
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
    RUN_DIR="${MODIPSL_DIR}/simulations/offline_1day"
fi

if [[ ! -d "${RUN_DIR}" ]]; then
    echo "ERREUR : Le répertoire d'exécution n'existe pas : ${RUN_DIR}" >&2
    echo "Initialisez-le d'abord avec : bash ${MODIPSL_DIR}/scripts/setup_offline_1day.sh ${RUN_DIR}" >&2
    exit 1
fi

cd "${RUN_DIR}"

REF_DIR_PHYSICAL="$(cd "${MODIPSL_DIR}/simulations/offline_1day" && pwd -P)"
RUN_DIR_PHYSICAL="$(pwd -P)"

# 3. Protection absolue du répertoire de référence partagé
# Toute nouvelle exécution y est strictement interdite, avec ou sans --force,
# même si les fichiers de résultats sont absents. Les simulations doivent utiliser
# des répertoires personnels.
if [[ "${RUN_DIR_PHYSICAL}" == "${REF_DIR_PHYSICAL}" ]]; then
    echo "=========================================================="
    echo "ERREUR DE SÉCURITÉ : Répertoire de référence partagé protégé !"
    echo "Le répertoire ${MODIPSL_DIR}/simulations/offline_1day est réservé"
    echo "à la référence officielle validée (Job 8228332)."
    echo "Toute nouvelle exécution y est strictement interdite."
    echo ""
    echo "Veuillez exécuter votre simulation dans un répertoire personnel :"
    echo "  bash ${MODIPSL_DIR}/scripts/setup_offline_1day.sh /chemin/vers/mon_repertoire"
    echo "=========================================================="
    exit 1
fi

# 4. Vérification de la présence de résultats antérieurs dans le répertoire personnel
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
        echo "Pour effectuer une nouvelle simulation, utilisez un autre répertoire :"
        echo "  bash ${MODIPSL_DIR}/scripts/setup_offline_1day.sh /chemin/vers/nouveau_repertoire"
        echo ""
        echo "Ou relancez volontairement dans ce répertoire personnel avec --force :"
        echo "  sbatch ${SCRIPT_DIR}/run_offline_1day.sh ${RUN_DIR} --force"
        echo "=========================================================="
        exit 1
    fi

    # Dans un répertoire personnel avec --force : suppression ciblée des sorties connues
    echo "Relance volontaire autorisée (--force) dans le répertoire personnel."
    echo "Suppression strictement limitée aux fichiers de résultats connus..."
    for f in "${KNOWN_OUTPUTS[@]}"; do
        rm -f "${f}"
    done
    rm -f out_orchidee_*
fi

# 4. Chargement de l'environnement Toubkal
ENV_FILE="${MODIPSL_DIR}/modipsl/modeles/ORCHIDEE/arch/arch-X64_TOUBKAL.env"
if [[ ! -f "${ENV_FILE}" ]]; then
    echo "ERREUR : Fichier d'environnement introuvable : ${ENV_FILE}" >&2
    exit 1
fi
# Désactiver temporairement set -u pour Lmod (/opt/ohpc/admin/lmod/lmod/init/bash)
set +u
source "${ENV_FILE}"
set -u

# 5. Paramétrage des limites système et threads
ulimit -s unlimited
export OMP_NUM_THREADS=1

# 6. Vérification de la présence de l'exécutable local
if [[ ! -x "./orchidee_ol" ]]; then
    echo "ERREUR : Exécutable ./orchidee_ol manquant ou non exécutable dans ${RUN_DIR}." >&2
    exit 1
fi

echo "=========================================================="
echo "Démarrage de la simulation ORCHIDEE Offline (1 jour)"
echo "Date                   : $(date)"
echo "Host                   : $(hostname)"
echo "Répertoire d'exécution : ${PWD}"
echo "Nombre de tâches MPI   : ${SLURM_NTASKS:-8}"
echo "Emplacement logs Slurm : slurm-${SLURM_JOB_ID:-manual}.out dans le dossier de soumission"
echo "=========================================================="

# 7. Lancement du modèle avec OpenMPI
mpirun -np ${SLURM_NTASKS:-8} ./orchidee_ol

echo "=========================================================="
echo "Fin de la simulation : $(date)"
echo "=========================================================="
