#!/bin/bash
# ==============================================================================
# setup_offline_1day_gcc12.sh
# ==============================================================================
# Automatisation de la configuration d'un répertoire d'exécution pour le test
# ORCHIDEE Offline GCC 12 d'une journée (1979-02-01).
#
# Usage:
#   bash scripts/setup_offline_1day_gcc12.sh [TARGET_DIR]
#
# Par défaut : TARGET_DIR = simulations/offline_1day_gcc12
# ==============================================================================

set -euo pipefail

# 1. Découverte dynamique de la racine modipsl
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -d "${SCRIPT_DIR}/../modeles" ]]; then
    MODIPSL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
elif [[ -d "${SCRIPT_DIR}/modeles" ]]; then
    MODIPSL_DIR="${SCRIPT_DIR}"
elif [[ -d "${PWD}/modeles" ]]; then
    MODIPSL_DIR="${PWD}"
else
    echo "ERREUR : Impossible de localiser la racine modipsl contenant le dossier modeles/." >&2
    exit 1
fi

TARGET_DIR="${1:-${MODIPSL_DIR}/simulations/offline_1day_gcc12}"

echo "=========================================================="
echo "Configuration du répertoire d'exécution ORCHIDEE Offline (GCC 12)"
echo "Répertoire racine modipsl  : ${MODIPSL_DIR}"
echo "Répertoire cible           : ${TARGET_DIR}"
echo "=========================================================="

mkdir -p "${TARGET_DIR}"
cd "${TARGET_DIR}"
TARGET_DIR_PHYSICAL="$(pwd -P)"
CONFLICT_COUNT=0

safe_link() {
    local target="$1"
    local link_name="$2"

    if [[ ! -e "${target}" ]]; then
        echo "ERREUR : Cible introuvable : ${target}" >&2
        return 1
    fi

    if [[ -L "${link_name}" ]]; then
        local current_target
        current_target="$(readlink "${link_name}")"
        if [[ "${current_target}" == "${target}" || "$(realpath -m "${current_target}")" == "$(realpath -m "${target}")" ]]; then
            return 0
        else
            echo "CONFLIT : Le lien '${link_name}' existe déjà mais pointe vers '${current_target}'." >&2
            CONFLICT_COUNT=$((CONFLICT_COUNT + 1))
            return 0
        fi
    elif [[ -e "${link_name}" ]]; then
        echo "CONFLIT : Un fichier personnel '${link_name}' existe déjà. Conservé." >&2
        CONFLICT_COUNT=$((CONFLICT_COUNT + 1))
        return 0
    else
        ln -s "${target}" "${link_name}"
        echo "Lien créé : ${link_name} -> ${target}"
    fi
}

# 1. Lien vers le nouvel exécutable GCC 12
EXEC="${MODIPSL_DIR}/bin/orchidee_ol_gcc12"
if [[ ! -x "${EXEC}" ]]; then
    echo "ERREUR : Exécutable GCC 12 introuvable ou non exécutable : ${EXEC}" >&2
    echo "Veuillez d'abord compiler ORCHIDEE Offline avec :" >&2
    echo "  ./compile_orchidee_gcc12.sh" >&2
    exit 1
fi
safe_link "${EXEC}" "orchidee_ol"

# 2. Liens vers les fichiers de configuration de référence
CONFIGS_DIR="${MODIPSL_DIR}/configs/offline_1day"
for cfg in run.def iodef.xml context_orchidee.xml context_input_orchidee.xml field_def_orchidee.xml file_def_orchidee.xml; do
    if [[ ! -f "${CONFIGS_DIR}/${cfg}" ]]; then
        echo "ERREUR : Fichier de configuration manquant : ${CONFIGS_DIR}/${cfg}" >&2
        exit 1
    fi
    safe_link "${CONFIGS_DIR}/${cfg}" "${cfg}"
done

# 3. Lien vers le script Slurm d'exécution GCC 12
safe_link "${SCRIPT_DIR}/run_offline_1day_gcc12.sh" "run_offline_1day.sh"

# 4. Liens vers les données de surface et conditions aux limites partagées (lecture seule)
IGCM_SRF="/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/igcmg/IGCM/SRF"

safe_link "${IGCM_SRF}/SOIL/soils_param.nc" "soils_param.nc"
safe_link "${IGCM_SRF}/reftemp.nc" "reftemp.nc"
safe_link "${IGCM_SRF}/ROUTING/routing.nc" "routing.nc"
safe_link "${IGCM_SRF}/ROUTING/floodplains.nc" "floodplains.nc"
safe_link "${IGCM_SRF}/cartepente2d_15min.nc" "cartepente2d_15min.nc"
safe_link "${IGCM_SRF}/albedo/alb_bg_modisopt_2D_ESA_v2.nc" "alb_bg.nc"
safe_link "${IGCM_SRF}/PFTMAPS/CMIP6/ESA-LUH2v2/historical/15PFT.v1/PFTmap_1979.nc" "PFTmap.nc"

# 5. Lien vers le restart initial SECHIBA (31 janvier 1979)
HIST_RESTART="/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist/SRF/Restart/sim35km-HIST_19790131_sechiba_rest.nc"
if [[ ! -f "${HIST_RESTART}" ]]; then
    echo "ERREUR : Restart historique introuvable : ${HIST_RESTART}" >&2
    exit 1
fi
safe_link "${HIST_RESTART}" "sechiba_rest_in.nc"

# 6. Forçage atmosphérique journalier (1er février 1979)
FORCING_FILE="forcing_19790201.nc"
LOCAL_FORCING="${MODIPSL_DIR}/simulations/offline_1day/${FORCING_FILE}"
SHARED_FORCING="/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/simulations/offline_1day/${FORCING_FILE}"

if [[ -f "${LOCAL_FORCING}" ]]; then
    safe_link "${LOCAL_FORCING}" "${FORCING_FILE}"
elif [[ -f "${SHARED_FORCING}" ]]; then
    safe_link "${SHARED_FORCING}" "${FORCING_FILE}"
else
    echo "ERREUR : Forçage de référence introuvable." >&2
    echo "Emplacements vérifiés :" >&2
    echo "  - ${LOCAL_FORCING}" >&2
    echo "  - ${SHARED_FORCING}" >&2
    exit 1
fi

echo "----------------------------------------------------------"
if [[ ${CONFLICT_COUNT} -gt 0 ]]; then
    echo "ATTENTION : ${CONFLICT_COUNT} conflit(s) détecté(s)."
else
    echo "Succès : Répertoire ${TARGET_DIR} vérifié et prêt pour GCC 12."
    echo "Pour lancer la simulation :"
    echo "  cd ${TARGET_DIR} && sbatch run_offline_1day.sh"
fi
echo "=========================================================="
