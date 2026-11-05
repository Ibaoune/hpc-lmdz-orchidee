#!/bin/bash
# ==============================================================================
# setup_offline_1day.sh
# ==============================================================================
# Automatisation non destructive de la configuration d'un répertoire
# d'exécution pour le test ORCHIDEE Offline d'une journée (1979-02-01).
#
# Règles de sécurité :
#   - Crée uniquement les liens symboliques manquants.
#   - Conserve les liens déjà corrects.
#   - Ne remplace et n'écrase JAMAIS un fichier régulier personnel existant.
#   - Signale clairement tout conflit ou fichier conservé.
#
# Usage:
#   bash scripts/setup_offline_1day.sh [TARGET_DIR]
#
# Par défaut : TARGET_DIR = simulations/offline_1day
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/../bin/orchidee_ol" ]]; then
    MODIPSL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
elif [[ -f "${SCRIPT_DIR}/../modipsl/bin/orchidee_ol" ]]; then
    MODIPSL_DIR="$(cd "${SCRIPT_DIR}/../modipsl" && pwd)"
else
    MODIPSL_DIR="/home/mohammad.elaabaribao/lustre/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl"
fi

TARGET_DIR="${1:-${MODIPSL_DIR}/simulations/offline_1day}"

echo "=========================================================="
echo "Configuration du répertoire d'exécution ORCHIDEE Offline"
echo "Répertoire racine modipsl  : ${MODIPSL_DIR}"
echo "Répertoire cible           : ${TARGET_DIR}"
echo "=========================================================="

mkdir -p "${TARGET_DIR}"
cd "${TARGET_DIR}"
TARGET_DIR_PHYSICAL="$(pwd -P)"
REF_DIR_PHYSICAL="$(cd "${MODIPSL_DIR}/simulations/offline_1day" && pwd -P)"
CONFLICT_COUNT=0

# Fonction non destructive de création de liens symboliques
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
        # Vérification si le lien pointe déjà vers la cible attendue
        if [[ "${current_target}" == "${target}" || "$(realpath -m "${current_target}")" == "$(realpath -m "${target}")" ]]; then
            # Lien déjà correct, préservation silencieuse
            return 0
        else
            echo "CONFLIT : Le lien '${link_name}' existe déjà mais pointe vers '${current_target}'." >&2
            echo "          Le lien existant est CONSERVÉ sans modification." >&2
            CONFLICT_COUNT=$((CONFLICT_COUNT + 1))
            return 0
        fi
    elif [[ -e "${link_name}" ]]; then
        echo "CONFLIT : Un fichier personnel '${link_name}' existe déjà." >&2
        echo "          Ce fichier personnel est STRICTEMENT PRÉSERVÉ (aucun écrasement)." >&2
            CONFLICT_COUNT=$((CONFLICT_COUNT + 1))
        return 0
    else
        ln -s "${target}" "${link_name}"
        echo "Lien créé : ${link_name} -> ${target}"
    fi
}

# 1. Lien vers l'exécutable compilé
EXEC="${MODIPSL_DIR}/bin/orchidee_ol"
if [[ ! -x "${EXEC}" ]]; then
    echo "ERREUR : Exécutable introuvable ou non exécutable : ${EXEC}" >&2
    exit 1
fi
safe_link "${EXEC}" "orchidee_ol"

# 2. Liens vers les fichiers de configuration de référence (configs/offline_1day/)
CONFIGS_DIR="${MODIPSL_DIR}/configs/offline_1day"
for cfg in run.def iodef.xml context_orchidee.xml context_input_orchidee.xml field_def_orchidee.xml file_def_orchidee.xml; do
    if [[ ! -f "${CONFIGS_DIR}/${cfg}" ]]; then
        echo "ERREUR : Fichier de configuration de référence manquant : ${CONFIGS_DIR}/${cfg}" >&2
        exit 1
    fi
    safe_link "${CONFIGS_DIR}/${cfg}" "${cfg}"
done

# 3. Lien vers le script Slurm d'exécution
safe_link "${SCRIPT_DIR}/run_offline_1day.sh" "run_offline_1day.sh"

# 4. Liens vers les données de surface et conditions aux limites partagées
IGCM_SRF="/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/igcmg/IGCM/SRF"

safe_link "${IGCM_SRF}/SOIL/soils_param.nc" "soils_param.nc"
safe_link "${IGCM_SRF}/reftemp.nc" "reftemp.nc"
safe_link "${IGCM_SRF}/ROUTING/routing.nc" "routing.nc"
safe_link "${IGCM_SRF}/ROUTING/floodplains.nc" "floodplains.nc"
safe_link "${IGCM_SRF}/cartepente2d_15min.nc" "cartepente2d_15min.nc"
safe_link "${IGCM_SRF}/albedo/alb_bg_modisopt_2D_ESA_v2.nc" "alb_bg.nc"
safe_link "${IGCM_SRF}/PFTMAPS/CMIP6/ESA-LUH2v2/historical/15PFT.v1/PFTmap_1979.nc" "PFTmap.nc"

# 5. Lien vers le restart initial SECHIBA (31 janvier 1979, strictement en lecture)
HIST_RESTART="/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/ARCHIVE/Simulations/LMDZ/r35/hist/SRF/Restart/sim35km-HIST_19790131_sechiba_rest.nc"
if [[ ! -f "${HIST_RESTART}" ]]; then
    echo "ERREUR : Restart historique introuvable : ${HIST_RESTART}" >&2
    exit 1
fi
safe_link "${HIST_RESTART}" "sechiba_rest_in.nc"

# 6. Forçage atmosphérique journalier
FORCING_FILE="forcing_19790201.nc"
if [[ -e "${FORCING_FILE}" ]]; then
    echo "Fichier de forçage déjà présent : ${FORCING_FILE} ($(du -h "${FORCING_FILE}" | cut -f1))"
else
    REF_FORCING="${MODIPSL_DIR}/simulations/offline_1day/${FORCING_FILE}"
    if [[ -f "${REF_FORCING}" && "${TARGET_DIR_PHYSICAL}" != "${REF_DIR_PHYSICAL}" ]]; then
        safe_link "${REF_FORCING}" "${FORCING_FILE}"
    else
        echo "NOTE : Le fichier de forçage ${FORCING_FILE} n'est pas encore présent dans ${TARGET_DIR}."
        echo "       Générez-le via la commande :"
        echo "         python3 ${MODIPSL_DIR}/scripts/prepare_forcing_1day.py --output ${TARGET_DIR}/${FORCING_FILE}"
    fi
fi

# 7. Bilan et instructions
echo "----------------------------------------------------------"
if [[ ${CONFLICT_COUNT} -gt 0 ]]; then
    echo "ATTENTION : ${CONFLICT_COUNT} conflit(s) détecté(s) !"
    echo "Certains liens pointent vers des cibles personnalisées ou des fichiers"
    echo "personnels ont été conservés. La configuration n'est pas standard."
    echo "Vérifiez vos fichiers dans ${TARGET_DIR} avant de lancer la simulation."
else
    echo "Succès : Répertoire ${TARGET_DIR} vérifié et prêt (configuration standard)."
    echo "Pour lancer la simulation Slurm dans ce répertoire :"
    echo "  cd ${TARGET_DIR} && sbatch run_offline_1day.sh"
fi
echo "=========================================================="
