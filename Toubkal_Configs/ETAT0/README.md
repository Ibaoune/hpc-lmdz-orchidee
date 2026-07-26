# Création de l'État Initial (ETAT0)

Ce dossier contient la configuration utilisée pour générer les états initiaux (`start.nc`, `startphy.nc`) et les forçages aux limites (`limit.nc`) pour tes simulations AMIP LMDZ.

## Configuration Clé (`config.card`)
- **ExpType** : `LMDZ/CREATE_amip`
- **Exécutable** : `ce0l_${ResolAtm}_${OptMode}.e` (au lieu du gcm classique)
- **PeriodLength** : `1Y`
- **Restart** : `n`

## Fonctionnement
Lors de l'exécution, libIGCM va :
1. Télécharger les forçages globaux (SST, banquise, etc.).
2. Adapter ces champs à la grille locale spécifiée (`384x220x79`).
3. Stocker les résultats finaux dans l'espace d'archivage `IGCM` défini par `libIGCM` (sur le cluster Toubkal, c'était dans `igcmg/`).

## Utilisation
Ces fichiers (`start.nc` et `limit.nc`) devront être renseignés dans la section `[InitialStateFiles]` et `[BoundaryFiles]` du fichier `COMP/lmdz.card` de tes futures simulations (ex. HIST).
