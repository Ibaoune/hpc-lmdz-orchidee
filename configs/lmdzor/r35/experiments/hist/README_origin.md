# Simulation Historique (HIST)

Ce dossier contient la configuration de référence pour ta simulation historique principale sur Toubkal, qui s'étend de **1979 à 2014**.

## Configuration Clé (`config.card`)
- **ExpType** : `LMDZOR/amip`
- **Exécutable** : `gcm_${ResolAtm}_${OptMode}.e`
- **DateBegin** : `1979-01-01`
- **DateEnd** : `2014-12-31`
- **PeriodLength** : `1M` (avance mois par mois)
- **Restart** : `n` (on démarre avec les états initiaux créés par l'étape ETAT0)

## Spécificité (`COMP/lmdz.card`)
Dans `COMP/lmdz.card`, la section `[InitialStateFiles]` doit absolument pointer vers les fichiers `start.nc` et `startphy.nc` générés par l'étape de création de l'état initial (`CREATE_amip`). De même pour la section `[BoundaryFiles]` avec les `limit.nc`.

## Exécution
libIGCM lance les mois les uns après les autres. À la fin de chaque mois, il génère des fichiers de Restart, qui sont réinjectés automatiquement pour le mois suivant. Le dernier Restart généré (2014-12-31) servira de point de départ pour les simulations SSP.
