# Simulation Future (SSP585)

Ce dossier contient la configuration de référence pour ta simulation future scénario SSP585, qui s'étend de **2015 à 2100**.

## Configuration Clé (`config.card`)
- **ExpType** : `LMDZOR/amip` (ou apparenté)
- **DateBegin** : `2015-01-01`
- **DateEnd** : `2100-12-31`
- **PeriodLength** : `1M` (avance mois par mois)

## Continuité avec HIST (`Restart=y`)
Le point crucial des simulations scénarios est de garantir la parfaite continuité climatique avec le dernier jour de la simulation Historique.

Dans la section `[Restarts]` de `config.card` :
```ini
Restart=y
RestartDate=2014-12-31           # La date de fin de la simulation HIST
RestartJobName=sim35km-HIST      # Le JobName exact de la simulation HIST
RestartPath=...                  # Le chemin où se trouvent les restarts de sim35km-HIST
```

Tout comme pour le scénario SSP245, libIGCM ira chercher l'état final de l'année 2014 de la simulation HIST, et commencera son intégration temporelle en appliquant le forçage extrême SSP585.
