# Configurations de Production Pluriannuelles (`experiments/`)

Ce dossier conserve les fichiers de configuration de référence utilisés pour les simulations longues multi-décennales sur Toubkal, importés depuis la branche `origin/master` (commit `61124b2`).

## Sous-dossiers :
1. **[`etat0/`](etat0/)** : Génération des états initiaux (`start.nc`, `startphy.nc`) et conditions aux limites mensuelles (`limit.nc`) pour le domaine régional à partir des forçages globaux AMIP (`ce0l.e`).
2. **[`hist/`](hist/)** : Simulation historique continue de 1979 à 2014 (`PeriodLength = 1M`, forçages solaires et gaz à effet de serre historiques CMIP6, sorties HF activées).
3. **[`ssp245/`](ssp245/)** : Scénario climatique futur moyen 2015–2100 (`PeriodLength = 1M`, forçages SSP2-4.5).
4. **[`ssp585/`](ssp585/)** : Scénario climatique futur extrême 2015–2100 (`PeriodLength = 1M`, forçages SSP5-8.5).

Toutes les valeurs originales sont rigoureusement préservées pour assurer la traçabilité historique des résultats de l'équipe.
