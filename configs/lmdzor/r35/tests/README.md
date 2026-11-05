# Bancs d'Essai Couplés 3 Jours (`tests/`)

Ce dossier contient trois bancs d'essai de courte durée (3 jours simulés, `PeriodLength = 1D`) préparés pour valider le fonctionnement du modèle couplé LMDZ–ORCHIDEE 35 km sur Toubkal avant d'engager des calculs longs.

> [!IMPORTANT]
> **Statut des tests :** Ces tests sont **préparés et prêts pour Slurm**, mais **non encore soumis** ni exécutés. Ils ne doivent pas être présentés comme validés tant qu'un job réel n'a pas tourné avec succès.

## 1. Description des trois tests :
- **[`hist_3days/`](hist_3days/)** : Test historique du **1980-01-01 au 1980-01-03** démarrant depuis les fichiers restart certifiés du 1979-12-31.
- **[`ssp245_3days/`](ssp245_3days/)** : Test scénario SSP2-4.5 du **2015-01-01 au 2015-01-03** démarrant depuis le restart HIST du 2014-12-31.
- **[`ssp585_3days/`](ssp585_3days/)** : Test scénario SSP5-8.5 du **2015-01-01 au 2015-01-03** démarrant depuis le restart HIST du 2014-12-31.

## 2. Procédure de soumission Slurm (lorsqu'autorisé) :
```bash
# Exemple pour hist_3days :
cd configs/lmdzor/r35/tests/hist_3days
sbatch Job_test-3days-HIST-r35
```
