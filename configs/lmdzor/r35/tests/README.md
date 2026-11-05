# Bancs d'Essai Couplés 3 Jours (`tests/`)

Ce dossier contient trois bancs d'essai de courte durée (3 jours simulés, `PeriodLength = 1D`) préparés pour valider le fonctionnement du modèle couplé LMDZ–ORCHIDEE 35 km sur Toubkal avant d'engager des calculs longs.

> [!IMPORTANT]
> **Statut des tests :** Ces tests sont **pleinement configurés et prêts pour la soumission Slurm**, mais **non encore soumis** conformément aux consignes. Ils ne doivent pas être présentés comme validés scientifiquement tant qu'un job réel n'a pas tourné avec succès.

## 1. Description des trois tests :
- **[`hist_3days/`](hist_3days/)** : Test historique du **1980-01-01 au 1980-01-03** démarrant depuis les fichiers restart certifiés du 1979-12-31 (`sim35km-HIST`).
- **[`ssp245_3days/`](ssp245_3days/)** : Test scénario SSP2-4.5 nominal du **2015-01-01 au 2015-01-03** démarrant depuis le restart HIST du 2014-12-31 (`sim35km-HIST-after2002`) avec forçages réels SSP245.
- **[`ssp585_3days/`](ssp585_3days/)** : Test scénario SSP5-8.5 du **2015-01-01 au 2015-01-03** démarrant depuis le restart HIST du 2014-12-31 (`sim35km-HIST-after2002`).

## 2. Accès aux Restarts via l'Espace de Liaison (Stage)
libIGCM requiert rigoureusement la structure d'arborescence :
`${RestartPath}/${RestartJobName}/${Comp}/Restart/${RestartJobName}_${Date}_${File}`

Pour concilier cette exigence avec la réorganisation des archives historiques sur Lustre (qui ne contenaient pas le dossier intermédiaire `${RestartJobName}`), un espace de liaison hors Git est utilisé :
- **Emplacement :** `modipsl/simulations/lmdzor_restart_stage/` (exclu de Git par `.gitignore`).
- **Liens symboliques :**
  - `sim35km-HIST -> /srv/lustre01/.../ARCHIVE/Simulations/LMDZ/r35/hist`
  - `sim35km-HIST-after2002 -> /srv/lustre01/.../ARCHIVE/Simulations/LMDZ/r35/hist`
- **Avantage :** Aucun fichier d'archive n'est modifié, aucun gros volume NetCDF n'est copié, et libIGCM résout nativement les 4 restarts requis (ATM restart + restartphy, SRF sechiba, SBG stomate).

## 3. Procédure de soumission Slurm (lorsqu'autorisé) :
```bash
# Exemple pour hist_3days :
cd configs/lmdzor/r35/tests/hist_3days
sbatch Job_test-3days-HIST-r35
```
