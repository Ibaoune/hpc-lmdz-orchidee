# Configurations des Simulations LMDZ pour Toubkal (35 km)

Ce répertoire `Toubkal_Configs/` sauvegarde les fichiers de configuration de référence que tu as utilisés pour lancer tes simulations climatiques historiques et futures (SSP245, SSP585) sur le cluster Toubkal. 

Ces configurations te permettent de reproduire exactement le même environnement de simulation de manière propre et versionnée.

## Structure

- `ETAT0/` : Configuration pour générer les états initiaux (`start.nc`, `startphy.nc`) et les fichiers de forçage (`limit.nc`) depuis les données globales (`ExpType=CREATE_amip`).
- `HIST/` : Configuration de la simulation historique principale de 1979 à 2014 (`ExpType=LMDZOR/amip`).
- `SSP245/` : Configuration pour le scénario SSP245 de 2015 à 2100.
- `SSP585/` : Configuration pour le scénario SSP585 de 2015 à 2100.

## Comment relancer une simulation ?

1. Assure-toi que ton modèle est compilé dans `modipsl/config/LMDZOR_v6/`.
2. Crée un nouveau dossier de simulation en t'inspirant du template : 
   ```bash
   cp -r Toubkal_Configs/HIST modipsl/config/LMDZOR_v6/sim35km-HIST-new
   ```
3. Adapte les paramètres de date, de Restart et les chemins d'accès au besoin dans `config.card`.
4. Lance la simulation via la commande habituelle (ex. `ins_job`).

## Gestion du stockage (Rappel)

**Attention** : Les sorties générées par ces simulations (mensuelles, journalières, HF) remplissent très vite l'espace scratch de libIGCM (`IGCM_OUT`). Il est fortement recommandé d'archiver régulièrement ces données vers ton lustre de grande capacité :
`~/lustre/.../ARCHIVE/Simulations/LMDZ/r35/`
