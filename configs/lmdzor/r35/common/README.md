# Paramètres communs LMDZ–ORCHIDEE 35 km (`common/`)

Ce dossier rassemble les fichiers de configuration physique, dynamique et diagnostique canoniques partagés par l'ensemble des expériences couplées régionales sur le domaine Maroc (résolution 35 km, grille $384 \times 220 \times 79$).

## 1. Fichiers de paramètres (`PARAM/`)
- **`gcm.def_384x220_NPv6.2`** : Paramètres de la grille atmosphérique étirée et zoomée :
  - Centre du zoom : longitude $-6.5^\circ$, latitude $24.0^\circ$.
  - Facteur de grossissement : $grossismx = grossismy = 2.8$ (résolution de maille fine d'environ 35 km sur le Maroc et l'Afrique du Nord-Ouest).
  - Pas de temps dynamique : `day_step = 5760` (intégration leapfrog-Matsuno).
  - Pas de temps physique : `iphysiq = 60`.
- **`physiq.def_NPv6.2`** : Nouvelle physique LMDZ6 (schémas thermiques, convection, rayonnement, microphysique des nuages).
- **`vert.def_L79`** : Discrétisation verticale sur 79 niveaux hybrides $\sigma-p$.
- **`orchidee.def_CWRR`** : Paramétrage du modèle de surface continentale ORCHIDEE avec schéma hydrologique sol à 11 couches (CWRR).
- **`iodef.xml`, `ping_LMDZ.xml`, `ping_orchidee.xml`** : Définition des contextes et flux d'E/S parallèles XIOS 2.
- **`run.def`, `guide.def`, `traceur.def`** : Paramètres généraux et traceurs d'advection.

## 2. Post-traitement et monitoring (`POST/`)
- `monitoring01_lmdz.cfg`, `monitoring01_orchidee.cfg`, `monitoring01_stomate.cfg` : configurations des séries temporelles et diagnostics automatisés par libIGCM.

## 3. Provenance
Issus de la configuration validée historiquement sur Toubkal :
`models/modipsl/config/LMDZOR_v6/Sim35km-Maroc-T01/PARAM/`.
