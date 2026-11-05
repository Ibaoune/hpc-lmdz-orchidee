# Test Couplé SSP2-4.5 3 Jours (`ssp245_3days`)

## 1. Objectif
Tester la transition entre la simulation historique (terminée au 2014-12-31) et le scénario climatique futur moyen SSP2-4.5 sur 3 jours nominaux (2015-01-01 au 2015-01-03).

## 2. Configuration (`config.card`)
- **Période :** 2015-01-01 au 2015-01-03 (`PeriodLength = 1D`, 3 pas journaliers).
- **Restarts utilisés (date 2014-12-31, fin de HIST) :**
  - ATM : `sim35km-HIST-after2002_20141231_restart.nc` et `restartphy.nc`.
  - SRF : `sim35km-HIST-after2002_20141231_sechiba_rest.nc`.
  - SBG : `sim35km-HIST-after2002_20141231_stomate_rest.nc`.
- **Forçages atmosphériques SSP2-4.5 :**
  - GES : `CO2_CMIP6_2015_2500.txt`, `CH4`, `N2O`, `CFC` sous `ATM/GHG/CMIP6/ssp/ssp245/`.
  - Aérosols : `aerosols2015_from_inca.nc` sous `ATM/AEROSOLS/CMIP6/Scenarios/ssp245/384x220/L79/`.
  - Conditions aux limites : `ELI-384x220x79-ssp245-T01_2015_limit.nc`.

## 3. Forçages de végétation et distinction avec la simulation historique
- **Forçages nominaux utilisés dans ce test (`COMP/orchidee.card`) :**
  - `PFTmap` : `${R_IN}/SRF/PFTMAPS/CMIP6/ESA-LUH2v2/scenarios/ssp245/PFTmap_${year}.nc` (fichier 2015 présent, ~60 Mo).
  - `woodharvest` : `${R_IN}/SRF/WOODHARVEST/LUH2v2/scenarios/ssp245/woodharvest_${year}.nc` (fichier 2015 présent, ~8 Mo).
- **Différence avec la simulation historique archivée :**
  - La simulation archivée de production (`experiments/ssp245/`, `sim35km-HIST-future-ssp245`) pointait par erreur vers `scenarios/ssp585/PFTmap_${year}.nc` et `scenarios/ssp585/woodharvest_${year}.nc`.
  - Pour préserver la traçabilité historique stricte, `experiments/ssp245/COMP/orchidee.card` est conservé inchangé avec ses chemins d'origine `ssp585`.
  - En revanche, le présent test nominal `tests/ssp245_3days/` utilise les véritables forçages du scénario SSP245.
