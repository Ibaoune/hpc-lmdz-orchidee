# Test Couplé SSP5-8.5 3 Jours (`ssp585_3days`)

## 1. Objectif
Tester la transition entre la simulation historique (2014-12-31) et le scénario climatique futur extrême SSP5-8.5 sur 3 jours (2015-01-01 au 2015-01-03).

## 2. Configuration (`config.card`)
- **Période :** 2015-01-01 au 2015-01-03 (`PeriodLength = 1D`).
- **Restarts utilisés (date 2014-12-31, fin de HIST) :**
  - ATM, SRF et SBG identiques au test SSP245 (même point de départ climatique 2014-12-31).
- **Forçages SSP5-8.5 :**
  - GES : forçages extrêmes sous `ATM/GHG/CMIP6/ssp/ssp585/`.
  - Aérosols : `aerosols2015_from_inca.nc` sous `ATM/AEROSOLS/CMIP6/Scenarios/ssp585/384x220/L79/`.
  - Conditions aux limites : `ELI-384x220x79-ssp585-T01_2015_limit.nc`.
  - Végétation : `PFTmap_2015.nc` sous `SRF/PFTMAPS/CMIP6/ESA-LUH2v2/scenarios/ssp585/`.
