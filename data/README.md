# Sources and measurements

Apart from the adapted configuration in `config/template/`, the repository includes no third-party files. `scripts/01_download.py` downloads every source listed in `sources.json` and verifies each file against its SHA-256. All measurements are from Falling Creek Reservoir, Virginia, USA, site 50 (deepest point, 37.30325° N, 79.8373° W).

## Sources

| Source | Version | Identifier | Licence |
|---|---|---|---|
| GLM-AED configuration, executable (glm+ 3.3.1a2, AED 2.1.0) and forcing (Kurucz et al., 2024b) | v1.0.1 | doi:10.5281/zenodo.11183590 | CC BY 4.0 |
| Laboratory filtered chlorophyll-a, file `filt-chla_2014_2024.csv` (Carey et al., 2025) | EDI data package edi.555.5 | doi:10.6073/pasta/7cc4e7541c75738e5c4d6bb948c2df8e | CC BY 4.0 |
| CTD profiles of temperature, oxygen and photosynthetically active radiation, and the CTD maintenance log (Carey et al., 2026a) | EDI data package edi.200.16 | doi:10.6073/pasta/96955da18a4bb1d01137837e73e79d5e | CC BY-NC 4.0 |
| Secchi depth, file `secchi_2013_2025.csv` (Carey et al., 2026b) | EDI data package edi.198.14 | doi:10.6073/pasta/b20ba95d3667532e4dd8bede8a00ab52 | CC BY-NC 4.0 |
| Intel Fortran runtime libraries required by glm+ (`intel-fortran-rt`, `intel-cmplr-lib-rt`) | 2023.1.0 | anaconda.org/anaconda | Intel End User License Agreement for Developer Tools |
| GLM3 source headers `glm.h` and `glm_types.h` | commit 8453b1d | github.com/AquaticEcoDynamics/GLM3 | GPL-3.0 |
| Ubuntu base image | 20.04, digest in `docker/Dockerfile` | Docker Hub | Ubuntu licences |

The measurement files are downloaded from the GitHub repository of the Carey Lab, CareyLabVT/Reservoirs, at commit d12c7dd3fbf0e2610f29ba66b6650fc4021e3913. The publishing files in that repository build the three EDI data packages above. The 2018 CTD profiles are published in edi.200.16 as `ctd_2013_2025.csv`. This analysis reads the individual cast files from which that file is built (repository folder `Data/DataNotYetUploadedToEDI/Raw_CTD/csv_outputs`). The download is about 270 MB, mostly the Zenodo archive.

Citations as recommended by EDI:

- Carey, C.C., A. Breef-Pilz, K.K. Hoffman, D.W. Howard, S. Kenny, E. Tipper, and S. Tannheiser. 2025. Filtered chlorophyll a time series for Beaverdam Reservoir, Carvins Cove Reservoir, Claytor Lake, Falling Creek Reservoir, Gatewood Reservoir, Smith Mountain Lake, Spring Hollow Reservoir in southwestern Virginia, and Lake Sunapee in Sunapee, New Hampshire, USA during 2014-2024 ver 5. Environmental Data Initiative. https://doi.org/10.6073/pasta/7cc4e7541c75738e5c4d6bb948c2df8e (Accessed 2026-10-13).
- Carey, C.C., A. Breef-Pilz, and K.K. Hoffman. 2026a. Time series of high-frequency profiles of depth, temperature, dissolved oxygen, conductivity, specific conductance, chlorophyll a, turbidity, pH, oxidation-reduction potential, photosynthetically active radiation, colored dissolved organic matter, phycocyanin, phycoerythrin, and descent rate for Beaverdam Reservoir, Carvins Cove Reservoir, Falling Creek Reservoir, Gatewood Reservoir, and Spring Hollow Reservoir in southwestern Virginia, USA 2013-2025 ver 16. Environmental Data Initiative. https://doi.org/10.6073/pasta/96955da18a4bb1d01137837e73e79d5e (Accessed 2026-10-13).
- Carey, C.C., A. Breef-Pilz, D.X. Edwards, Q. Campbell, S.E. Kenny, and E. Walters. 2026b. Secchi depth data and discrete depth profiles of water temperature, dissolved oxygen, conductivity, specific conductance, photosynthetic active radiation, oxidation-reduction potential, and pH for Beaverdam Reservoir, Carvins Cove Reservoir, Falling Creek Reservoir, Gatewood Reservoir, and Spring Hollow Reservoir in southwestern Virginia, USA 2013-2025 ver 14. Environmental Data Initiative. https://doi.org/10.6073/pasta/b20ba95d3667532e4dd8bede8a00ab52 (Accessed 2026-10-13).

## Measurements used

All times are converted to model hours: hours since 1 December 2016 00:00 local standard time (UTC−5).

**Assimilated temperature.** 18 CTD casts on 17 dates in 2018: 7, 11, 14, 17, 28, 30 and 31 May, 4, 7, 11, 14, 18 (two casts), 21 and 25 June, and 2, 5 and 9 July. The value is the mean of the 0.1 m bin means over the top 0.5 m. Each value enters at the cast time rounded up to the next full hour, with an observation error of 0.2 °C.

**Assimilated Chla.** Laboratory Chla at 0.1 and 1.6 m (flag 0, replicate mean), on the 7 dates with both depths: 7, 21 and 28 May and 4, 11, 18 and 25 June 2018. Each date gives one value, `y = w1·Chla(0.1 m) + w2·Chla(1.6 m)`. The weights follow `exp(−2 K_d z)` down to `min(ln 100 / (2 K_d), 3 m)`, with the 0.1 m sample for 0 to 0.85 m and the 1.6 m sample below. Laboratory values carry a date only. They enter at the following local midnight, with an observation error of 5% of the value.
- K<sub>d</sub> comes from the CTD profile of photosynthetically active radiation nearest in time (within 3 days, daylight, at least 5 points, R² ≥ 0.9), fitted from the shallowest reading to the depth of 1% of its value or 3 m.
- Where no profile qualifies (7 and 21 May), `K_d = 1.7 / Secchi depth` (Secchi depth 3.14 m and 1.30 m).

**Withheld Chla.** Laboratory Chla at 3.8 and 5 m (flag 0, replicate mean): 11 samples on 28 May, 4, 11, 18 and 25 June and 2 July 2018. Model values come from the state at local midnight of the sampling day.

**Withheld oxygen.** CTD oxygen (mg L⁻¹) from 0 to 9.2 m in 0.2 m bins: 14 casts on 13 dates in 2018, 644 bins. The casts are those of 28, 30 and 31 May, 4, 7, 11, 14, 18 (two casts), 21 and 25 June, and 2, 5 and 9 July. Model values come from the forecast state at the cast time rounded up to the next full hour.

The withheld measurements enter only the evaluation (`scripts/05_evaluate.py`), never the filter.
