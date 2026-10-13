# Assimilating chlorophyll-a near the surface into GLM-AED for Falling Creek Reservoir

Preliminary analysis by Mortimer Werther (Eawag, ORCID 0000-0002-0775-9285).

## Summary

This repository assimilates chlorophyll-a (Chla) near the surface, measured in situ and combined as
a satellite sensor observes the surface layer, into an ensemble of 128 members of the General Lake
Model coupled to the Aquatic EcoDynamics library (GLM-AED) for Falling Creek Reservoir, Virginia,
USA. Compared with surface temperature updates alone, assimilating Chla lowered the continuous
ranked probability score (CRPS) of simulated Chla at 3.8 and 5 m, below the assimilated depths
(0.1 and 1.6 m), by 37% (from 1.80 to 1.14 mg m⁻³, 95% interval 21 to 51%). Information from the surface thus reached
the interior through the process model.

## Configuration

- Model: the published GLM-AED configuration for Falling Creek Reservoir (`Deepm2_exm_weight2`),
  calibrated against temperature and oxygen profiles from December 2016 to December 2019, which
  include the 2018 profiles used here (Kurucz et al., 2024a, 2024b).
- Ensemble: 128 members that differ in six parameters, drawn by Latin hypercube sampling: background
  light extinction (×0.75 to 1.25), initial dissolved nitrogen and phosphorus (×0.5 to 1.5),
  phytoplankton sinking speed (×0.5 to 1.5), and phytoplankton growth, respiration and light
  saturation (×1/3 to 3, log-uniform). Forcing and initial temperatures are identical. Every member runs from 1 December 2016.
- Period: 7 May to 10 July 2018. Temperature updates on 17 dates (18 casts). Chla updates on 7 dates
  from 7 May to 25 June, about weekly.
- Filter: a localised ensemble transform Kalman filter (Bishop et al., 2001; Hunt et al., 2007),
  applied to one observation at a time, with Gaspari–Cohn localisation over a half-width of 1 m
  (Gaspari and Cohn, 1999) and relaxation to prior spread with a factor of 0.5 (Whitaker and Hamill,
  2012). Where an update would take Chla below the model minimum, its increments are scaled down by
  one common factor. The filter updates temperature and phytoplankton. Each member keeps its own
  parameter values throughout the run.
- Temperature observations: surface temperature measured by CTD, the mean of the 0.1 m bins over the
  top 0.5 m, with an observation error of 0.2 °C.
- Chla observations: Chla measured at 0.1 and 1.6 m, combined into one value per date and weighted
  by the measured light attenuation, as a satellite sensor observes the surface layer
  (Gordon and McCluney, 1975). Observation error: 5% of the value.

## Evaluation

Two configurations share the model, forcing and ensemble: surface temperature updates alone, and
surface temperature updates with Chla. Both are evaluated against Chla measured at 3.8 and 5 m
(11 samples on 6 dates), below the assimilated depths (0.1 and 1.6 m), and against oxygen over the
profile (0 to 9.2 m, 14 CTD casts on 13 dates), with the CRPS. The deepest model layer changed by any update
ends at 3.78 m. The 95% intervals come from a paired bootstrap over dates (10,000 draws).

| Variable | Temperature updates | Temperature and Chla updates | Change (95% interval) |
|---|---|---|---|
| Chla at 3.8 and 5 m (CRPS, mg m⁻³) | 1.80 | 1.14 | −37% (21 to 51%) |
| Oxygen over the profile (CRPS, mg L⁻¹) | 0.939 | 0.925 | −1.5% (0.2 to 2.9%) |

The filter updated phytoplankton but neither oxygen nor the rates of phytoplankton growth and
sinking and of organic matter decomposition. An oxygenation system operated in the bottom water
throughout this period. In the configuration it supplies about 7.4 kg O₂ per day (15 kg per day in
the inflow file, scaled by the calibrated inflow factor of 0.49). The oxygen error over the profile
therefore changed by less than 2%. Constraining these rates with a surface record spanning several
years is the subject of ongoing work.

## Reproducing the results

Requirements: conda (or Python 3.12), Docker with linux/amd64 images (emulated on Apple silicon),
about 2 GB of free disk space and an internet connection.

```sh
conda env create -f environment.yml
conda activate fcr-chla-assimilation
python scripts/01_download.py              # model, runtime libraries and measurements, about 270 MB
python scripts/02_prepare_observations.py  # assimilated and withheld measurements
python scripts/03_build_runtime.py         # Docker image, state exchange library, 128 members
python scripts/04_assimilate.py            # both configurations, about 10 minutes each
python scripts/05_evaluate.py              # writes results/crps_summary.csv
```

Each configuration runs every member from 1 December 2016 to 10 July 2018 with two members
integrating at a time. All members pause at 83 event hours, where the filter updates them.
`config/README.md` describes the configuration and the ensemble. `data/README.md` lists the
sources and the selection of measurements.

## Sources and licences

- GLM-AED configuration, executable and forcing: Kurucz et al. (2024b), version 1.0.1,
  doi:10.5281/zenodo.11183590, CC BY 4.0. Downloaded. The adapted template member is in
  `config/template/` (CC BY 4.0, changes listed in `config/README.md`).
- Laboratory chlorophyll-a: EDI data package edi.555.5 (Carey et al., 2025), CC BY 4.0.
- CTD profiles and CTD maintenance log: EDI data package edi.200.16 (Carey et al., 2026a),
  CC BY-NC 4.0.
- Secchi depth: EDI data package edi.198.14 (Carey et al., 2026b), CC BY-NC 4.0.
- The three EDI packages are built from the GitHub repository CareyLabVT/Reservoirs. The scripts
  download the files from that repository at commit d12c7dd.
- Intel Fortran runtime libraries 2023.1.0 (Anaconda), Intel End User License Agreement. Downloaded.
- GLM3 headers, commit 8453b1d, GPL-3.0. Downloaded.

Code: MIT licence, except `native/native_exchange.c` (GPL-3.0-or-later). Derived results and
figures: CC BY-NC 4.0, because two input packages (edi.200.16, edi.198.14) are CC BY-NC 4.0.
Details are in `LICENSE-data.md`.

## How to cite

Werther, M., 2026. Assimilating chlorophyll-a near the surface into GLM-AED for Falling Creek
Reservoir (v1.0.0). GitHub. https://github.com/mowerther/fcr-chla-assimilation/releases/tag/v1.0.0

## References

Bishop, C.H., Etherton, B.J., Majumdar, S.J., 2001. Adaptive sampling with the ensemble transform
Kalman filter. Part I: theoretical aspects. Monthly Weather Review 129, 420–436.
https://doi.org/10.1175/1520-0493(2001)129<0420:ASWTET>2.0.CO;2

Carey, C.C., Breef-Pilz, A., Hoffman, K.K., Howard, D.W., Kenny, S., Tipper, E., Tannheiser, S.,
2025. Filtered chlorophyll a time series for Beaverdam Reservoir, Carvins Cove Reservoir, Claytor
Lake, Falling Creek Reservoir, Gatewood Reservoir, Smith Mountain Lake, Spring Hollow Reservoir in
southwestern Virginia, and Lake Sunapee in Sunapee, New Hampshire, USA during 2014-2024. Version 5.
Environmental Data Initiative. https://doi.org/10.6073/pasta/7cc4e7541c75738e5c4d6bb948c2df8e

Carey, C.C., Breef-Pilz, A., Hoffman, K.K., 2026a. Time series of high-frequency profiles of depth,
temperature, dissolved oxygen, conductivity, specific conductance, chlorophyll a, turbidity, pH,
oxidation-reduction potential, photosynthetically active radiation, colored dissolved organic
matter, phycocyanin, phycoerythrin, and descent rate for Beaverdam Reservoir, Carvins Cove
Reservoir, Falling Creek Reservoir, Gatewood Reservoir, and Spring Hollow Reservoir in southwestern
Virginia, USA 2013-2025. Version 16. Environmental Data Initiative.
https://doi.org/10.6073/pasta/96955da18a4bb1d01137837e73e79d5e

Carey, C.C., Breef-Pilz, A., Edwards, D.X., Campbell, Q., Kenny, S.E., Walters, E., 2026b. Secchi
depth data and discrete depth profiles of water temperature, dissolved oxygen, conductivity,
specific conductance, photosynthetic active radiation, oxidation-reduction potential, and pH for
Beaverdam Reservoir, Carvins Cove Reservoir, Falling Creek Reservoir, Gatewood Reservoir, and Spring
Hollow Reservoir in southwestern Virginia, USA 2013-2025. Version 14. Environmental Data Initiative.
https://doi.org/10.6073/pasta/b20ba95d3667532e4dd8bede8a00ab52

Gaspari, G., Cohn, S.E., 1999. Construction of correlation functions in two and three dimensions.
Quarterly Journal of the Royal Meteorological Society 125, 723–757. https://doi.org/10.1002/qj.49712555417

Gordon, H.R., McCluney, W.R., 1975. Estimation of the depth of sunlight penetration in the sea for
remote sensing. Applied Optics 14, 413–416. https://doi.org/10.1364/AO.14.000413

Hunt, B.R., Kostelich, E.J., Szunyogh, I., 2007. Efficient data assimilation for spatiotemporal
chaos: a local ensemble transform Kalman filter. Physica D 230, 112–126.
https://doi.org/10.1016/j.physd.2006.11.008

Kurucz, K., Carey, C.C., Huang, P., de Sousa, E.R., White, J.T., Hipsey, M.R., 2024a. Using
system-inspired metrics to improve water quality prediction in stratified lakes. Water Resources
Research 60, e2023WR036570. https://doi.org/10.1029/2023WR036570

Kurucz, K., Hipsey, M., Carey, C., Huang, P., de Sousa, E., White, J., 2024b. Data and software for
the GLM-AED simulation and model calibration of the Falling Creek Reservoir (v1.0.1). Zenodo.
https://doi.org/10.5281/zenodo.11183590

Whitaker, J.S., Hamill, T.M., 2012. Evaluating methods to account for system errors in ensemble
data assimilation. Monthly Weather Review 140, 3078–3089. https://doi.org/10.1175/MWR-D-11-00276.1

## Acknowledgements

Falling Creek Reservoir measurements come from the Carey Lab at Virginia Tech (Carey et al., 2025,
2026a, 2026b). The Western Virginia Water Authority owns and operates the reservoir. The GLM-AED
configuration comes from Kurucz et al. (2024b).
