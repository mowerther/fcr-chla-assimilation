# Model configuration

## Source

The GLM-AED configuration `Deepm2_exm_weight2` for Falling Creek Reservoir from Kurucz et al. (2024b): Zenodo record 10.5281/zenodo.11183590, version 1.0.1, licence CC BY 4.0. The executable (glm+ 3.3.1a2 with AED 2.1.0) and the forcing files come from the same record. `scripts/01_download.py` downloads them.

## Template member (`template/`)

`template/` holds `glm3.nml`, `aed/aed2_4zones.nml` and `aed/aed2_phyto_pars_2May2022_RQT.nml`, adapted from the published files. The namelists are written in a normalised format: keys in lower case, new indentation and numbers at full precision. The values differ from the published files as follows.

`glm3.nml`:
- `stop` is 2018-07-11 (published: 2019-12-31). The start stays at 2016-12-01.
- `out_fn` is `trajectory.partial` and `nsave` is 1 (published: `output`, 24). The model output is not written in these runs. The state exchange stores the states instead.
- `kw` is replaced in every member (see below).
- `bsn_wid_outl`: the malformed separator in `226.03. 226.03` is a comma. Both values are unchanged.
- `the_temps`: the published list has 28 values for 29 depths. The missing value at 9.2 m, 8.005 °C, comes from the published temperature profile of 2 December 2016.
- `wq_init_vals`: the rows of `NIT_amm`, `NIT_nit` and `PHS_frp` are multiplied by 0.6074009986278577, the nutrient multiplier of the template.

`aed2_4zones.nml`:
- `don_initial`, `donr_initial`, `dop_initial` and `dopr_initial` are multiplied by the same 0.6074009986278577.

`aed2_phyto_pars_2May2022_RQT.nml`:
- `pd%w_p` (sinking speed) is multiplied by 0.733813353129992, the sinking multiplier of the template.

`template.json` records these two multipliers and the published `kw` of 0.4136778243201735.

## Members (`members.csv`)

128 members drawn by Latin hypercube sampling (seed 20260928, SciPy `qmc.LatinHypercube`), with six multipliers:

| Multiplier | Range | Applied to |
|---|---|---|
| Background light extinction | 0.75 to 1.25, uniform | published `kw` |
| Initial dissolved nitrogen and phosphorus | 0.5 to 1.5, uniform | `NIT_amm`, `NIT_nit`, `PHS_frp`, DON, DOP (refractory pools included) |
| Phytoplankton sinking speed | 0.5 to 1.5, uniform | `pd%w_p` of the three groups |
| Phytoplankton growth rate | 1/3 to 3, log-uniform | `pd%r_growth` |
| Phytoplankton respiration rate | 1/3 to 3, log-uniform | `pd%r_resp` |
| Light saturation | 1/3 to 3, log-uniform | `pd%i_k` |

The growth, respiration and light saturation multipliers are shared by cyanobacteria, green algae and diatoms. The nutrient and sinking values of the template are rescaled by the ratio of the member multiplier to the template multiplier. Forcing, initial temperatures and all other parameters are identical in every member. The filter does not update any parameter. `src/fcr_assimilation/members.py` writes the members and regenerates the Latin hypercube to confirm `members.csv`.

## Event calendar (`events.csv`)

83 model hours at which every member pauses and its state is stored: the 65 local midnights (04:00 UTC) from 7 May to 10 July 2018 and the 18 hours of the temperature casts (each cast time rounded up to the next full hour). Model hour 0 is 1 December 2016 00:00 local standard time (UTC−5). `scripts/02_prepare_observations.py` confirms the calendar from the measurements.
