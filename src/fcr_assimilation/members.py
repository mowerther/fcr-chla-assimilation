"""Write the 128 member configurations from the template member and config/members.csv.

Six multipliers differ between members (Latin hypercube, seed 20260928):
background light extinction (Kw), initial dissolved nitrogen and phosphorus pools,
phytoplankton sinking speed, and phytoplankton growth, respiration and light
saturation (I_K, shared by the three groups). The template already carries its own
nutrient and sinking multipliers (config/template.json). Members rescale these by the
ratio of the member multiplier to the template multiplier.
"""
from pathlib import Path
import csv
import json
import re
import shutil
import numpy as np
from scipy.stats import qmc

CONFIG = Path(__file__).resolve().parents[2] / 'config'
NAMES = ['background_Kw_multiplier', 'initial_dissolved_NP_multiplier', 'phyto_sinking_magnitude_multiplier',
         'growth_multiplier', 'respiration_multiplier', 'light_IK_multiplier']


def sample(n=128, seed=20260928):
    """Latin hypercube sample of the six multipliers: three uniform, three log-uniform from 1/3 to 3."""
    u = qmc.LatinHypercube(d=6, seed=seed).random(n)
    v = np.empty_like(u)
    v[:, :3] = np.array([.75, .5, .5]) + u[:, :3] * np.array([.5, 1, 1])
    v[:, 3:] = np.exp(np.log(1 / 3) + u[:, 3:] * np.log(9))
    return u, v


def _num(s):
    return [float(x.strip().replace('D', 'e').replace('d', 'e')) for x in s.split(',')]


def _value(s, key):
    m = re.search(r'(?mi)^\s*' + re.escape(key) + r'\s*=([^!\n]+)', s)
    assert m, key
    return m[1].strip()


def _replace(s, key, new):
    pattern = r'(?mi)^(\s*' + re.escape(key) + r'\s*=)[^!\n]+'
    s, n = re.subn(pattern, lambda m: m[1] + ' ' + new, s)
    assert n == 1, (key, n)
    return s


def _fmt(x):
    return ', '.join(format(t, '.17g') for t in np.atleast_1d(x))


def read_members():
    return list(csv.DictReader((CONFIG / 'members.csv').open()))


def write_members(target):
    """Write one folder per member (glm3.nml and two AED files) under target."""
    target = Path(target)
    template = json.loads((CONFIG / 'template.json').read_text())
    src = CONFIG / 'template'
    nmlref = (src / 'glm3.nml').read_text()
    aedref = (src / 'aed/aed2_4zones.nml').read_text()
    phyref = (src / 'aed/aed2_phyto_pars_2May2022_RQT.nml').read_text()
    members = read_members()
    u, v = sample()
    for i, r in enumerate(members):
        np.testing.assert_array_equal([float(r[k]) for k in NAMES], v[i])
        np.testing.assert_array_equal([float(r[f'unit_{k + 1}']) for k in range(6)], u[i])
    for r in members:
        m = {k: float(r[k]) for k in NAMES}
        folder = target / r['member_id']
        shutil.copytree(src, folder)
        nml = _replace(nmlref, 'kw', _fmt(template['published_kw'] * m['background_Kw_multiplier']))
        a = np.array(_num(_value(nml, 'wq_init_vals'))).reshape(5, 29)
        a[1:4] *= m['initial_dissolved_NP_multiplier'] / template['initial_dissolved_NP_multiplier']
        nml = _replace(nml, 'wq_init_vals', _fmt(a.ravel()))
        (folder / 'glm3.nml').write_text(nml)
        aed = aedref
        for key in ['don_initial', 'donr_initial', 'dop_initial', 'dopr_initial']:
            aed = _replace(aed, key, _fmt(float(_value(aedref, key)) * m['initial_dissolved_NP_multiplier'] / template['initial_dissolved_NP_multiplier']))
        (folder / 'aed/aed2_4zones.nml').write_text(aed)
        phy = _replace(phyref, 'pd%w_p', _fmt(np.array(_num(_value(phyref, 'pd%w_p'))) * m['phyto_sinking_magnitude_multiplier'] / template['phyto_sinking_magnitude_multiplier']))
        for key, mult in [('pd%r_growth', m['growth_multiplier']), ('pd%r_resp', m['respiration_multiplier']), ('pd%i_k', m['light_IK_multiplier'])]:
            phy = _replace(phy, key, _fmt(np.array(_num(_value(phyref, key))) * mult))
        (folder / 'aed/aed2_phyto_pars_2May2022_RQT.nml').write_text(phy)
        assert all(x < 0 for x in _num(_value(phy, 'pd%w_p')))
    return [r['member_id'] for r in members]
