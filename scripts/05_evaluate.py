"""Evaluate both configurations against the withheld Chla (3.8 and 5 m) and oxygen (0-9.2 m) measurements.

Model values come from the forecast state (before any update at that hour) of every member:
Chla in the native layer that contains the sample depth, oxygen in the native layer at each CTD
reading, averaged within 0.2 m bins. The continuous ranked probability score (CRPS) is the empirical
ensemble form. Chla scores are averaged over the points of each date, then over dates. Oxygen scores
are averaged over the bins of each cast, the casts of each date, then over dates. The 95% intervals
come from a paired bootstrap over dates (10,000 draws, seed 20260928). Writes results/crps_summary.csv
and output/evaluation/.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from fcr_assimilation.archive import load_archive  # noqa: E402
from fcr_assimilation.analysis import FACTORS  # noqa: E402

OBS = ROOT / 'output/observations'
RUNS = ROOT / 'output/runs'
CONFIGURATIONS = ['temperature', 'temperature_chla']
N, SEED, DRAWS = 128, 20260928, 10000


def predict(configuration, event, depths, variable):
    out = []
    for s in load_archive(RUNS / configuration / 'native' / f'event_{event:03d}.npz', 'forecast'):
        h = s['geometry'][:, 0]
        top, bottom = h[-1] - h, h[-1] - np.r_[0, h[:-1]]
        hit = (np.asarray(depths)[:, None] >= top) & (np.asarray(depths)[:, None] < bottom)
        assert np.all(hit.sum(1) == 1)
        v = s['wq'][s['ids']].T @ FACTORS if variable == 'chla' else s['wq'][3] * .032   # oxygen: mmol m-3 to mg L-1
        out.append(v[np.argmax(hit, axis=1)])
    assert len(out) == N
    return np.array(out)


def chla_support():
    support = pd.read_csv(OBS / 'withheld_chla.csv', float_precision='round_trip')
    predictions = {c: np.concatenate([predict(c, int(e), g.depth_m, 'chla') for e, g in support.groupby('event_id', sort=False)], axis=1)
                   for c in CONFIGURATIONS}
    return support, support['observed'].to_numpy(), predictions


def oxygen_support():
    readings = pd.read_csv(OBS / 'withheld_oxygen.csv', float_precision='round_trip')
    rows, predictions = [], {c: [] for c in CONFIGURATIONS}
    for cast, g in readings.groupby('cast', sort=True):
        event = int(g['event_id'].iloc[0])
        bins = sorted(g['bin'].unique())
        for b in bins:
            part = g[g['bin'].eq(b)]
            rows.append(dict(cast=cast, date=g['date'].iloc[0], bin=int(b), observed=float(part['observed'].mean())))
        for c in CONFIGURATIONS:
            x = predict(c, event, g['depth_m'], 'oxygen')
            predictions[c].append(np.array([x[:, g['bin'].eq(b).to_numpy()].mean(1) for b in bins]).T)
    support = pd.DataFrame(rows)
    return support, support['observed'].to_numpy(), {c: np.concatenate(p, axis=1) for c, p in predictions.items()}


def point_crps(x, y):
    """Empirical CRPS of members x (N, points) for observations y (points)."""
    n = x.shape[0]
    mean_distance = np.abs(x - y).mean(axis=0)
    half_pair_mean = ((2 * np.arange(1, n + 1) - n - 1)[:, None] * np.sort(x, axis=0)).sum(axis=0) / (n * n)
    return mean_distance - half_pair_mean


def by_date(support, values, variable):
    t = support.assign(value=values)
    dates = sorted(support['date'].unique())
    if variable == 'oxygen':
        return t.groupby(['date', 'cast'])['value'].mean().groupby('date').mean().reindex(dates)
    return t.groupby('date')['value'].mean().reindex(dates)


def main():
    (ROOT / 'output/evaluation').mkdir(parents=True, exist_ok=True)
    rows = []
    for variable, label, unit, make in [('chla', 'Chla at 3.8 and 5 m', 'mg m-3', chla_support),
                                        ('oxygen', 'Oxygen over 0-9.2 m', 'mg L-1', oxygen_support)]:
        support, y, predictions = make()
        support.to_csv(ROOT / f'output/evaluation/{variable}_support.csv', index=False, float_format='%.17g')
        np.savez_compressed(ROOT / f'output/evaluation/{variable}_predictions.npz', observed=y, **predictions)
        daily = {c: by_date(support, point_crps(predictions[c], y), variable).to_numpy() for c in CONFIGURATIONS}
        n_dates = support['date'].nunique()
        ix = np.random.default_rng(SEED).integers(0, n_dates, size=(DRAWS, n_dates))
        ref, new = daily['temperature'], daily['temperature_chla']
        gain = (ref - new)[ix].mean(axis=1)
        percent = 100 * gain / ref[ix].mean(axis=1)
        low, high = np.percentile(percent, [2.5, 97.5])
        rows.append(dict(variable=label, unit=unit, points=len(y), dates=n_dates,
                         crps_temperature_updates=ref.mean(), crps_temperature_and_chla_updates=new.mean(),
                         reduction_percent=100 * (ref.mean() - new.mean()) / ref.mean(),
                         reduction_percent_low=low, reduction_percent_high=high))
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / 'results/crps_summary.csv', index=False, float_format='%.6f')
    print(table.to_string(index=False))


if __name__ == '__main__':
    main()
