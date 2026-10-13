"""Lossless storage of the forecast and accepted native states at each event."""
import numpy as np


def archive(path, forecast, accepted):
    counts = np.array([s['n'] for s in forecast], dtype=np.int16)
    assert all(s['n'] == a['n'] for s, a in zip(forecast, accepted))
    arrays = dict(counts=counts, ids=forecast[0]['ids'],
                  geometry_forecast=np.concatenate([s['geometry'] for s in forecast]),
                  accepted_temperature_density=np.concatenate([s['geometry'][:, [2, 4]] for s in accepted]),
                  wq_forecast=np.concatenate([s['wq'].T for s in forecast]),
                  accepted_phytoplankton=np.concatenate([s['wq'][s['ids']].T for s in accepted]))
    np.savez(path, **arrays)
    # Read the file back before the exchange files are removed.
    with np.load(path) as z:
        for k, v in arrays.items():
            np.testing.assert_array_equal(z[k], v)
    recovered = load_archive(path)
    for a, b in zip(accepted, recovered):
        np.testing.assert_array_equal(a['geometry'], b['geometry'])
        np.testing.assert_array_equal(a['wq'], b['wq'])


def compress(path):
    with np.load(path) as z:
        arrays = {k: z[k] for k in z.files}
    tmp = path.with_suffix('.tmp')
    with tmp.open('wb') as f:
        np.savez_compressed(f, **arrays)
    tmp.replace(path)


def load_archive(path, phase='accepted'):
    with np.load(path) as source:
        z = {k: source[k] for k in source.files}
    offsets = np.r_[0, np.cumsum(z['counts'])]
    ids = z['ids']
    out = []
    for i, n in enumerate(z['counts']):
        g = z['geometry_forecast'][offsets[i]:offsets[i + 1]].copy()
        w = z['wq_forecast'][offsets[i]:offsets[i + 1]].T.copy()
        if phase == 'accepted':
            g[:, [2, 4]] = z['accepted_temperature_density'][offsets[i]:offsets[i + 1]]
            w[ids] = z['accepted_phytoplankton'][offsets[i]:offsets[i + 1]].T
        out.append(dict(n=int(n), ids=ids, geometry=g, wq=w, values=np.column_stack([g[:, 2], w[ids].T])))
    return out
