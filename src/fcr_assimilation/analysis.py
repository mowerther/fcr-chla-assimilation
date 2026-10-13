"""Localised ensemble transform Kalman filter update for one observation.

The update acts on physical temperature and total chlorophyll-a (Chla) on a common
200-cell grid. Chla observations update total Chla and keep the share of each
phytoplankton group within each native layer. Temperature observations update
temperature only. Observations are applied one at a time (serial update), with
Gaspari-Cohn localisation of the observation error and bounded relaxation to
prior spread (RTPS) for Chla.
"""
from dataclasses import dataclass
import numpy as np

# Chla per unit phytoplankton carbon (12/Xcc) for cyanobacteria, green algae and diatoms.
FACTORS = np.array([.15, .4, .2])
PHY_FLOOR = .03
CHL_FLOOR = float(FACTORS.sum() * PHY_FLOOR)


class AnalysisError(RuntimeError):
    pass


@dataclass(frozen=True)
class Configuration:
    chl_depth_m: float = .1
    temperature_layer_m: float = .5
    localization_halfwidth_m: float = 1.
    grid_cells: int = 200
    observation_consistency_tolerance_sigma: float = 1e-8
    variant: str = 'bounded_rtps'
    rtps_alpha: float = .5
    rtps_max_sd_multiplier: float = 1.5


def relax_spread(prior, posterior, alpha=.5, max_multiplier=1.5):
    """Bounded RTPS: keep the analysis mean and relax the spread towards the forecast spread.

    Factor per row = min(1 + alpha*(sf/sa - 1), cap), at least 1. A common row factor also
    keeps total Chla at or above its native minimum. No member is clipped.
    """
    if not 0 <= alpha <= 1 or max_multiplier < 1:
        raise ValueError('Invalid RTPS settings.')
    if alpha == 0:
        return np.array(posterior, copy=True)
    prior = np.asarray(prior)
    posterior = np.asarray(posterior)
    sf = prior.std(axis=-1, ddof=1)
    sa = posterior.std(axis=-1, ddof=1)
    ratio = np.divide(sf, sa, out=np.ones_like(sf), where=sa > 1e-12)
    scale = np.minimum(max_multiplier, np.maximum(1, 1 + alpha * (ratio - 1)))
    mean = posterior.mean(axis=-1, keepdims=True)
    anomalies = posterior - mean
    low = anomalies.min(axis=-1)
    feasible = np.divide(mean[..., 0] - CHL_FLOOR, -low, out=np.full_like(low, np.inf), where=low < 0)
    if np.any(posterior < CHL_FLOOR - 1e-10):
        raise AnalysisError('RTPS input is physically infeasible.')
    # Round inward so that a floating-point boundary cannot give negative biomass.
    scale = np.minimum(scale, np.maximum(1, feasible * (1 - 1e-12)))
    return mean + scale[..., None] * anomalies


def geometry(s):
    h = s['geometry'][:, 0]
    return h[-1] - h, h[-1] - np.r_[0, h[:-1]]


def chlorophyll(v):
    return v[:, 1:] @ FACTORS


def expected(s, config=Configuration()):
    top, bottom = geometry(s)
    hit = np.flatnonzero((top <= config.chl_depth_m) & (bottom > config.chl_depth_m))
    if len(hit) != 1:
        raise AnalysisError('Sample depth is outside the native column.')
    length = np.maximum(0, np.minimum(bottom, config.temperature_layer_m) - np.maximum(top, 0))
    if abs(length.sum() - config.temperature_layer_m) > 1e-8:
        raise AnalysisError('Temperature layer is uncovered.')
    return np.array([length @ s['values'][:, 0] / length.sum(), chlorophyll(s['values'])[hit[0]]])


def gaspari_cohn(distance, halfwidth):
    if halfwidth <= 0:
        raise ValueError('Localization halfwidth must be positive.')
    x = np.abs(np.asarray(distance, dtype=float)) / halfwidth
    v = np.zeros_like(x)
    a = x <= 1
    b = (x > 1) & (x < 2)
    r = x[a]
    v[a] = 1 - 5 * r * r / 3 + 5 * r ** 3 / 8 + r ** 4 / 2 - r ** 5 / 4
    r = x[b]
    v[b] = 4 - 5 * r + 5 * r * r / 3 + 5 * r ** 3 / 8 - r ** 4 / 2 + r ** 5 / 12 - 2 / (3 * r)
    return np.maximum(0, np.minimum(1, v))


def scalar_transform(x, h, y, sigma, rho=1.):
    """Symmetric rank-one ETKF for one observation. rho implements a local R = R/rho.

    Rows of x are state entries; the last axis holds the members.
    """
    x = np.asarray(x, float)
    h = np.asarray(h, float)
    rho = np.asarray(rho, float)
    if sigma <= 0 or not np.isfinite(sigma) or not np.isfinite(y):
        raise AnalysisError('Invalid observation/error.')
    n = h.size
    b = h - h.mean()
    norm = float(b @ b)
    if norm <= 1e-24:
        return x.copy()
    rho = np.broadcast_to(rho, x.shape[:-1])
    fraction = np.zeros_like(rho)
    active = rho > 0
    fraction[active] = norm / (norm + (n - 1) * sigma * sigma / rho[active])
    projection = ((x - x.mean(-1, keepdims=True)) @ b) / norm
    direction = fraction[..., None] * (y - h.mean()) + (np.sqrt(1 - fraction)[..., None] - 1) * b
    return x + projection[..., None] * direction


def common_physical(states, config):
    values = []
    maps = []
    centres = []
    for s in states:
        top, bottom = geometry(s)
        h = s['geometry'][-1, 0]
        edges = np.linspace(0, h, config.grid_cells + 1)
        length = np.maximum(0, np.minimum(bottom[None, :], edges[1:, None]) - np.maximum(top[None, :], edges[:-1, None]))
        volume = length / (bottom - top)[None, :] * s['geometry'][:, 1][None, :]
        if (volume.sum(1) <= 0).any():
            raise AnalysisError('Empty analysis cell.')
        field = np.column_stack([s['values'][:, 0], chlorophyll(s['values'])])
        values.append(volume @ field / volume.sum(1)[:, None])
        maps.append(volume)
        centres.append((edges[:-1] + edges[1:]) / 2)
    return np.stack(values, axis=-1), maps, np.mean(centres, axis=0)


def allocate_total(prior, target):
    """Distribute a total Chla target over the groups, keeping each group above its native minimum."""
    if np.any(target < CHL_FLOOR - 1e-10):
        raise AnalysisError('Total chlorophyll proposal below the native minimum.')
    excess = prior[:, 1:] - PHY_FLOOR
    if np.any(excess < -1e-9):
        raise AnalysisError('Prior phytoplankton below its declared native minimum.')
    total = excess @ FACTORS
    out = prior.copy()
    changed = target != chlorophyll(prior)
    active = (total > 1e-12) & changed
    out[active, 1:] = PHY_FLOOR + excess[active] * ((target[active] - CHL_FLOOR) / total[active])[:, None]
    empty = changed & ~active
    out[empty, 1:] = (target[empty] / FACTORS.sum())[:, None]
    return out


def _report(variable, observation, y, sigma, h, realized, linear, before_projection, projection_max, config, states, outputs):
    report = dict(variable=variable, observation_id=observation.get('observation_id'), value=y, sigma=sigma,
                  prior_mean=float(h.mean()), analysis_mean=float(realized.mean()), linear_analysis_mean=float(linear.mean()),
                  operator_max_abs_error=float(np.max(np.abs(realized - linear))),
                  localization_halfwidth_m=config.localization_halfwidth_m,
                  remap_operator_error_before_projection=float(np.max(np.abs(before_projection - linear))),
                  native_operator_projection_max_change=projection_max,
                  prior_sd=float(h.std(ddof=1)), posterior_sd=float(realized.std(ddof=1)), variant=config.variant)
    report['mean_chlorophyll_inventory_increment_kg'] = float(np.mean([
        (chlorophyll(a['values']) - chlorophyll(b['values'])) @ b['geometry'][:, 1] / 1e6
        for b, a in zip(states, outputs)]))
    return report


def update(states, observation, config=Configuration()):
    """Apply one measured temperature or point Chla observation to all members."""
    if observation['variable'] == 'Chla_weighted':
        return update_weighted(states, observation, config)
    variable = observation['variable']
    column = 0 if variable == 'T' else 1
    if variable not in ['T', 'Chla']:
        raise AnalysisError('Unsupported observation variable.')
    y = float(observation['value'])
    sigma = float(observation['sigma'])
    if variable == 'Chla' and y <= 0:
        raise AnalysisError('Nonpositive source value is ineligible.')
    x, maps, depth = common_physical(states, config)
    h = np.array([expected(s, config)[column] for s in states])
    # Distance from the support of the measurement: 0-0.5 m for temperature, the sample depth for Chla.
    distance = np.maximum(0, depth - config.temperature_layer_m) if variable == 'T' else np.abs(depth - config.chl_depth_m)
    rho = gaspari_cohn(distance, config.localization_halfwidth_m)
    if config.variant not in ['baseline', 'bounded_rtps']:
        raise AnalysisError('Unregistered variant.')
    # Unit weight within all native cells that support the observation operator.
    if variable == 'T':
        rho[depth <= config.temperature_layer_m] = 1
    else:
        bottoms = []
        for s in states:
            top, bottom = geometry(s)
            hit = np.flatnonzero((top <= config.chl_depth_m) & (bottom > config.chl_depth_m))[0]
            bottoms.append(bottom[hit])
        if config.chl_depth_m == .1:
            rho[depth <= max(bottoms)] = 1
        else:
            support = np.zeros(len(depth), dtype=bool)
            for s in states:
                top, bottom = geometry(s)
                hit = np.flatnonzero((top <= config.chl_depth_m) & (bottom > config.chl_depth_m))[0]
                support |= (depth >= top[hit]) & (depth <= bottom[hit])
            rho[support] = 1
    proposed = scalar_transform(x[:, column, :], h, y, sigma, rho)
    if variable == 'Chla' and config.variant == 'bounded_rtps':
        relaxed = relax_spread(x[:, column, :], proposed, config.rtps_alpha, config.rtps_max_sd_multiplier)
        proposed = np.where(rho[:, None] > 0, relaxed, proposed)
    delta = proposed - x[:, column, :]
    outputs = []
    for j, s in enumerate(states):
        weights = maps[j].T / maps[j].sum(0)[:, None]
        increment = weights @ delta[:, j]
        v = s['values'].copy()
        if variable == 'T':
            v[:, 0] += increment
            if (v[:, 0] < -2).any() or (v[:, 0] > 45).any():
                raise AnalysisError('Temperature proposal outside physical bounds.')
        else:
            v = allocate_total(v, chlorophyll(v) + increment)
        if not np.isfinite(v).all() or (v[:, 1:] > 10000).any():
            raise AnalysisError('Nonfinite or excessive native state.')
        outputs.append(dict(s, values=v))
    linear = scalar_transform(h, h, y, sigma)
    if variable == 'Chla' and config.variant == 'bounded_rtps':
        linear = relax_spread(h, linear, config.rtps_alpha, config.rtps_max_sd_multiplier)
    before_projection = np.array([expected(s, config)[column] for s in outputs])
    # Member geometries differ, so conservative remapping need not preserve the observation
    # operator. The analysed operator ensemble is restored by the minimum-norm correction
    # confined to the observed native cells. The target is the Kalman analysis.
    projection_max = 0.
    for j, s in enumerate(outputs):
        residual_j = linear[j] - before_projection[j]
        top, bottom = geometry(s)
        v = s['values'].copy()
        if variable == 'T':
            weights = np.maximum(0, np.minimum(bottom, config.temperature_layer_m) - np.maximum(top, 0)) / config.temperature_layer_m
            correction = weights * residual_j / (weights @ weights)
            v[:, 0] += correction
            projection_max = max(projection_max, float(np.max(np.abs(correction))))
            if (v[:, 0] < -2).any() or (v[:, 0] > 45).any():
                raise AnalysisError('Operator projection violates temperature bounds.')
        else:
            hit = np.flatnonzero((top <= config.chl_depth_m) & (bottom > config.chl_depth_m))[0]
            target = chlorophyll(v)
            target[hit] = linear[j]
            v = allocate_total(v, target)
            projection_max = max(projection_max, abs(float(residual_j)))
        outputs[j] = dict(s, values=v)
        if not np.isfinite(v).all() or (v[:, 1:] > 10000).any():
            raise AnalysisError('Operator projection produced a nonfinite or excessive state.')
    realized = np.array([expected(s, config)[column] for s in outputs])
    residual = float(np.max(np.abs(realized - linear)))
    if residual > config.observation_consistency_tolerance_sigma * sigma:
        raise AnalysisError(f'Native observation operator mismatch: {residual:g} > {config.observation_consistency_tolerance_sigma * sigma:g}')
    return outputs, _report(variable, observation, y, sigma, h, realized, linear, before_projection, projection_max, config, states, outputs)


# Light-weighted Chla: one value per date, H(x) = w1 Chla(0.1 m) + w2 Chla(1.6 m).

def weighted_cells(s, observation):
    """Native cells and weights of the operator for one member. Weights of two depths in one cell add."""
    top, bottom = geometry(s)
    cells = {}
    for z, w in ((float(observation['depth_1_m']), float(observation['w1'])), (float(observation['depth_2_m']), float(observation['w2']))):
        if w == 0:
            continue
        hit = np.flatnonzero((top <= z) & (bottom > z))
        if len(hit) != 1:
            raise AnalysisError('Sample depth is outside the native column.')
        cells[int(hit[0])] = cells.get(int(hit[0]), 0.) + w
    return cells


def weighted_expected(s, observation):
    c = chlorophyll(s['values'])
    return sum(w * c[i] for i, w in weighted_cells(s, observation).items())


def update_weighted(states, observation, config=Configuration()):
    """Apply one light-weighted Chla value.

    Localisation: Gaspari-Cohn on the distance to the nearer weighted depth, with unit weight on
    the native operator cells of the weighted depths. The ETKF, bounded RTPS and remapping are
    those of the point Chla update. The analysed operator ensemble is restored in the weighted
    native cells: directly for one cell, by the minimum-norm correction for two.
    """
    if config.variant not in ['baseline', 'bounded_rtps']:
        raise AnalysisError('Unregistered variant.')
    y = float(observation['value'])
    sigma = float(observation['sigma'])
    if y <= 0:
        raise AnalysisError('Nonpositive source value is ineligible.')
    w1, w2 = float(observation['w1']), float(observation['w2'])
    z1, z2 = float(observation['depth_1_m']), float(observation['depth_2_m'])
    if not (w1 > 0 and w2 >= 0 and abs(w1 + w2 - 1) <= 1e-12 and 0 < z1 < z2 <= 3):
        raise AnalysisError('Invalid light-weighted operator.')
    column = 1
    x, maps, depth = common_physical(states, config)
    h = np.array([weighted_expected(s, observation) for s in states])
    weighted = [z for z, w in ((z1, w1), (z2, w2)) if w > 0]
    distance = np.min([np.abs(depth - z) for z in weighted], axis=0)
    rho = gaspari_cohn(distance, config.localization_halfwidth_m)
    for z in weighted:
        if z == .1:
            bottoms = []
            for s in states:
                top, bottom = geometry(s)
                hit = np.flatnonzero((top <= z) & (bottom > z))[0]
                bottoms.append(bottom[hit])
            rho[depth <= max(bottoms)] = 1
        else:
            support = np.zeros(len(depth), dtype=bool)
            for s in states:
                top, bottom = geometry(s)
                hit = np.flatnonzero((top <= z) & (bottom > z))[0]
                support |= (depth >= top[hit]) & (depth <= bottom[hit])
            rho[support] = 1
    proposed = scalar_transform(x[:, column, :], h, y, sigma, rho)
    if config.variant == 'bounded_rtps':
        relaxed = relax_spread(x[:, column, :], proposed, config.rtps_alpha, config.rtps_max_sd_multiplier)
        proposed = np.where(rho[:, None] > 0, relaxed, proposed)
    delta = proposed - x[:, column, :]
    outputs = []
    for j, s in enumerate(states):
        weights = maps[j].T / maps[j].sum(0)[:, None]
        increment = weights @ delta[:, j]
        # Chla of the C-ordered copy, as in the point update (BLAS rounding depends on memory layout).
        v = s['values'].copy()
        v = allocate_total(v, chlorophyll(v) + increment)
        if not np.isfinite(v).all() or (v[:, 1:] > 10000).any():
            raise AnalysisError('Nonfinite or excessive native state.')
        outputs.append(dict(s, values=v))
    linear = scalar_transform(h, h, y, sigma)
    if config.variant == 'bounded_rtps':
        linear = relax_spread(h, linear, config.rtps_alpha, config.rtps_max_sd_multiplier)
    before_projection = np.array([weighted_expected(s, observation) for s in outputs])
    projection_max = 0.
    for j, s in enumerate(outputs):
        residual_j = linear[j] - before_projection[j]
        v = s['values'].copy()
        cells = weighted_cells(s, observation)
        target = chlorophyll(v)
        if len(cells) == 1:
            (i, w), = cells.items()
            target[i] = linear[j] / w
            projection_max = max(projection_max, abs(float(residual_j)))
        else:
            norm = sum(w * w for w in cells.values())
            for i, w in cells.items():
                target[i] += w * residual_j / norm
                projection_max = max(projection_max, abs(float(w * residual_j / norm)))
        v = allocate_total(v, target)
        outputs[j] = dict(s, values=v)
        if not np.isfinite(v).all() or (v[:, 1:] > 10000).any():
            raise AnalysisError('Operator projection produced a nonfinite or excessive state.')
    realized = np.array([weighted_expected(s, observation) for s in outputs])
    residual = float(np.max(np.abs(realized - linear)))
    if residual > config.observation_consistency_tolerance_sigma * sigma:
        raise AnalysisError(f'Native observation operator mismatch: {residual:g} > {config.observation_consistency_tolerance_sigma * sigma:g}')
    report = _report('Chla_weighted', observation, y, sigma, h, realized, linear, before_projection, projection_max, config, states, outputs)
    report.update(w1=w1, w2=w2, weighted_depths_m=weighted)
    return outputs, report
