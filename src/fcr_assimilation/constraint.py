"""Common feasible scale on the Chla increments, without clipping any state.

When the full update would take total Chla below its native minimum in any member,
the mean and perturbation corrections of all members are reduced by one common
factor, found by bisection (26 steps).
"""
import numpy as np
from . import analysis

Configuration = analysis.Configuration
_transform = analysis.scalar_transform
_PHYSICAL = ('RTPS input is physically infeasible.', 'Total chlorophyll proposal below the native minimum.')


def update(states, observation, config=Configuration()):
    attempts = []

    def trial(scale):
        def scaled(x, h, y, sigma, rho=1.):
            target = _transform(x, h, y, sigma, rho)
            return target if scale == 1. else np.asarray(x) + scale * (target - np.asarray(x))
        analysis.scalar_transform = scaled
        try:
            outputs, report = analysis.update(states, observation, config)
            minimum = float(min((s['values'][:, 1:] @ analysis.FACTORS).min() for s in outputs))
            phy_min = float(min(s['values'][:, 1:].min() for s in outputs))
            if minimum < analysis.CHL_FLOOR - 1e-10 or phy_min < analysis.PHY_FLOOR - 1e-9:
                raise analysis.AnalysisError('Total chlorophyll proposal below the native minimum.')
            attempts.append(dict(scale=float(scale), feasible=True, minimum_chlorophyll=minimum))
            return outputs, report
        except analysis.AnalysisError as exc:
            if str(exc) not in _PHYSICAL:
                raise
            attempts.append(dict(scale=float(scale), feasible=False, error=str(exc)))
            raise
        finally:
            analysis.scalar_transform = _transform

    if observation['variable'] not in ('Chla', 'Chla_weighted'):
        out, report = analysis.update(states, observation, config)
        report.update(increment_scale=1., positivity_constraint_active=False)
        return out, report
    scale = 1.
    try:
        out, report = trial(scale)
    except analysis.AnalysisError:
        trial(0.)
        low, high = 0., 1.
        for _ in range(26):
            middle = (low + high) / 2
            try:
                trial(middle)
                low = middle
            except analysis.AnalysisError:
                high = middle
        scale = low * (1 - 1e-10)
        out, report = trial(scale)
    report.update(increment_scale=float(scale), positivity_constraint_active=scale < 1., positivity_trials=attempts,
                  minimum_accepted_chlorophyll=float(min((s['values'][:, 1:] @ analysis.FACTORS).min() for s in out)))
    return out, report
