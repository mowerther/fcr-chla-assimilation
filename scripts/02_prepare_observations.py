"""Prepare the assimilated and withheld measurements from the downloaded Carey Lab files.

Writes to output/observations/:
  temperature.csv          surface temperature (mean of the 0.1 m bins over the top 0.5 m), 18 casts listed below
  chla_samples.csv         laboratory Chla at 0.1 and 1.6 m, 7 May to 9 July 2018
  kd.csv                   diffuse attenuation of light (K_d) for the dates with both depths
  chla_weighted.csv        one light-weighted Chla value per date, 7 dates
  withheld_chla.csv        laboratory Chla at 3.8 and 5 m, 28 May to 9 July 2018 (evaluation only)
  withheld_oxygen.csv      CTD oxygen readings 0-9.2 m, 28 May to 9 July 2018 (evaluation only)
Model hour = hours since 2016-12-01 00:00 local standard time (UTC-5), the model start.
"""
from pathlib import Path
import math
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CAREY = ROOT / 'downloads/carey_lab'
OUT = ROOT / 'output/observations'
ORIGIN = pd.Timestamp('2016-12-01T05:00:00Z')
LAT, LON = 37.30325, -79.8373          # Falling Creek Reservoir, site 50
CTD_ORDER = ['Temp_C', 'DO_mgL', 'DOsat_percent', 'Cond_uScm', 'SpCond_uScm', 'Chla_ugL', 'Turbidity_NTU', 'pH', 'ORP_mV',
             'PAR_umolm2s', 'Salinity', 'DescRate_ms']
# CTD casts whose surface temperature is assimilated (7 May to 9 July 2018).
TEMPERATURE_CASTS = ['050718_fcr50.csv', '051118_fcr50.csv', '051418_fcr50.csv', '051718_fcr50.csv', '052818_fcr50.csv',
                     '053018_fcr50.csv', '053118_fcr50.csv', '060418_fcr50.csv', '060718_fcr50.csv', '061118_fcr50.csv',
                     '061418_fcr50.csv', '061818_fcr50.csv', '061818_fcr50_b.csv', '062118_fcr50.csv', '062518_fcr50.csv',
                     '070218_fcr50.csv', '070518_fcr50.csv', '070918_FCR50.csv']


def model_hour(t, ceil=False):
    x = (t - ORIGIN).total_seconds() / 3600
    return int(np.ceil(x) if ceil else np.floor(x))


def read_cast(path):
    d = pd.read_csv(path)
    d.columns = d.columns.str.strip()
    assert d.DateTime.nunique() == 1
    return d, pd.to_datetime(d.DateTime.iloc[0], utc=True)


def flagged(log, local, variable):
    """True if a CTD maintenance log entry with a nonzero flag covers this variable at this local time."""
    for _, q in log.iterrows():
        if pd.notna(q.Reservoir) and str(q.Reservoir).strip() != 'FCR':
            continue
        if pd.notna(q.Site) and float(q.Site) != 50:
            continue
        if pd.notna(q.SN) and float(q.SN) != 7809:
            continue
        start, end = str(q.start_parameter).strip(), str(q.end_parameter).strip()
        if start not in CTD_ORDER or end not in CTD_ORDER:
            continue
        if not CTD_ORDER.index(start) <= CTD_ORDER.index(variable) <= CTD_ORDER.index(end):
            continue
        parse = lambda v: pd.Timestamp(str(v).replace(' EST', '').strip(), tz='America/New_York') if pd.notna(v) else None
        st, en = parse(q.TIMESTAMP_start), parse(q.TIMESTAMP_end)
        if (st is None or local >= st) and (en is None or local <= en) and int(q.flag) != 0:
            return True
    return False


def solar_elevation(t):
    """NOAA approximate solar position, degrees above the horizon, for a UTC timestamp."""
    doy = t.dayofyear
    hour = t.hour + t.minute / 60 + t.second / 3600
    g = 2 * math.pi / 365 * (doy - 1 + (hour - 12) / 24)
    eqtime = 229.18 * (0.000075 + 0.001868 * math.cos(g) - 0.032077 * math.sin(g) - 0.014615 * math.cos(2 * g) - 0.040849 * math.sin(2 * g))
    decl = (0.006918 - 0.399912 * math.cos(g) + 0.070257 * math.sin(g) - 0.006758 * math.cos(2 * g) + 0.000907 * math.sin(2 * g)
            - 0.002697 * math.cos(3 * g) + 0.00148 * math.sin(3 * g))
    tst = hour * 60 + eqtime + 4 * LON
    ha = math.radians(tst / 4 - 180)
    lat = math.radians(LAT)
    cosz = math.sin(lat) * math.sin(decl) + math.cos(lat) * math.cos(decl) * math.cos(ha)
    return 90 - math.degrees(math.acos(max(-1, min(1, cosz))))


def temperature(casts, log, events):
    rows = []
    for f in casts:
        if f.name not in TEMPERATURE_CASTS:
            continue
        d, t = read_cast(f)
        local = t.tz_convert('America/New_York')
        assert not flagged(log, local, 'Temp_C')
        surface = d[d.Depth_m.between(0, .5, inclusive='left') & d.Temp_C.notna()].copy()
        surface['bin'] = np.floor(surface.Depth_m / .1).astype(int)
        means = surface.groupby('bin').Temp_C.mean()
        assert set(means.index) == set(range(5)), f.name
        rows.append(dict(observation_id='T_' + f.stem, variable='T', sample_date=local.strftime('%Y-%m-%d'), measured_utc=t.isoformat(),
                         event_hour=model_hour(t, True), value=float(means.mean()), depth_m=.5, source_file=f.name))
    assert len(rows) == len(TEMPERATURE_CASTS)
    return pd.DataFrame(rows).sort_values('event_hour').reset_index(drop=True)


def laboratory(path):
    lab = pd.read_csv(path)
    lab['date'] = pd.to_datetime(lab.DateTime).dt.strftime('%Y-%m-%d')
    lab = lab[lab.Reservoir.eq('FCR') & lab.Site.eq(50) & lab.Flag_Chla_ugL.eq(0) & lab.Chla_ugL.gt(0)]
    rows = []
    for (date, depth), g in lab[lab.date.between('2018-05-07', '2018-07-09')].groupby(['date', 'Depth_m']):
        assert g.Flag_DateTime.eq(1).all()
        rows.append(dict(sample_date=date, depth_m=float(depth), value=float(g.Chla_ugL.mean()), replicates=len(g)))
    return pd.DataFrame(rows)


def light_attenuation(casts, secchi_path, dates):
    """K_d per date from CTD PAR profiles, with K_d = 1.7 / Secchi depth where no profile qualifies."""
    fits = []
    for f in casts:
        d = pd.read_csv(f, usecols=lambda c: c.strip() in ('DateTime', 'Depth_m', 'PAR', 'Flag'))
        d.columns = [c.strip() for c in d.columns]
        start = pd.Timestamp(d.DateTime.iloc[0])
        down = d.loc[:d.Depth_m.idxmax()]
        valid = down[(down.Depth_m > 0) & (down.PAR > 0) & (down.Flag == 0)]
        row = dict(cast=f.name, cast_date=str(start.tz_convert('America/New_York').date()), start_utc=start.isoformat(),
                   solar_elevation_deg=solar_elevation(start))
        if len(valid):
            top = valid.loc[valid.Depth_m.idxmin()]
            z0, par0 = float(top.Depth_m), float(top.PAR)
            below = valid[(valid.Depth_m >= z0) & (valid.PAR <= .01 * par0)]
            z1 = float(below.Depth_m.min()) if len(below) else np.inf
            pts = valid[(valid.Depth_m >= z0) & (valid.Depth_m <= min(z1, 3.0))]
            x, y = pts.Depth_m.to_numpy(float), np.log(pts.PAR.to_numpy(float))
            if len(pts) >= 2 and np.ptp(x) > 0:
                slope, intercept = np.polyfit(x, y, 1)
                r2 = 1 - ((y - (slope * x + intercept)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
            else:
                slope, r2 = np.nan, np.nan
            row.update(n_points=len(pts), K_d=-slope, R2=r2)
        row['qualifies'] = bool(row['solar_elevation_deg'] > 0 and row.get('n_points', 0) >= 5 and row.get('R2', 0) >= .9 and row.get('K_d', 0) > 0)
        fits.append(row)
    F = pd.DataFrame(fits)
    out = []
    for date in dates:
        c = F[F.qualifies].copy()
        c['offset_days'] = (pd.to_datetime(c.cast_date) - pd.Timestamp(date)).dt.days
        c = c[c.offset_days.abs() <= 3]
        if not len(c):
            sd = pd.read_csv(secchi_path)
            sd = sd[sd.Reservoir.eq('FCR') & sd.Site.eq(50) & sd.Secchi_m.gt(0)].copy()
            sd['offset_days'] = (pd.to_datetime(sd.DateTime.str[:10]) - pd.Timestamp(date)).dt.days
            sd = sd[sd.offset_days.abs() <= 3].assign(abs_offset=lambda x: x.offset_days.abs()).sort_values(['abs_offset', 'offset_days', 'DateTime'])
            r = sd.iloc[0]
            out.append(dict(date=date, K_d=1.7 / r.Secchi_m, source='Secchi depth', profile_or_reading_date=r.DateTime[:10], R2=np.nan))
            continue
        r = c.assign(abs_offset=c.offset_days.abs()).sort_values(['abs_offset', 'offset_days', 'start_utc']).iloc[0]
        out.append(dict(date=date, K_d=r.K_d, source='PAR profile ' + r.cast, profile_or_reading_date=r.cast_date, R2=r.R2))
    out = pd.DataFrame(out)
    # K_d is kept to 10 significant digits, as in the original analysis.
    out['K_d'] = [float(f'{v:.10g}') for v in out.K_d]
    return out


def weighted(samples, kd):
    """y = w1 Chla(0.1 m) + w2 Chla(1.6 m) with weights exp(-2 K_d z) down to min(ln 100 / (2 K_d), 3 m)."""
    lab = samples.set_index(['sample_date', 'depth_m'])
    rows = []
    for r in kd.itertuples():
        z_max = min(math.log(100) / (2 * r.K_d), 3.)
        W1 = 1 - math.exp(-1.7 * r.K_d)
        W2 = 0. if z_max <= .85 else math.exp(-1.7 * r.K_d) - math.exp(-2 * r.K_d * z_max)
        w1, w2 = W1 / (W1 + W2), W2 / (W1 + W2)
        # Laboratory values carry a date only. They become available at the end of the local day.
        available = pd.Timestamp(r.date, tz='America/New_York') + pd.Timedelta(days=1)
        rows.append(dict(observation_id='chla_weighted_' + r.date, variable='Chla_weighted', sample_date=r.date,
                         event_hour=model_hour(available.tz_convert('UTC'), True),
                         value=w1 * lab.loc[(r.date, .1)].value + w2 * lab.loc[(r.date, 1.6)].value,
                         w1=w1, w2=w2, depth_1_m=.1, depth_2_m=1.6, K_d=r.K_d, z_max=z_max))
    return pd.DataFrame(rows)


def withheld_chla(samples, events):
    rows = []
    for r in samples[samples.depth_m.isin([3.8, 5.0]) & samples.sample_date.between('2018-05-28', '2018-07-09')].itertuples():
        hour = model_hour(pd.Timestamp(r.sample_date, tz='America/New_York').tz_convert('UTC'))
        match = events[events.hour.eq(hour)]
        assert len(match) == 1
        rows.append(dict(date=r.sample_date, depth_m=r.depth_m, observed=r.value, event_id=int(match.event_id.iloc[0])))
    return pd.DataFrame(rows).sort_values(['date', 'depth_m']).reset_index(drop=True)


def withheld_oxygen(casts, log, events):
    times = pd.to_datetime(events.time_utc)
    rows = []
    for f in casts:
        d, t = read_cast(f)
        local = t.tz_convert('America/New_York')
        date = local.strftime('%Y-%m-%d')
        if not '2018-05-28' <= date <= '2018-07-09' or flagged(log, local, 'DO_mgL'):
            continue
        i = (times - t).abs().idxmin()
        if abs((times.iloc[i] - t).total_seconds()) / 60 > 60:
            continue
        g = d[d.Depth_m.ge(0) & d.Depth_m.le(9.4) & d.DO_mgL.notna() & d.DO_mgL.ge(0)]
        g = g[g.Depth_m.ge(0) & g.Depth_m.lt(9.2)]
        for r in g.itertuples():
            rows.append(dict(cast=f.name, date=date, sample_utc=t.isoformat(), event_id=int(events.event_id.iloc[i]),
                             depth_m=float(r.Depth_m), bin=int(np.floor((float(r.Depth_m) + 1e-10) / .2)), observed=float(r.DO_mgL)))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    events = pd.read_csv(ROOT / 'config/events.csv')
    log = pd.read_csv(CAREY / 'ctd_maintenancelog_2013_2025.csv', skipinitialspace=True)
    casts = sorted((CAREY / 'ctd').glob('*.csv'))
    T = temperature(casts, log, events)
    # The event calendar: local midnights from 7 May to 10 July 2018 and the hour of each temperature input.
    midnights = range(model_hour(pd.Timestamp('2018-05-07T04:00Z')), model_hour(pd.Timestamp('2018-07-10T04:00Z')) + 1, 24)
    if sorted(set(midnights) | set(T.event_hour)) != events.hour.tolist():
        sys.exit('config/events.csv differs from the event calendar of the measurements.')
    samples = laboratory(CAREY / 'filt-chla_2014_2024.csv')
    both = sorted(set(samples[samples.depth_m.eq(.1)].sample_date) & set(samples[samples.depth_m.eq(1.6)].sample_date))
    kd = light_attenuation(casts, CAREY / 'secchi_2013_2025.csv', both)
    Y = weighted(samples, kd)
    assert set(Y.event_hour) <= set(events.hour)
    for name, table in [('temperature', T), ('chla_samples', samples), ('kd', kd), ('chla_weighted', Y),
                        ('withheld_chla', withheld_chla(samples, events)), ('withheld_oxygen', withheld_oxygen(casts, log, events))]:
        table.to_csv(OUT / f'{name}.csv', index=False, float_format='%.17g')
    print(f'{len(T)} temperature inputs, {len(Y)} light-weighted Chla values on {", ".join(Y.sample_date)}')


if __name__ == '__main__':
    main()
