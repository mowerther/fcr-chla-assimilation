"""Run the 128-member ensemble with sequential updates at the scheduled event hours.

At most two members integrate at a time. All members pause at each event hour, the
filter updates the paused states, and the members continue to the next event. The
forecast and accepted states of every member are stored at every event.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv
import gzip
import json
import time
import numpy as np
from . import native
from .archive import archive, compress
from .constraint import update, Configuration

ROOT = Path(__file__).resolve().parents[2]
TEMPERATURE_SD = .2
CHLA_RELATIVE_SD = .05


def read_csv(path):
    return list(csv.DictReader(Path(path).open()))


def run(configuration, temperature, chla, members_folder, out):
    """configuration: name of the run; temperature, chla: observation rows to assimilate."""
    events = read_csv(ROOT / 'config/events.csv')
    members = [r['member_id'] for r in read_csv(ROOT / 'config/members.csv')]
    hours = [int(e['hour']) for e in events]
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'native').mkdir()
    jobs, used, reports = [], [], []
    started = time.monotonic()
    pool = ThreadPoolExecutor(max_workers=2)
    io = ThreadPoolExecutor(max_workers=8)

    def progress(phase, event=None):
        print(json.dumps(dict(configuration=configuration, phase=phase, event=event, events=len(events),
                              elapsed_seconds=round(time.monotonic() - started))), flush=True)

    try:
        progress('spin-up from 1 December 2016')

        def start(member):
            j = native.launch(out / member, members_folder / member, hours)
            jobs.append(j)
            return j, native.wait(j, 'forecast', 0)
        pairs = list(pool.map(start, members))
        jobs = [j for j, s in pairs]
        states = [s for j, s in pairs]
        for e in events:
            i = int(e['event_id'])
            hour = int(e['hour'])
            forecast = states
            proposed = states
            selected = [r for r in temperature if int(r['event_hour']) == hour] + [r for r in chla if int(r['event_hour']) == hour]
            for r in selected:
                q = dict(r)
                q['value'] = float(r['value'])
                q['sigma'] = TEMPERATURE_SD if r['variable'] == 'T' else CHLA_RELATIVE_SD * q['value']
                proposed, report = update(proposed, q, Configuration(variant='bounded_rtps', chl_depth_m=.1))
                used.append(q)
                reports.append(dict(report, event_id=i))
            for state in proposed:
                assert np.min(state['values'][:, 1:]) >= .03 - 1e-9
            accepted = list(io.map(lambda item: native.write_update(*item), [(j, s, a['values']) for j, s, a in zip(jobs, forecast, proposed)]))
            path = out / 'native' / f'event_{i:03d}.npz'
            archive(path, forecast, accepted)
            compress(path)

            def advance(j):
                native.resume(j, i)
                return native.wait(j, 'forecast', i + 1) if i < len(events) - 1 else native.finish(j)
            result = list(pool.map(advance, jobs))
            states = result if i < len(events) - 1 else []
            for j in jobs:
                for kind in ['forecast', 'accepted', 'update', 'continue']:
                    (j['folder'] / 'exchange' / f'{kind}_{i:03d}.bin').unlink(missing_ok=True)
            if i % 10 == 0 or i == len(events) - 1:
                progress('integrating', i + 1)
        counts = {v: sum(r['variable'] == v for r in used) for v in ['T', 'Chla_weighted']}
        manifest = dict(status='COMPLETED', configuration=configuration, members=len(members), events=len(events),
                        observation_counts=counts, elapsed_seconds=time.monotonic() - started, maximum_integrating_members=2)
        for name, data in [('manifest.json', manifest), ('observation_use.json', used), ('analysis_reports.json', reports)]:
            (out / name).write_text(json.dumps(data, indent=2) + '\n')
        progress('completed', len(events))
        return manifest
    except BaseException as error:
        (out / 'failure.json').write_text(json.dumps(dict(status='FAILED', error=repr(error), elapsed_seconds=time.monotonic() - started), indent=2) + '\n')
        for j in jobs:
            if j['process'].poll() is None:
                native.stop(j)
        raise
    finally:
        pool.shutdown(wait=True)
        io.shutdown(wait=True)
        for j in jobs:
            log = j['folder'] / 'model.log'
            if log.exists():
                with gzip.open(log.with_suffix('.log.gz'), 'wb') as f:
                    f.write(log.read_bytes())
                log.unlink()
