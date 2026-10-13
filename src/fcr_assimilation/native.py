"""Run the published GLM-AED executable in Docker and exchange its state at event hours.

Each member runs in its own container. The library native/native_exchange.c (loaded with
LD_PRELOAD) pauses the model at each event hour, writes the forecast state, waits for the
updated temperature and phytoplankton, writes the accepted state and waits for a continue
signal. All other model memory stays in the paused process.
"""
from pathlib import Path
import json
import shutil
import subprocess
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUNTIME = ROOT / 'runtime'
IMAGE = 'fcr-chla-assimilation:glm'
MAGIC = 0x46435258


def launch(folder, member_source, event_hours):
    """Start one member in folder, from the configuration files in member_source."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    (folder / 'exchange').mkdir()
    (folder / 'inputs').mkdir()
    shutil.copytree(Path(member_source) / 'aed', folder / 'aed')
    shutil.copy2(Path(member_source) / 'glm3.nml', folder / 'glm3.nml')
    (folder / 'event_hours.txt').write_text('\n'.join(str(h) for h in event_hours) + '\n')
    rel = folder.relative_to(ROOT)
    name = 'fcr-' + str(rel).replace('/', '-').replace('_', '-')
    cmd = ['docker', 'run', '--rm', '--network', 'none', '--platform', 'linux/amd64', '--cpus', '1', '--memory', '192m',
           '--name', name, '-e', 'OMP_NUM_THREADS=1', '-e', 'LD_LIBRARY_PATH=/work/runtime/lib',
           '-v', f'{ROOT}:/work',
           '--mount', f'type=bind,source={RUNTIME}/forcing,target=/work/{rel}/inputs,readonly',
           '-w', f'/work/{rel}',
           '-e', 'LD_PRELOAD=/work/runtime/native_exchange.so',
           '-e', f'GLM_EXCHANGE_DIR=/work/{rel}/exchange',
           '-e', f'GLM_EVENT_HOURS=/work/{rel}/event_hours.txt',
           '-e', 'GLM_SKIP_OUTPUT=1',
           IMAGE, '/work/runtime/glm+', '--nml', 'glm3.nml', '--quiet', '2']
    (folder / 'command.json').write_text(json.dumps(cmd, indent=2) + '\n')
    with (folder / 'model.log').open('w') as log:
        process = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
    return dict(process=process, folder=folder, start=time.monotonic(), name=name)


def read(path):
    with Path(path).open('rb') as f:
        h = np.fromfile(f, dtype='<i4', count=8)
        assert h[0] == MAGIC
        n, nv, ml = map(int, h[2:5])
        g = np.fromfile(f, dtype='<f8', count=n * 7).reshape(n, 7)
        wq = np.fromfile(f, dtype='<f8', count=n * nv).reshape(nv, n)
    ids = h[5:8].astype(int)
    return dict(record=int(h[1]), n=n, ids=ids, geometry=g, wq=wq, values=np.column_stack([g[:, 2], wq[ids].T]))


def wait(job, kind, rec, timeout=500):
    p = job['folder'] / 'exchange' / f'{kind}_{rec:03d}.bin'
    start = time.monotonic()
    while not p.exists():
        if job['process'].poll() is not None:
            raise RuntimeError(job['name'] + ' exited: ' + (job['folder'] / 'model.log').read_text()[-2500:])
        if time.monotonic() - start > timeout:
            raise TimeoutError(str(p))
        time.sleep(.001)
    return read(p)


def write_update(job, state, values):
    p = job['folder'] / 'exchange' / f'update_{state["record"]:03d}.bin'
    tmp = p.with_suffix('.partial')
    with tmp.open('wb') as f:
        np.array([state['record'], state['n']], dtype='<i4').tofile(f)
        np.asarray(values, dtype='<f8').tofile(f)
    tmp.replace(p)
    accepted = wait(job, 'accepted', state['record'])
    # Only temperature (and its density) and the three phytoplankton pools may change.
    other = np.delete(np.arange(state['wq'].shape[0]), state['ids'])
    np.testing.assert_array_equal(state['wq'][other], accepted['wq'][other])
    np.testing.assert_array_equal(state['geometry'][:, [0, 1, 3, 5, 6]], accepted['geometry'][:, [0, 1, 3, 5, 6]])
    np.testing.assert_array_equal(values, accepted['values'])
    return accepted


def resume(job, rec):
    (job['folder'] / 'exchange' / f'continue_{rec:03d}.bin').touch()


def finish(job):
    rc = job['process'].wait(timeout=500)
    if rc:
        raise RuntimeError(f'{job["name"]} exited {rc}: ' + (job['folder'] / 'model.log').read_text()[-1200:])
    return dict(name=job['name'], wall_seconds=time.monotonic() - job['start'], returncode=rc)


def stop(job):
    subprocess.run(['docker', 'stop', '-t', '1', job['name']], capture_output=True)
