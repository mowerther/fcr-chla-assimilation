"""Build the Docker image, compile the state exchange library and write the 128 member configurations.

Requires Docker and the files from scripts/01_download.py. Writes runtime/native_exchange.so and
output/members/member_000 ... member_127.
"""
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from fcr_assimilation.members import write_members  # noqa: E402
from fcr_assimilation.native import IMAGE  # noqa: E402


def run(cmd):
    print(' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    run(['docker', 'build', '--platform', 'linux/amd64', '-t', IMAGE, '-f', str(ROOT / 'docker/Dockerfile'), str(ROOT / 'docker')])
    run(['docker', 'run', '--rm', '--platform', 'linux/amd64', '-v', f'{ROOT}:/work', '-w', '/work', IMAGE,
         'sh', 'native/build.sh', '/work/downloads/glm3', '/work/runtime/native_exchange.so'])
    target = ROOT / 'output/members'
    if target.exists():
        shutil.rmtree(target)
    ids = write_members(target)
    print(f'Wrote {len(ids)} member configurations to {target.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
