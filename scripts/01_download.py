"""Download the model, runtime libraries and measurements listed in data/sources.json.

Every file is verified against its SHA-256. Archives are unpacked only for the listed
members, which are verified again. Files go to downloads/ and runtime/ (both ignored by git).
"""
from pathlib import Path
import bz2
import hashlib
import io
import json
import sys
import tarfile
import time
import urllib.error
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = json.loads((ROOT / 'data/sources.json').read_text())


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url, digest, target):
    target = Path(target)
    if target.exists() and sha256(target.read_bytes()) == digest:
        return target.read_bytes()
    target.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={'User-Agent': 'fcr-chla-assimilation'})
    for attempt in range(1, 7):
        print('Downloading', url, flush=True)
        try:
            with urllib.request.urlopen(request, timeout=900) as response:
                data = response.read()
            break
        except (urllib.error.URLError, TimeoutError) as error:
            if attempt == 6:
                sys.exit(f'Download failed after {attempt} attempts: {url} ({error})')
            print(f'  attempt {attempt} failed ({error}), retrying in {30 * attempt} s', flush=True)
            time.sleep(30 * attempt)
    if sha256(data) != digest:
        sys.exit(f'SHA-256 mismatch for {url}')
    target.write_bytes(data)
    return data


def members_of(kind, data):
    if kind == 'zip':
        z = zipfile.ZipFile(io.BytesIO(data))
        return lambda name: z.read(name)
    if kind == 'tar.bz2':
        t = tarfile.open(fileobj=io.BytesIO(bz2.decompress(data)))
        return lambda name: t.extractfile(name).read()
    raise ValueError(kind)


def main():
    for source in SOURCES['sources']:
        for item in source['files']:
            data = fetch(item['url'], item['sha256'], ROOT / 'downloads' / item['path'])
            for member in item.get('extract', []):
                content = members_of(item['archive'], data)(member['member'])
                if sha256(content) != member['sha256']:
                    sys.exit(f'SHA-256 mismatch for {member["member"]}')
                target = ROOT / member['target']
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(content)
                if member.get('executable'):
                    target.chmod(0o755)
    print('All files downloaded and verified.')


if __name__ == '__main__':
    main()
