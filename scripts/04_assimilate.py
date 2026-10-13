"""Run the two configurations of the 128-member ensemble from 1 December 2016 to 10 July 2018.

  temperature        surface temperature updates (18 casts, 7 May to 9 July 2018)
  temperature_chla   the same temperature updates and 7 light-weighted Chla values (7 May to 25 June 2018)

Usage: python scripts/04_assimilate.py [temperature|temperature_chla|all]
Each configuration takes about 10 minutes with two members integrating at a time.
Outputs go to output/runs/<configuration>/.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from fcr_assimilation.runner import run, read_csv  # noqa: E402

OBS = ROOT / 'output/observations'


def main():
    choice = sys.argv[1] if len(sys.argv) > 1 else 'all'
    temperature = read_csv(OBS / 'temperature.csv')
    chla = read_csv(OBS / 'chla_weighted.csv')
    configurations = {'temperature': [], 'temperature_chla': chla}
    for name in configurations if choice == 'all' else [choice]:
        out = ROOT / 'output/runs' / name
        if (out / 'manifest.json').exists():
            print(f'{name}: already completed')
            continue
        run(name, temperature, configurations[name], ROOT / 'output/members', out)


if __name__ == '__main__':
    main()
