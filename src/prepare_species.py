"""Prepare real available atomic data by species name; report all data gaps.

python prepare_species.py --config elements_H_Li.json
python prepare_species.py --species H Li "Li II" --out atomic_data_trial
This command imports data; it does not run a multi-element plasma simulation.
"""
import argparse
import json
from pathlib import Path
import sys
from atomic_database import prepare_species, ROOT
from progress import progress_from_config


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path)
    parser.add_argument('--species', nargs='+')
    parser.add_argument('--out', '--output', dest='output', type=Path)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--progress', choices=('text', 'json', 'off'))
    args = parser.parse_args(argv)
    root = args.config.resolve().parent if args.config else ROOT
    config = json.loads(args.config.read_text(encoding='utf-8-sig')) if args.config else {}
    if config.get('schema', 'AtomicData-1.0') == 'AtomicData-1.0' and set(config) - {'schema', 'atomic_data', 'progress', 'output_directory'}:
        raise ValueError('Неизвестные поля AtomicData-1.0')
    options = config.setdefault('atomic_data', {})
    if args.species:
        options['species'] = args.species
    if args.output:
        options['prepared_directory'] = str(args.output.resolve())
    if args.offline:
        options['offline'] = True
    progress = progress_from_config(config, args.progress)
    try:
        rows = prepare_species(config, root, progress)
        for row in rows:
            print(f"{row['species']}: уровней {row['levels']}, линий {row['lines']}, σ(E)-каналов {row['cross_section_channels']}; {row['status']}", flush=True)
        print('Это подготовка данных. Полнота и пробелы: process_coverage.csv. Расчёт смеси: python run.py --config configs/mixture_H_He.json.', flush=True)
        return 0 if all(row['levels'] and not row['download_errors'] for row in rows) else 2
    except Exception as exc:
        progress.fail(str(exc))
        raise


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        sys.exit('Ошибка: ' + str(exc))
