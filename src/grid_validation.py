"""Input-only validation shared by the CLI, launcher and grid editor.

No numerical dependencies, coefficient changes or automatic grid repairs.
The string ``auto`` is accepted only before the model resolves its bounds.
"""
import math
from numbers import Real


def validate_phase_grid(config, label='phase_grid', allow_auto=False):
    """Reject invalid fields together, naming their actual supplied values."""
    if not isinstance(config, dict):
        raise ValueError(f'{label}: нужен JSON-объект с параметрами сетки.')
    allowed = {'energy_min_eV', 'energy_max_eV', 'energy_points', 'angular_points'}
    errors = []
    unknown = set(config) - allowed
    if unknown:
        errors.append('неизвестные поля: ' + ', '.join(sorted(unknown)))
    lo, hi = config.get('energy_min_eV', 1e-4), config.get('energy_max_eV', 50.)
    ne, na = config.get('energy_points', 100), config.get('angular_points', 8)
    def positive(value):
        return isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value) and value > 0
    if not positive(lo):
        errors.append(f'energy_min_eV={lo!r}: требуется конечное число > 0')
    if isinstance(hi, str) and hi == 'auto':
        if not allow_auto:
            errors.append('energy_max_eV="auto" не разрешён в готовую численную границу; создавайте расчёт через MixtureSlab/run.py')
    elif not positive(hi):
        errors.append(f'energy_max_eV={hi!r}: требуется конечное число > 0 или "auto" во входном JSON')
    elif positive(lo) and hi <= lo:
        errors.append(f'energy_max_eV={hi!r} должен быть больше energy_min_eV={lo!r}')
    for key, value, minimum in [('energy_points', ne, 16), ('angular_points', na, 4)]:
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            errors.append(f'{key}={value!r}: требуется целое число >= {minimum}')
    if errors:
        values = f'Emin={lo!r}, Emax={hi!r}, energy_points={ne!r}, angular_points={na!r}'
        raise ValueError(f'{label}: неверная энергетическая сетка ({values}). ' + '; '.join(errors) +
                         '. Исправьте «Параметры → EEDF / ионы → Сетки кинетической энергии» или соответствующий блок JSON. Значения не изменены.')


def validate_configuration_grids(config):
    """Check user grids before data preparation; leave automatic bounds intact."""
    if config.get('schema') not in ('MixtureKineticSlab-1.8', 'MixtureKineticSlab-2.1'):
        return
    problems = []
    for label in ('phase_grid', 'heavy_phase_grid'):
        try:
            validate_phase_grid(config.get(label, {}), label, allow_auto=True)
        except ValueError as exc:
            problems.append(str(exc))
    if problems:
        raise ValueError('\n'.join(problems))
