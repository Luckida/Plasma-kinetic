"""Element-independent names, tabulated cross sections and stoichiometry.

No database queries, helium constants or implicit rate approximations here.
Level energies in imported packs are relative to that ion's ground state;
joining different charge stages requires explicit ionization energies.
"""
from dataclasses import dataclass
import re
try:
    import numpy as np
except ImportError:
    # Species names are also used by the standard-library-only launcher.
    np=None

SYMBOLS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn '
           'Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce '
           'Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn '
           'Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()
ALIASES = {'водород': 'H', 'hydrogen': 'H', 'гелий': 'He', 'helium': 'He', 'литий': 'Li', 'lithium': 'Li',
           'натрий':'Na','sodium':'Na','углерод':'C','carbon':'C','кислород':'O','oxygen':'O','азот':'N','nitrogen':'N','аргон':'Ar','argon':'Ar'}


def roman(number):
    result = ''
    for value, letters in ((100, 'C'), (90, 'XC'), (50, 'L'), (40, 'XL'), (10, 'X'), (9, 'IX'), (5, 'V'), (4, 'IV'), (1, 'I')):
        count, number = divmod(number, value)
        result += letters * count
    return result


@dataclass(frozen=True)
class Species:
    element: str
    charge: int = 0

    def __post_init__(self):
        if self.element not in SYMBOLS or isinstance(self.charge, bool) or not isinstance(self.charge, int) or not 0 <= self.charge <= SYMBOLS.index(self.element) + 1:
            raise ValueError('Неверный атом/заряд')

    @property
    def key(self):
        return f'{self.element}_{self.charge}'

    @property
    def label(self):
        return f'{self.element} {roman(self.charge + 1)}'

    @classmethod
    def parse(cls, value):
        if not isinstance(value, str):
            raise ValueError('Название вида должно быть строкой: H, Li I, Li+, He II')
        text = value.strip()
        text = ALIASES.get(text.casefold(), text)
        # A lone Cyrillic Н is a frequent keyboard-layout substitution for H.
        if text == 'Н':
            text = 'H'
        match = re.fullmatch(r'([A-Z][a-z]?)(?:\s+([IVXC]+)|\s*(\+\d*|\d+\+|\+{2,}))?', text)
        if not match:
            raise ValueError(f'Неизвестное обозначение {value!r}; используйте H, Li I, Li+, Li II')
        element, stage, ion = match.groups()
        if element not in SYMBOLS:
            raise ValueError('Неизвестный элемент: ' + element)
        if stage:
            choices = {roman(q + 1): q for q in range(SYMBOLS.index(element) + 2)}
            if stage not in choices:
                raise ValueError('Неверная спектроскопическая стадия: ' + stage)
            charge = choices[stage]
        elif ion:
            charge = int(ion.replace('+', '')) if ion.replace('+', '') else len(ion)
        else:
            charge = 0
        return cls(element, charge)


def selected_species(config):
    options = config.get('atomic_data', {})
    if not isinstance(options, dict):
        raise ValueError('atomic_data должен быть объектом')
    names = options.get('species', ['He'])
    if not isinstance(names, list) or not names:
        raise ValueError('atomic_data.species: непустой список названий')
    species = [Species.parse(name) for name in names]
    if len(set(species)) != len(species):
        raise ValueError('Повтор атома/иона в atomic_data.species')
    return species


class CrossSection:
    """sigma(E), cm2: strict validity interval, no guessed high-energy tail.

    E is the energy in the coordinate system declared by the data source.
    A missing interval above threshold is an error, not an implicit zero.
    """
    def __init__(self, energy_eV, sigma_cm2, threshold_eV=0., source=''):
        self.energy = np.asarray(energy_eV, dtype=float)
        self.sigma = np.asarray(sigma_cm2, dtype=float)
        self.threshold = float(threshold_eV)
        self.source = source
        if self.energy.ndim != 1 or self.sigma.shape != self.energy.shape or len(self.energy) < 2 or not np.all(np.isfinite(self.energy)) or not np.all(np.isfinite(self.sigma)) or np.any(np.diff(self.energy) <= 0) or np.any(self.energy < 0) or np.any(self.sigma < 0) or not np.isfinite(self.threshold) or self.threshold < 0:
            raise ValueError('Неверная таблица sigma(E)')

    def at(self, energy_eV):
        energy = np.asarray(energy_eV, dtype=float)
        if not np.all(np.isfinite(energy)) or np.any(energy < 0):
            raise ValueError('Неверная энергия')
        active = energy >= self.threshold
        if np.any(active & ((energy < self.energy[0]) | (energy > self.energy[-1]))):
            raise ValueError(f'Энергия вне области sigma(E) {self.source}: [{self.energy[0]}, {self.energy[-1]}] eV')
        return np.where(active, np.interp(energy, self.energy, self.sigma), 0.)

    def rate(self, energy_eV, electron_nodes_cm3, electron_mass_g=9.1093837139e-28):
        """Return volume event rate / target density, s^-1, from node counts."""
        energy, counts = np.asarray(energy_eV), np.asarray(electron_nodes_cm3)
        if not np.isfinite(electron_mass_g) or electron_mass_g <= 0:
            raise ValueError('Нужна положительная масса электрона')
        if counts.shape != energy.shape or np.any(counts < 0) or not np.all(np.isfinite(counts)):
            raise ValueError('Неверные электронные узлы')
        occupied = counts > 0
        values = self.at(energy[occupied])
        return float(np.sum(counts[occupied] * values * np.sqrt(2 * energy[occupied] * 1.602176634e-12 / electron_mass_g)))


class Network:
    """Generic stoichiometric matrix including explicit electrons and photons.

    Each state declares nuclei={element: count}, charge and a readable label.
    Reaction rate laws are provided separately; this object audits conservation
    and converts event rates to population derivatives for any element mixture.
    """
    def __init__(self, states, reactions):
        self.states, self.reactions = states, reactions
        if not states:
            raise ValueError('Нужен непустой список состояний')
        for state in states:
            if isinstance(state['charge'],bool) or not isinstance(state['charge'],int):
                raise ValueError('Заряд состояния должен быть целым числом e')
            for element,count in state['nuclei'].items():
                if element not in SYMBOLS or isinstance(count,bool) or not isinstance(count,int) or count<0:
                    raise ValueError('Число ядер каждого элемента должно быть целым и неотрицательным')
        ids = [state['id'] for state in states]
        if len(set(ids)) != len(ids):
            raise ValueError('Повтор id состояния')
        self.index = {name: k for k, name in enumerate(ids)}
        self.elements = sorted({element for state in states for element in state['nuclei']})
        self.nuclei = np.array([[state['nuclei'].get(element, 0) for state in states] for element in self.elements], dtype=float).reshape(len(self.elements), len(states))
        self.charge = np.array([state['charge'] for state in states], dtype=float)
        self.matrix = np.zeros((len(states), len(reactions)))
        self.labels = []
        if not np.all(np.isfinite(self.charge)) or np.any(self.nuclei < 0) or not np.all(np.isfinite(self.nuclei)):
            raise ValueError('Неверный заряд/состав состояния')
        for k, reaction in enumerate(reactions):
            sides = []
            for side, sign in (('reactants', -1), ('products', 1)):
                labels = []
                if not reaction[side]:
                    raise ValueError('Пустая сторона реакции')
                for name, count in reaction[side].items():
                    if name not in self.index or isinstance(count, bool) or not isinstance(count, int) or count <= 0:
                        raise ValueError('Неизвестный компонент/стехиометрия: ' + name)
                    j = self.index[name]
                    self.matrix[j, k] += sign * count
                    labels.append((str(count) + ' ' if count != 1 else '') + states[j]['label'])
                sides.append(' + '.join(labels))
            self.labels.append(' → '.join(sides))
        if np.any(self.nuclei @ self.matrix) or np.any(self.charge @ self.matrix):
            raise ValueError('Реакция нарушает сохранение ядер или заряда')

    def derivative(self, event_rates):
        rates = np.asarray(event_rates, dtype=float)
        if rates.shape != (len(self.reactions),) or np.any(rates < 0) or not np.all(np.isfinite(rates)):
            raise ValueError('Неверные скорости событий')
        return self.matrix @ rates
