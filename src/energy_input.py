"""Explicit energy-PDF input for either particle mass; no collision fit here.

The common convention is integral P(E)dE=1, with E in eV. An EEPF input
uses integral sqrt(E)f(E)dE=1 and is converted explicitly, never guessed.
"""
from pathlib import Path
import numpy as np
from spectrum_input import _number,_text_rows,_xlsx_rows


def read_energy_table(path,form='probability'):
    path=Path(path)
    if form not in ('probability','eepf'):raise ValueError('Форма распределения: probability или eepf.')
    source=((row,False) for row in _xlsx_rows(path)) if path.suffix.lower()=='.xlsx' else _text_rows(path)
    pairs={}
    for lineno,(row,comment) in enumerate(source,1):
        if comment or len(row)<2:continue
        try:energy=_number(row[0])
        except (ValueError,TypeError):continue
        try:value=_number(row[1])
        except (ValueError,TypeError):raise ValueError(f'Распределение, строка {lineno}: второй столбец должен быть числом.') from None
        if not np.isfinite(energy+value) or energy<0 or value<0:raise ValueError('Энергии и значения распределения должны быть конечными и неотрицательными.')
        if energy in pairs:raise ValueError(f'Повтор энергии {energy:g} эВ; удалите дубликат.')
        pairs[energy]=value
    if len(pairs)<2:raise ValueError('Нужны два столбца: энергия [эВ] и P(E) [эВ⁻¹] или явно выбранная EEPF [эВ⁻³ᐟ²].')
    x=np.array(sorted(pairs));y=np.array([pairs[e] for e in x])
    if form=='eepf':y=y*np.sqrt(x)
    area=float(np.trapezoid(y,x))
    if not np.isfinite(area) or area<=0:raise ValueError('Интеграл распределения должен быть положительным.')
    return x,y/area


def energy_moments(x,y):
    """Integral and first moment of the piecewise-linear input, exactly."""
    dx=np.diff(x);area=float(np.sum(dx*(y[:-1]+y[1:])/2))
    moment=float(np.sum(dx*((2*x[:-1]+x[1:])*y[:-1]+(x[:-1]+2*x[1:])*y[1:])/6))
    return area,moment/area


def energy_phase(grid,path,density,form='probability'):
    """Conservative initial-condition projection, isotropic in the lab frame.

    Two-point Gauss quadrature integrates P(E) and E*P(E) on each linear
    segment exactly. Deposit then preserves their number and energy moments.
    This is an input projection, never a renormalisation of an evolved state.
    """
    from particle_kinetics import GridOverflow
    if not np.isfinite(density) or density<0:raise ValueError('Концентрация >=0.')
    x,y=read_energy_table(path,form)
    if (0<x[0]<grid.e[0]) or x[-1]>grid.e[-1]:raise GridOverflow(f'Сетка частиц должна покрывать весь файл: {x[0]:g}…{x[-1]:g} эВ. Измените Emin/Emax во вкладке распределений.')
    # Split at calculation nodes as well as source knots, so the discrete
    # deposit is also exact for a linear PDF in every energy bin.
    knots=np.unique(np.r_[x,grid.e[(grid.e>x[0])&(grid.e<x[-1])]])
    lo,hi=knots[:-1],knots[1:];mid=(hi+lo)/2;half=(hi-lo)/2
    points=(mid[:,None]+half[:,None]*np.array([-1,1])/np.sqrt(3)).ravel()
    amounts=np.repeat(half,2)*np.interp(points,x,y)*density
    phase,defect=grid.deposit(amounts,points)
    if abs(defect)>1e-10*max(float(amounts@points),1e-300):raise GridOverflow('Сетка обрезает энергию входного распределения: уменьшите Emin или увеличьте Emax.')
    return phase


def distribution_kind(options):
    kind=options.get('kind') or ('energy_table' if options.get('energy_distribution_file') else 'phase_table' if options.get('distribution_file') else 'maxwellian')
    if kind not in ('maxwellian','energy_table','phase_table'):raise ValueError('Вид распределения: maxwellian, energy_table, phase_table.')
    return kind


def energy_support(options,root):
    """Exact file bounds used by automatic electron energy-grid sizing."""
    if distribution_kind(options)=='energy_table':
        x,_=read_energy_table(Path(root)/options['energy_distribution_file'],options.get('energy_form','probability'))
        return float(x[0]),float(x[-1])
    if distribution_kind(options)=='phase_table':
        import csv
        with (Path(root)/options['distribution_file']).open(encoding='utf-8-sig',newline='') as stream:
            values=[float(row['energy_eV']) for row in csv.DictReader(stream,delimiter=';')]
        if not values:raise ValueError('Пустой фазовый CSV.')
        return min(values),max(values)
    return None
