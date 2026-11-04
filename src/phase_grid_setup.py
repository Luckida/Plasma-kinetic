"""Resolve numerical heavy-particle energy bounds from the supplied inputs.

No populations, temperatures, sources or rate domains are changed. The
thermal coverage criterion is the existing 8 kT check, with a 10% margin.
"""
from copy import deepcopy
from pathlib import Path
import numpy as np
from energy_input import distribution_kind,energy_support
from particle_kinetics import EV,KB,GridOverflow

AMU=1.66053906892e-24


def resolve_heavy_grid(config,data,root):
    requested=deepcopy(config.get('heavy_phase_grid',dict(energy_min_eV=1e-6,energy_max_eV=5.,energy_points=16,angular_points=4)))
    if requested.get('energy_max_eV',50.)!='auto':return requested
    required=0.;file_minima=[];phase_ranges=[];initial=config.get('heavy_initial',{});level_options={}
    for name,options in initial.items():
        if name in data.index or ':' in name or name.startswith('HeI_'):level_options[data.resolve(name)]=options

    def include(options,mass):
        nonlocal required
        kind=distribution_kind(options)
        if kind!='maxwellian':
            low,high=energy_support(options,root)
            if not np.isfinite(low+high) or not 0<=low<high:raise ValueError('Неверный энергетический диапазон распределения тяжёлых частиц.')
            if kind=='phase_table':phase_ranges.append((low,high))
            else:
                required=max(required,high)
                if low>0:file_minima.append(low)
            return
        beams=options.get('beams',[dict(fraction=1.,Ti_K=options.get('Ti_K',300.),uz_cm_s=options.get('uz_cm_s',0.))])
        for beam in beams:
            if float(beam.get('fraction',0.))<=0:continue
            temp=float(beam.get('Ti_K',300.));speed=float(beam.get('uz_cm_s',0.))
            if not np.isfinite(temp+speed) or temp<=0:raise ValueError('Для автоматической сетки нужны конечная скорость и Ti>0 K.')
            thermal=temp*KB/EV;bulk=mass*speed*speed/(2*EV)
            # A Galilean shift of the thermal speed bound, not a collision model.
            required=max(required,(np.sqrt(bulk)+np.sqrt(8*thermal))**2)

    for j,state in enumerate(data.states):
        options=level_options.get(j,initial.get(state['species'],initial.get(state['element'],{'Ti_K':300.})))
        include(options,float(state['mass_u'])*AMU)
    for source in config.get('spatial_sources',[]):
        if source.get('enabled',True) and source.get('mode') in ('ions','particles'):
            state=data.states[data.resolve(source['species'])];include(source.get('distribution',{}),float(state['mass_u'])*AMU)
    boundary=config.get('particle_boundary',{})
    if boundary.get('kind')=='reservoir':
        for side in ('left_heavy','right_heavy'):
            for name,options in boundary.get(side,{}).items():
                state=data.states[data.resolve(name)];include(options,float(state['mass_u'])*AMU)
    low=float(requested.get('energy_min_eV',1e-6))
    if file_minima:low=min(low,min(file_minima))
    # Retain the previous normal 50 eV grid as the lower default, including
    # room for reaction products. Hot inputs can only increase that bound.
    high=max(50.,1.1*required,10*low)
    if phase_ranges:
        # A phase CSV must cover the entire shared grid. Extending it beyond
        # its supplied support would violate the existing input contract.
        floor=max(a for a,b in phase_ranges);cap=min(b for a,b in phase_ranges)
        low=max(low,floor)
        if required>cap or low>=cap or (file_minima and low>min(file_minima)):
            raise GridOverflow('Фазовые файлы не покрывают общую автоматическую сетку тяжёлых частиц. Согласуйте диапазоны файлов с температурами и источниками.')
        high=cap
    requested.update(energy_min_eV=low,energy_max_eV=float(high))
    return requested
