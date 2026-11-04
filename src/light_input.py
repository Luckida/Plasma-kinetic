"""A spectral shape plus an explicit normal power; no new kinetic operator."""
from math import isfinite
from pathlib import Path
from spectrum_input import Spectrum,load_spectrum,HC_EV_A,EV_ERG


def integral(x,y):return sum((x2-x1)*(y1+y2)/2 for x1,x2,y1,y2 in zip(x,x[1:],y,y[1:]))


def shape_spectrum(cfg,root,probe=None):
    mode=cfg.get('spectrum_mode','file')
    if mode=='probe':
        p=probe or {};mode=p.get('mode','off')
        if mode=='off':raise ValueError('Для световой накачки выберите пробную форму flat/file или отдельный файл.')
        cfg=p
    if mode=='flat':
        low=float(cfg.get('minimum_nm',1082.));high=float(cfg.get('maximum_nm',1084.));count=int(cfg.get('points',501))
        if not 0<low<high or not 2<=count<=200000:raise ValueError('Форма света: 0<λ_min<λ_max, от 2 до 200000 точек.')
        x=tuple(10*(low+(high-low)*i/(count-1)) for i in range(count));return Spectrum(x,(1.,)*count,dict(input_x_unit='nm',input_y_kind='relative',warnings=[]))
    if mode!='file':raise ValueError('Форма света: file, flat или probe.')
    options=dict(cfg.get('spectrum_options',{}))
    # Relative spectral densities are per input coordinate. Apply its Jacobian
    # before normalising, including when X is energy or frequency.
    if options.get('y_kind','relative') in ('auto','relative'):options['y_kind']='energy'
    return load_spectrum(Path(root)/cfg.get('spectrum_file',''),options,purpose='probe')


def light_spectrum(cfg,root,probe=None):
    """Canonical 4pi J_lambda; multiplier_W_cm2 is normal beam power."""
    normalised=cfg.get('normalization')=='integral'
    if not normalised:return load_spectrum(Path(root)/cfg['spectrum_file'],cfg.get('spectrum_options'),purpose='driving')
    power=float(cfg.get('multiplier_W_cm2',0.));mu=float(cfg.get('mu',cfg.get('beam_mu',1.)))
    if not isfinite(power+mu) or power<0 or not 0<mu<=1:raise ValueError('Множитель света — мощность >=0 Вт/см²; 0<μ<=1.')
    data=shape_spectrum(cfg,root,probe);x=data.wavelength_A;y=data.values
    if any(v<0 or not isfinite(v) for v in y):raise ValueError('Форма световой накачки должна быть неотрицательной; исправьте фон спектра.')
    area=integral(x,y)
    if not isfinite(area) or area<=0:raise ValueError('Интеграл формы спектра должен быть положительным.')
    values=tuple(v/area*power*1e7/mu for v in y)
    photons=integral(x,tuple(v/(HC_EV_A/lam*EV_ERG) for lam,v in zip(x,values)))*mu
    info=dict(data.metadata,normalization='integral',multiplier_W_cm2=power,normal_power_W_cm2=power,photon_flux_cm2_s=photons,
        minimum_photon_eV=HC_EV_A/max(x),maximum_photon_eV=HC_EV_A/min(x),mean_photon_eV=power*1e7/photons/EV_ERG if photons else 0.,
        spectrum_convention='normal energy flux; integral shape = 1, before multiplication',mu=mu)
    return Spectrum(x,values,info)


def light_rows(cfg,root,probe=None):return light_spectrum(cfg,root,probe).canonical_rows()


def light_enabled(cfg):
    return bool(cfg.get('enabled',True)) and (float(cfg.get('multiplier_W_cm2',0.))>0 if cfg.get('normalization')=='integral' else float(cfg.get('scale',1.))>0 and bool(cfg.get('spectrum_file')))


def particle_inputs(item):
    """Beam concentration n and signed normal speed u define flux n*abs(u)."""
    from copy import deepcopy
    row=deepcopy(item)
    if row.get('mode')=='radiation' or 'density_cm3' not in row:return row
    if 'velocity_cm_s' not in row:raise ValueError('Для пучка с концентрацией задайте velocity_cm_s, см/с.')
    n=float(row['density_cm3']);u=float(row['velocity_cm_s'])
    if not isfinite(n+u) or n<0:raise ValueError('Концентрация пучка >=0 см⁻³, скорость — конечное число в см/с.')
    strength=n*abs(u)
    if 'rate_cm2_s' in row and abs(float(row['rate_cm2_s'])-strength)>1e-10*max(strength,1.):raise ValueError('Не согласованы концентрация×скорость и старое поле потока; удалите rate_cm2_s.')
    row['rate_cm2_s']=strength;row.setdefault('distribution',{})['uz_cm_s']=u
    return row
