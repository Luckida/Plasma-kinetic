"""Real He I 2^3P dark test: nonionizing decay versus electron ionization."""
import json
from pathlib import Path
import sys
import numpy as np
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'src'))
from mixture_kinetics import MixtureSlab,LINE_WING_WIDTHS
from project_version import VERSION


def configuration(electrons=0.,metastable_decay=False):
    cfg=json.loads((ROOT/'configs/gas_probe.json').read_text(encoding='utf-8'))
    cfg['atomic_data']['species']=['He']
    cfg['gas']['components']=[dict(species='He_0:L004',density_cm3=1e8)]
    if electrons:cfg['gas']['components'].append(dict(species='He+',density_cm3=electrons))
    cfg['initial_cm3']={row['species']:row['density_cm3'] for row in cfg['gas']['components']}
    cfg['heavy_initial']={'He':{'Ti_K':300.}}
    cfg['diagnostics']['reaction_levels']=['He_0:L004','He_0:L002','He+']
    cfg['electron_initial']['density_cm3']=electrons
    cfg['cells']=1;cfg['level_selection']['mode']='all';cfg['chemistry']['enabled']=False
    cfg['processes']={key:False for key in ['electron_bound_bound','electron_ionization','spontaneous','two_photon','photo_bound_bound','photo_bound_free','radiative_recombination','stimulated_recombination']}
    cfg['processes'].update(spontaneous=True,photo_bound_free=True,electron_ionization=electrons>0)
    if not metastable_decay:cfg['processes']['disabled_ids']=['A:RAD_HeI_2_3S__HeI_1_1S']
    cfg['radiation'].update(enabled=False,scale=0.)
    cfg['spatial_sources']=[];cfg['probe']={'mode':'off'}
    cfg['solver'].update(t_end_s=1e-8,max_step_s=1e-9,samples=3)
    return cfg


def check(electrons=0.,metastable_decay=False):
    cfg=configuration(electrons,metastable_decay);model=MixtureSlab(cfg,ROOT/'configs')
    i=model.index['He_0:L004'];j=model.index['He_1:L001'];gap=model.internal[j]-model.internal[i]
    field=model.field(0.)['mean_4piJ'][0];rates=model.rates(0,field)[0]
    pi=next(k for k,row in enumerate(model.channels) if row['id']=='PI:He_0:L004')
    ion=next((k for k,row in enumerate(model.channels) if row['id']=='ION:He_0:L004'),None)
    initial_ions=float(model.n[0,j]);model.solve()
    result=dict(electrons_initial_cm3=electrons,metastable_decay=metastable_decay,configuration=cfg,ionization_threshold_eV=float(gap),
        external_light_max=float(model.source.max()),initial_photoionization_frequency_s=float(rates[pi]),
        initial_electron_ionization_frequency_s=float(rates[ion]) if ion is not None else 0.,
        electron_ionization_rate_coefficient_cm3_s=float(rates[ion]/electrons) if ion is not None else 0.,
        ions_initial_cm3=initial_ions,ions_final_cm3=float(model.n[0,j]),electrons_final_cm3=float(model.F.sum()),
        ion_events_from_2_3P_cm3={row['id']:float(model.integral[0,k]) for k,row in enumerate(model.channels) if row['i']==i and row['j']==j},
        energy_residual_relative=max(abs(s['energy']) for s in model.snapshots)/model.initial_energy,
        maximum_nuclear_residual=max(float(abs(s['nuclei']).max()) for s in model.snapshots),
        maximum_charge_residual=max(abs(s['charge']) for s in model.snapshots))
    assert model.source.max()==0.,'External source is not off'
    assert rates[pi]==0.,'Nonionizing line created photoionization'
    if electrons:
        assert rates[ion]>0. and model.integral[0,ion]>0.,'Electron impact was accidentally disabled'
    elif metastable_decay:
        assert model.n[0,j]>0.,'Real ionizing M1 photons were accidentally suppressed'
    else:
        assert model.n[0,j]==0. and model.F.sum()==0.,'Dark isolated decay created ions/electrons'
    assert result['energy_residual_relative']<1e-8,'Energy balance failed'
    return result


def main():
    with threadpool_limits(limits=1):cases=[check(),check(1e8),check(metastable_decay=True)]
    report=dict(version=VERSION,scope='Real stored He LS levels, A, CCC electron cross sections and photo tables; implementation regression, not an experimental validation.',
        line_wing_widths=LINE_WING_WIDTHS,cases=cases)
    folder=ROOT/'validation';folder.mkdir(exist_ok=True)
    path=folder/'dark_helium_2_8_3.json';path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    for row in cases:
        print(f"ne={row['electrons_initial_cm3']:.0e} cm^-3, M1={row['metastable_decay']}: initial PI={row['initial_photoionization_frequency_s']:.6g} s^-1; electron ION={row['initial_electron_ionization_frequency_s']:.6g} s^-1; final He+={row['ions_final_cm3']:.12g} cm^-3")
    print('Real helium dark regression passed; report: '+str(path.relative_to(ROOT)))


if __name__=='__main__':main()
