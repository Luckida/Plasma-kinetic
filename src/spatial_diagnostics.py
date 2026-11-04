"""Cell-resolved populations, reaction rates and actual calculated spectra."""
import csv
import numpy as np
from atomic_database import write_csv
from particle_kinetics import EV

HPLANCK=6.62607015e-27
CLIGHT=2.99792458e10


def export_spatial(model,folder):
    x=(np.arange(model.cells)+.5)*model.dz
    profiles=[];kinetics=[];lines=[]
    for snap in model.snapshots:
        time=snap['time_s']
        for z in range(model.cells):
            common=dict(time_s=time,cell=z+1,x_cm=float(x[z]),left_cm=float(z*model.dz),right_cm=float((z+1)*model.dz))
            profiles.append(dict(common,id='e',label='e⁻',element='',charge=-1,molecular=False,
                                 **model.ge.moments(snap['F'][z])))
            for j,state in enumerate(model.states):
                profiles.append(dict(common,id=state['id'],label=state['label'],element=state['element'],charge=state['charge'],
                    molecular=bool(state.get('molecular',False)),**model.gh[state['element']].moments(snap['P'][z,j])))
            for k,channel in enumerate(model.channels):
                value=float(snap['events'][z,k])
                if value==0 and snap['integral'][z,k]==0:continue
                kinetics.append(dict(common,id=channel['id'],kind=channel['kind'],reaction=channel['reaction'],
                    events_cm3_s=value,integral_events_cm3=float(snap['integral'][z,k]),source=channel['source']))
            local=[];n=snap['n'][z]
            for line in model.lines:
                i,j=line['i'],line['j'];gap=model.internal[j]-model.internal[i]
                row=dict(common,id=line['id'],lower=model.states[i]['id'],upper=model.states[j]['id'],
                         transition=model.states[j]['label']+' → '+model.states[i]['label'],wavelength_nm=float(line['wavelength_nm']),
                         wavelength_basis='imported LS energies; not NIST-calibrated',atomic_data_source=line['source'],
                         spectral_kind=line['spectral_kind'],emissivity_erg_cm3_s=float(n[j]*float(line['A_s'])*gap*EV) if line['sp'] is not None else 0.,
                         line_center_net_opacity_cm1=0.,maximum_cell_tau=0.,isolated_cell_equivalent_width_nm=0.)
                if line['optical']:
                    up,down,sigma,ratio=line['optical']
                    opacity=(n[i]*ratio if up is not None else 0.)*sigma-(n[j] if down is not None else 0.)*sigma
                    tau=opacity*model.dz;nu0=gap*EV/HPLANCK;center=int(np.argmin(abs(model.nu-nu0)))
                    row['line_center_net_opacity_cm1']=float(opacity[center]);row['maximum_cell_tau']=float(tau.max())
                    wavelength=CLIGHT/model.nu*1e7
                    row['isolated_cell_equivalent_width_nm']=float(-np.trapezoid(-np.expm1(-np.clip(tau,-50,700)),wavelength))
                local.append(row)
            for key,rank in [('emissivity_erg_cm3_s','emission_rank'),('isolated_cell_equivalent_width_nm','absorption_rank')]:
                for number,index in enumerate(sorted(range(len(local)),key=lambda k:-local[k][key]),1):local[index][rank]=number
            lines.extend(local)
    write_csv(folder/'spatial_profiles.csv',profiles,list(profiles[0]))
    write_csv(folder/'spatial_kinetics.csv',kinetics,list(kinetics[0]) if kinetics else ['time_s','cell','x_cm','left_cm','right_cm','id','kind','reaction','events_cm3_s','integral_events_cm3','source'])
    line_fields=list(lines[0]) if lines else ['time_s','cell','x_cm','id','wavelength_nm','emissivity_erg_cm3_s','isolated_cell_equivalent_width_nm','emission_rank','absorption_rank']
    write_csv(folder/'lines_by_position.csv',lines,line_fields)
    top_n=int(model.cfg.get('diagnostics',{}).get('spatial_top_lines',10))
    if top_n<1:raise ValueError('diagnostics.spatial_top_lines>=1')
    write_csv(folder/'strongest_lines_by_position.csv',[r for r in lines if (r['emission_rank']<=top_n and r['emissivity_erg_cm3_s']>0) or (r['absorption_rank']<=top_n and r['isolated_cell_equivalent_width_nm']>0)],line_fields)
    field=model.snapshots[-1]['field'];opacity=field['opacity'];emission=field['emission']
    tau=opacity*model.dz;left=np.cumsum(tau,axis=0)-tau/2;right=np.cumsum(tau[::-1],axis=0)[::-1]-tau/2
    fields=['cell','x_cm','wavelength_A','net_opacity_cm1','atomic_emission_erg_cm3_s_Hz','mean_4piJ_erg_cm2_s_Hz','cell_tau','tau_to_left','tau_to_right']
    with (folder/'spectra_by_position_final.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fields,delimiter=';');writer.writeheader()
        for z in range(model.cells):
            for k,nu in enumerate(model.nu):
                writer.writerow(dict(cell=z+1,x_cm=float(x[z]),wavelength_A=float(CLIGHT/nu*1e8),net_opacity_cm1=float(opacity[z,k]),
                    atomic_emission_erg_cm3_s_Hz=float(emission[z,k]),mean_4piJ_erg_cm2_s_Hz=float(field['mean_4piJ'][z,k]),
                    cell_tau=float(tau[z,k]),tau_to_left=float(left[z,k]),tau_to_right=float(right[z,k])))
    np.savez_compressed(folder/'spatial_radiation_history.npz',time_s=[s['time_s'] for s in model.snapshots],x_cm=x,
                        wavelength_A=CLIGHT/model.nu*1e8,opacity_cm1=[s['field']['opacity'] for s in model.snapshots],
                        atomic_emission_erg_cm3_s_Hz=[s['field']['emission'] for s in model.snapshots],mean_4piJ_erg_cm2_s_Hz=[s['field']['mean_4piJ'] for s in model.snapshots])
