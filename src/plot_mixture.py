"""Three compact plots from computed mixture outputs; no resimulation/fits."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def spatial_plots(folder,states,n,fe,note):
    with np.load(folder/'mixture_history.npz',allow_pickle=False) as data:
        if 'position_cm' not in data:return
        x=data['position_cm'];times=data['time_s']
    if len(x)<2:return
    fig,axes=plt.subplots(2,1,figsize=(10,8),sharex=True)
    groups=sorted({s['species'] for s in states})
    for name in groups:
        ids=[j for j,s in enumerate(states) if s['species']==name];values=n[-1,:,ids].sum(axis=0)
        if values.max()>0:axes[0].plot(x,values,marker='.',label=name)
    axes[0].plot(x,fe[-1].sum(axis=(1,2)),color='black',ls='--',label='e⁻')
    axes[0].set(ylabel='n, см⁻³',yscale='log',title='Концентрации компонентов в конце расчёта');axes[0].legend(fontsize=8,ncol=4);axes[0].grid(alpha=.2)
    metastables=[j for j,s in enumerate(states) if s['element']=='He' and s['charge']=='0' and s['n']=='2' and s['L']=='0' and s['multiplicity']=='3']
    chosen=metastables[0] if metastables else int(np.argmax(n[-1].max(axis=0)))
    for index in sorted({0,len(times)//2,len(times)-1}):axes[1].plot(x,n[index,:,chosen],marker='.',label=f't={times[index]:.3g} с')
    axes[1].set(xlabel='x, см',ylabel='n, см⁻³',yscale='log',title=states[chosen]['label']);axes[1].legend();axes[1].grid(alpha=.2)
    fig.suptitle(note);fig.tight_layout();fig.savefig(folder/'spatial_profiles.png',dpi=160);plt.close(fig)
    with (folder/'lines_by_position.csv').open(encoding='utf-8-sig',newline='') as f:lines=list(csv.DictReader(f,delimiter=';'))
    final=[row for row in lines if float(row['time_s'])==times[-1]]
    if final:
        fig,axes=plt.subplots(2,1,figsize=(10,8),sharex=True)
        for ax,key,label in [(axes[0],'emissivity_erg_cm3_s','ε, эрг см⁻³ с⁻¹'),(axes[1],'isolated_cell_equivalent_width_nm','W ячейки, нм')]:
            scores={ident:max(float(r[key]) for r in final if r['id']==ident) for ident in {r['id'] for r in final}}
            for ident in sorted(scores,key=lambda k:-scores[k])[:5]:
                if scores[ident]<=0:continue
                rows=sorted([r for r in final if r['id']==ident],key=lambda r:float(r['x_cm']))
                ax.plot([float(r['x_cm']) for r in rows],[float(r[key]) for r in rows],marker='.',label=f"{float(rows[0]['wavelength_nm']):.3f} нм · {ident}")
            if ax.lines:ax.set_yscale('log');ax.legend(fontsize=8)
            ax.set_ylabel(label);ax.grid(alpha=.2)
        axes[0].set_title('Сильнейшие линии модели LS · λ из импортированных энергий ADAS')
        axes[1].set_xlabel('x, см');fig.suptitle(note);fig.tight_layout();fig.savefig(folder/'spatial_lines.png',dpi=160);plt.close(fig)
    with (folder/'chemical_event_history.csv').open(encoding='utf-8-sig',newline='') as f:history=list(csv.DictReader(f,delimiter=';'))
    final=[row for row in history if float(row['time_s'])==times[-1]]
    if final:
        fig,ax=plt.subplots(figsize=(10,4.5));scores={ident:max(float(r['events_cm3_s']) for r in final if r['id']==ident) for ident in {r['id'] for r in final}}
        for ident in sorted(scores,key=lambda k:-scores[k])[:6]:
            if scores[ident]<=0:continue
            rows=sorted([r for r in final if r['id']==ident],key=lambda r:int(r['cell']))
            ax.plot([x[int(r['cell'])-1] for r in rows],[float(r['events_cm3_s']) for r in rows],marker='.',label=ident)
        if ax.lines:ax.set_yscale('log');ax.legend(fontsize=8)
        ax.set(xlabel='x, см',ylabel='События, см⁻³ с⁻¹',title='Ведущие межкомпонентные реакции');ax.grid(alpha=.2)
        fig.suptitle(note);fig.tight_layout();fig.savefig(folder/'spatial_chemistry.png',dpi=160);plt.close(fig)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path);args=parser.parse_args(argv)
    folder=args.output
    with (folder/'states.csv').open(encoding='utf-8-sig',newline='') as f:states=list(csv.DictReader(f,delimiter=';'))
    quality=json.loads((folder/'quality.json').read_text(encoding='utf-8-sig'));note='Неполные атомные данные — см. data_gaps.csv' if quality['atomic_data_status']=='incomplete' else 'Скалярная кинетическая модель'
    with np.load(folder/'mixture_history.npz',allow_pickle=False) as z:
        t=z['time_s'];n=z['populations_cm3'];fe=z['electrons_cm3'];p=z['heavy_phase_cm3'];e=z['energy_eV'];we=z['energy_weights_eV'];he=z['heavy_energy_eV'];mu=z['heavy_mu']
    cells=n.shape[1];fig,axes=plt.subplots(cells,1,figsize=(10,3.8*cells),squeeze=False)
    candidates=[j for j,s in enumerate(states) if s.get('active')=='True']
    scores=n.max(axis=(0,1));active=sorted(candidates,key=lambda j:-scores[j])[:12]
    for cell,ax in enumerate(axes[:,0]):
        for j in active:
            values=n[:,cell,j];good=values>0
            if good.any():ax.plot(t[good],values[good],label=states[j]['label'])
        ne=fe[:,cell].sum(axis=(1,2));good=ne>0
        if good.any():ax.plot(t[good],ne[good],color='black',ls='--',lw=2,label='e⁻')
        ax.set(xlabel='t, с',ylabel='n, см⁻³',yscale='log',title=f'Ячейка {cell+1}');ax.grid(alpha=.2)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.suptitle(note+' · основные 12 уровней; все уровни в CSV/NPZ',fontsize=10,y=.995)
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.977),fontsize=8,ncol=3)
    fig.tight_layout(rect=(0,0,1,.88));fig.savefig(folder/'populations.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,cells,figsize=(6*cells,4.2),squeeze=False)
    for cell,ax in enumerate(axes[0]):
        distributions=[]
        for index in [0,-1]:
            counts=fe[index,cell].sum(axis=1);distributions.append(counts/max(counts.sum(),1e-300)/we)
        peak=max(max(float(pdf.max()) for pdf in distributions),1e-300);floor=peak*1e-12
        for index,style,label in [(0,'--','начало'),(-1,'-','конец')]:
            counts=fe[index,cell].sum(axis=1);pdf=counts/max(counts.sum(),1e-300)/we;good=pdf>=floor
            ax.plot(e[good],pdf[good],style,label=label)
        # Display floor only; the exported distribution is never truncated.
        ax.set_ylim(floor,peak*2)
        ax.set(xlabel='E, эВ',ylabel='EEDF, эВ⁻¹',yscale='log',xscale='log',title=f'Ячейка {cell+1}');ax.grid(alpha=.2);ax.legend()
    fig.suptitle(note);fig.tight_layout();fig.savefig(folder/'electron_eedf.png',dpi=160);plt.close(fig)
    groups=sorted({s['species'] for s in states if int(s['charge'])>0 and any(n[:,:,j].max()>0 for j,r in enumerate(states) if r['species']==s['species'])})
    if groups:
        fig,axes=plt.subplots(len(groups),1,figsize=(9,3.3*len(groups)),squeeze=False)
        for name,ax in zip(groups,axes[:,0]):
            indices=[j for j,s in enumerate(states) if s['species']==name];mass=float(states[indices[0]]['mass_u'])*1.66053906892e-24
            velocity=(np.sqrt(2*he[:,None]*1.602176634e-12/mass)*mu[None,:]/1e5).ravel()
            combined=sum(p[index,0,indices].sum(axis=0) for index in [0,-1]).ravel();order=np.argsort(abs(velocity));cdf=np.cumsum(combined[order])
            cut=min(np.searchsorted(cdf,cdf[-1]*(1-1e-5)),len(order)-1) if cdf[-1]>0 else 0
            vmax=max(float(abs(velocity[order[cut]])),1e-6);bins=np.linspace(-vmax,vmax,61)
            for index,style,label in [(0,'--','начало'),(-1,'-','конец')]:
                counts=p[index,0,indices].sum(axis=0).ravel();hist,edges=np.histogram(velocity,bins,weights=counts)
                ax.stairs(hist,edges,linestyle=style,label=label)
            ax.set(xlabel='v_z, км/с',ylabel='n в интервале, см⁻³',title=name+' · ячейка 1');ax.grid(alpha=.2);ax.legend()
        fig.suptitle(note);fig.tight_layout();fig.savefig(folder/'ion_velocity.png',dpi=160);plt.close(fig)
    if (folder/'spatial_profiles.csv').exists():spatial_plots(folder,states,n,fe,note)
    from probe_spectrum import plot_probe
    plot_probe(folder)
    print('Графики сохранены:',folder,flush=True)


if __name__=='__main__':main()
