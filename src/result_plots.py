"""Plot only requested quantities. Stable reaction IDs, explicit top-five gaps."""
from dataclasses import dataclass,field,replace
import textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from result_store import top_indices

LABELS={'components':'Компоненты и реагенты','levels':'Населённости уровней','electrons':'Электроны','rates':'Скорости реакций','integrals':'Интегралы реакций','top5':'Top-5 реакций','mass':'Сохранение массы','charge':'Сохранение заряда','kinetic':'Проверка кинетики уровней','energy':'Сохранение энергии','eedf':'Энергии электронов','heavy_energy':'Средние энергии частиц','velocity':'Средние скорости частиц','absorption':'Поглощение пробного света','transmission':'Пропускание пробного света','emission':'Спектр излучения','probe_distance':'Пропускание по расстоянию','heavy_distribution':'Распределения скоростей частиц'}
LABELS['components']='Концентрации компонентов и реагентов'
LABELS['kinetic_e']='Проверка кинетики электронов'
LABELS['iedf']='Энергетические распределения ионов / атомов'


@dataclass
class PlotOptions:
    metrics:list=field(default_factory=lambda:['components'])
    levels:list|None=None
    components:list|None=None
    reactions:list|None=None
    reaction_filter:list|None=None
    axis:str='time'
    cell:int=-1
    time_index:int=-1
    target:str='all'
    direction:str='source'
    x_scale:str='linear'
    y_scale:str='linear'
    fontsize:int=16
    curves_per_figure:int=5
    mark_starts:bool=True
    fractions:bool=False
    time_min:float|None=None
    time_max:float|None=None
    wavelength_min:float|None=None
    wavelength_max:float|None=None
    probe_wavelength_nm:float=1083.3


@dataclass
class PlotPage:
    key:str
    title:str
    figure:object
    x:object
    curves:list
    size_inches:tuple=()


def reduce_cell(values,cell):return values.mean(axis=1) if cell<0 else values[:,cell]


def target_weights(data,target,direction):
    if target=='all':return np.ones(len(data.reactions))
    if direction=='participation':
        if target=='electron':members={'e'}
        elif target.startswith('level:'):members={target[6:]}
        else:members={data.states[j]['id'] for j in data.meta['components'][target[10:]]}
        return np.array([bool((set(r.get('reactants',{}))|set(r.get('products',{})))&members) for r in data.reactions],dtype=float)
    if target=='electron':stoich=data['electron_stoichiometry']
    elif target.startswith('level:'):
        ident=target[6:];j=next(j for j,s in enumerate(data.states) if s['id']==ident);stoich=data['stoichiometry'][j]
    elif target.startswith('component:'):
        ids=data.meta['components'][target[10:]];stoich=data['stoichiometry'][ids].sum(axis=0)
    else:raise ValueError('Неизвестная цель реакций: '+target)
    return abs(stoich) if direction=='both' else np.maximum(stoich if direction=='source' else -stoich,0)


def top_history(data,options):
    rates=reduce_cell(data['reaction_rates'],options.cell)*target_weights(data,options.target,options.direction)
    if options.reaction_filter is not None:
        allowed=set(options.reaction_filter);rates=rates*np.array([r['id'] in allowed for r in data.reactions])
    masks=np.zeros_like(rates,dtype=bool)
    for t,row in enumerate(rates):masks[t,top_indices(row)]=True
    selected=np.flatnonzero(masks.any(axis=0));selected=sorted(selected,key=lambda k:-float(rates[:,k].max()))
    return [(int(k),np.where(masks[:,k],rates[:,k],np.nan)) for k in selected]


def start_markers(ax,sources,options):
    if not options.mark_starts:return
    grouped={}
    for row in sources:
        strength=row.get('normal_power_W_cm2',row.get('rate_cm2_s',0.))
        if strength<=0:continue
        grouped.setdefault(float(row['start_s']),[]).append(row['id'])
    left,right=ax.get_xlim()
    for level,(time,names) in enumerate(sorted(grouped.items())):
        label='Старт: '+', '.join(names)+f' · {time:.4g} с'
        if options.x_scale=='log' and time<=0:
            ax.text(.015,.97-.08*level,label,transform=ax.transAxes,ha='left',va='top',fontsize=max(options.fontsize-3,11),bbox=dict(facecolor='white',alpha=.9,edgecolor='none'));continue
        if not left<=time<=right:continue
        line=ax.axvline(time,color='#666666',linestyle='--',linewidth=1.25,alpha=.8);line.set_gid('pump_start')
        ax.annotate(label,xy=(time,.98-.08*level),xycoords=('data','axes fraction'),xytext=(5,0),textcoords='offset points',ha='left',va='top',fontsize=max(options.fontsize-3,11),bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))


def build_plots(data,options=None):
    o=options or PlotOptions();nt=len(data['time_s']);ti=o.time_index%nt
    if o.cell>=len(data['position_cm']):raise ValueError('Выбранной ячейки нет в результате.')
    if not 12<=int(o.fontsize)<=28:raise ValueError('Размер текста должен быть от 12 до 28 pt.')
    if not 1<=o.curves_per_figure<=8:raise ValueError('От 1 до 8 кривых на изображении.')
    if o.time_min is not None and o.time_max is not None and o.time_min>=o.time_max:raise ValueError('Начало диапазона времени должно быть меньше конца.')
    if o.x_scale=='log' and o.time_min is not None and o.time_min<=0:raise ValueError('Для log X начало времени должно быть >0; оставьте поле пустым для автоматического диапазона.')
    if o.wavelength_min is not None and o.wavelength_max is not None and o.wavelength_min>=o.wavelength_max:raise ValueError('Начало спектрального диапазона должно быть меньше конца.')
    n=data['n_cm3'];x=data['time_s'] if o.axis=='time' else data['position_cm'];xlabel='Время, с' if o.axis=='time' else 'Расстояние, см'
    def context(values):return reduce_cell(values,o.cell) if o.axis=='time' else values[ti]
    def curve_label(s):return s.get('display_label',s['label'])
    states={s['id']:j for j,s in enumerate(data.states)}
    level_ids=[states[ident] for ident in o.levels if ident in states] if o.levels is not None else sorted(range(len(data.states)),key=lambda j:-float(n[:,:,j].max()))[:5]
    if not level_ids and any(key in o.metrics for key in ('levels','kinetic','heavy_energy','velocity','heavy_distribution','iedf')):raise ValueError('Выберите хотя бы один уровень во вкладке «Компоненты».')
    components=o.components if o.components is not None else [name for name,ids in data.meta['components'].items() if n[:,:,ids].max()>0]
    if not components and 'components' in o.metrics:raise ValueError('Выберите хотя бы один компонент.')
    reactions={r['id']:k for k,r in enumerate(data.reactions)};selected=[reactions[ident] for ident in o.reactions if ident in reactions] if o.reactions is not None else top_indices(data['reaction_rates'].max(axis=(0,1)))
    rates=data['reaction_rates']
    if not selected and o.reactions is not None and any(key in o.metrics for key in ('rates','integrals')):raise ValueError('Выберите реакции во вкладке «Реакции».')
    pages=[]
    def add(key,title,ylabel,curves,xx=None,xlab=None,time_axis=None,empty_message='Нет значений для выбранного масштаба'):
        xx=np.asarray(x if xx is None else xx);xlab=xlabel if xlab is None else xlab;time_axis=(o.axis=='time') if time_axis is None else time_axis
        cell_label='среднее по слою' if o.cell<0 else f'x={data["position_cm"][o.cell]:.6g} см'
        if time_axis:title+=' · '+cell_label
        elif xlab=='Расстояние, см' and 't=' not in title:title+=f' · t={data["time_s"][ti]:.5g} с'
        context=('time_' if time_axis else f't{data["time_s"][ti]:.7g}_')+('avg' if o.cell<0 else 'cell'+str(o.cell+1))
        source_note=''
        if not time_axis and o.mark_starts:
            starts=[f"{row['id']}: {row['start_s']:.4g} с" for row in data.meta.get('sources',[]) if row.get('normal_power_W_cm2',row.get('rate_cm2_s',0.))>0]
            if starts:source_note='\n'.join(textwrap.wrap('Начало накачек · '+'; '.join(starts),105))
        for part in range(0,max(len(curves),1),o.curves_per_figure):
            chunk=curves[part:part+o.curves_per_figure];wrap_width=max(35,int(920/o.fontsize));legend_height=.025*o.fontsize*sum(max(1,len(textwrap.wrap(label,wrap_width))) for label,_ in chunk)
            fig=Figure(figsize=(11,6.1+legend_height),dpi=125);FigureCanvasAgg(fig);ax=fig.add_subplot(111)
            lower=(1.+legend_height)/(6.1+legend_height);fig.subplots_adjust(left=.14,right=.97,top=.84 if source_note else .90,bottom=lower)
            if source_note:fig.text(.14,.985,source_note,ha='left',va='top',fontsize=max(o.fontsize-4,11),color='#555555')
            if o.x_scale!='linear':ax.set_xscale(o.x_scale)
            if o.y_scale!='linear':ax.set_yscale(o.y_scale)
            any_line=False;handles=[]
            for label,values in chunk:
                yy=np.asarray(values,dtype=float).copy();mask=np.ones(len(xx),dtype=bool)
                if o.x_scale=='log':mask&=xx>0
                if o.y_scale=='log':yy[yy<=0]=np.nan
                visible=np.isfinite(yy[mask]).any()
                caption=label+(' · нет положительных значений' if o.y_scale=='log' and not visible else '')
                # A labelled handle is retained for zero curves on a log axis.
                line,=ax.plot(xx[mask],yy[mask],lw=2.1,label='\n'.join(textwrap.wrap(caption,wrap_width)) or caption,marker='o' if len(xx)<12 else None,markersize=3.5)
                handles.append(line);any_line=any_line or visible
            if o.y_scale=='log':
                positive=[np.asarray(values,dtype=float).ravel() for _,values in chunk]
                positive=np.concatenate(positive) if positive else np.array([])
                positive=positive[np.isfinite(positive)&(positive>0)]
                if len(positive) and np.ptp(positive)<=1e-10*positive.max():
                    # Keep a constant concentration readable on log Y; data
                    # and zero values are unchanged.
                    ax.set_ylim(positive.min()*.8,positive.max()*1.25)
            if not any_line:ax.text(.5,.5,empty_message,transform=ax.transAxes,ha='center',va='center',fontsize=o.fontsize,wrap=True)
            ax.set_xlabel(xlab,fontsize=o.fontsize);ax.set_ylabel(ylabel,fontsize=o.fontsize);ax.set_title(title,fontsize=o.fontsize+1,pad=16)
            ax.tick_params(labelsize=o.fontsize-1);ax.xaxis.get_offset_text().set_fontsize(o.fontsize-2);ax.yaxis.get_offset_text().set_fontsize(o.fontsize-2);ax.grid(alpha=.18);ax.margins(x=.03)
            if time_axis:
                domain=data['time_s'];positive=domain[domain>0]
                left=float(positive[0]) if o.x_scale=='log' and len(positive) else float(domain[0]);right=float(domain[-1])
                if left<right:ax.set_xlim(left,right)
                if o.time_min is not None:ax.set_xlim(left=o.time_min)
                if o.time_max is not None:ax.set_xlim(right=o.time_max)
                start_markers(ax,data.meta.get('sources',[]),o)
            if handles:
                # Anchor the legend INSIDE the figure, measure its real text,
                # then reserve that space. Axes-relative anchors clipped long
                # reaction names both in the Tk canvas and exported PNGs.
                legend=fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.14,.025),bbox_transform=fig.transFigure,borderaxespad=0,frameon=False,fontsize=max(o.fontsize-2,12),ncol=1)
                fig.canvas.draw();height=legend.get_window_extent().height/fig.dpi
                figure_height=max(6.1+legend_height,height+5.5)
                fig.set_size_inches(11,figure_height,forward=False)
                fig.subplots_adjust(bottom=(height+.95)/figure_height)
            pages.append(PlotPage(key+'_'+context+'_'+str(part//o.curves_per_figure+1),title,fig,xx,chunk,tuple(fig.get_size_inches())))
    for metric in o.metrics:
        if metric in ('levels','components','electrons'):
            curves=[]
            if metric=='levels':curves=[(curve_label(data.states[j]),context(n[:,:,j])) for j in level_ids]
            elif metric=='components':
                for name in components:
                    if name=='e⁻':curves.append(('e⁻',context(data['electrons_cm3'].sum(axis=(2,3)))))
                    elif name in data.meta['components']:curves.append((name,context(n[:,:,data.meta['components'][name]].sum(axis=2))))
            else:curves=[('e⁻',context(data['electrons_cm3'].sum(axis=(2,3))))]
            if o.fractions:
                denominator=context(n.sum(axis=2));curves=[(label,np.divide(v,denominator,out=np.zeros_like(v),where=denominator>0)) for label,v in curves]
            add(metric,LABELS[metric],'Доля от числа тяжёлых частиц' if o.fractions else 'Концентрация, см⁻³',curves)
        elif metric in ('mass','charge','energy','kinetic','kinetic_e'):
            k={'mass':0,'charge':1,'kinetic':2,'kinetic_e':2,'energy':3}[metric]
            if metric=='kinetic':
                values=data['kinetic_residual_cm3'];curves=[(curve_label(data.states[j]),np.abs(reduce_cell(values,o.cell)[:,j])) for j in level_ids];ylabel='|Невязка кинетики|, см⁻³'
            elif metric=='kinetic_e':curves=[('e⁻',abs(reduce_cell(data['kinetic_residual_e_cm3'],o.cell)))];ylabel='|Невязка кинетики e⁻|, см⁻³'
            else:curves=[(LABELS[metric],abs(data['conservation_relative'][:,k]))];ylabel='Модуль относительной невязки'
            add(metric,LABELS[metric],ylabel,curves,xx=data['time_s'],xlab='Время, с',time_axis=True)
        elif metric in ('rates','integrals'):
            values=context(data['reaction_rates' if metric=='rates' else 'reaction_integrals'])
            add(metric,LABELS[metric],'Скорость событий, см⁻³ с⁻¹' if metric=='rates' else 'Интеграл событий, см⁻³',[(data.reactions[k]['reaction']+' · '+data.reactions[k]['id'],values[:,k]) for k in selected])
        elif metric=='top5':
            target='Все реакции' if o.target=='all' else 'e⁻' if o.target=='electron' else (curve_label(data.states[states[o.target[6:]]]) if o.target.startswith('level:') else o.target.split(':',1)[1])
            directions=('source','sink') if o.direction=='both' and o.target!='all' else (o.direction,)
            for direction in directions:
                curves=[(data.reactions[k]['reaction']+' · '+data.reactions[k]['id'],value) for k,value in top_history(data,replace(o,direction=direction))]
                suffix='' if o.target=='all' else ' · '+{'source':'заселение','sink':'расселение','participation':'участие'}.get(direction,'события')
                key=metric+'_'+direction if len(directions)>1 else metric
                add(key,'Top-5 · '+target+suffix,'Скорость для цели, см⁻³ с⁻¹',curves,xx=data['time_s'],xlab='Время, с',time_axis=True,empty_message='Нет ненулевых реакций для этой цели\nв сохранённом расчёте')
        elif metric=='eedf':
            f=data['electrons_cm3'].sum(axis=3);f=f.mean(axis=1) if o.cell<0 else f[:,o.cell]
            curves=[]
            for i in sorted({0,ti}):
                total=float(f[i].sum());curves.append((f't = {data["time_s"][i]:.5g} с',f[i]/max(total,1e-300)/data['energy_weights_eV']))
            add(metric,LABELS[metric],'Вероятность, эВ⁻¹',curves,xx=data['energy_eV'],xlab='Энергия электронов, эВ',time_axis=False)
        elif metric=='iedf':
            h=data['heavy_cm3'][ti];h=h.mean(axis=0) if o.cell<0 else h[o.cell]
            from radiation import weights
            width=weights(data['heavy_energy_eV']);curves=[]
            for j in level_ids:
                marginal=h[j].sum(axis=1);total=float(marginal.sum())
                curves.append((curve_label(data.states[j])+f' · t={data["time_s"][ti]:.5g} с',marginal/max(total,1e-300)/width))
            add(metric,LABELS[metric],'Плотность вероятности, эВ⁻¹',curves,xx=data['heavy_energy_eV'],xlab='Кинетическая энергия частицы, эВ',time_axis=False)
        elif metric in ('heavy_energy','velocity'):
            curves=[];energy=data['heavy_energy_eV'][:,None];mu=data['heavy_mu'][None,:];h=data['heavy_cm3']
            for j in level_ids:
                weight=energy if metric=='heavy_energy' else np.sqrt(2*energy*1.602176634e-12/(float(data.states[j]['mass_u'])*1.66053906892e-24))*mu
                numerator=context((h[:,:,j]*weight).sum(axis=(2,3)));den=context(n[:,:,j]);curves.append((curve_label(data.states[j]),np.divide(numerator,den,out=np.zeros_like(numerator),where=den>0)))
            add(metric,LABELS[metric],'Средняя энергия, эВ' if metric=='heavy_energy' else 'Средняя u_z, см/с',curves)
        elif metric=='heavy_distribution':
            h=data['heavy_cm3'][ti];h=h.mean(axis=0) if o.cell<0 else h[o.cell];curves=[];speed=np.linspace(-1,1,121)
            # One velocity axis per mass, avoiding an invalid shared E-to-v map.
            for j in level_ids:
                v=np.sqrt(2*data['heavy_energy_eV'][:,None]*1.602176634e-12/(float(data.states[j]['mass_u'])*1.66053906892e-24))*data['heavy_mu'][None,:]
                occupied=h[j].ravel()>max(float(h[j].sum())*1e-12,0.);vmax=max(float(abs(v.ravel()[occupied]).max()) if occupied.any() else 1.,1.)
                bins=np.linspace(-vmax,vmax,81);counts,edges=np.histogram(v.ravel(),bins,weights=h[j].ravel());density=counts/np.diff(edges)
                add(metric+'_'+str(j),curve_label(data.states[j])+f' · t={data["time_s"][ti]:.4g} с','dn/dv_z, см⁻³ / (см/с)',[(curve_label(data.states[j]),density)],xx=(edges[:-1]+edges[1:])/2,xlab='Скорость v_z, см/с',time_axis=False)
        elif metric in ('absorption','transmission'):
            if 'probe_wavelength_nm' not in data.arrays:continue
            wave=data['probe_wavelength_nm'];mask=np.ones(len(wave),dtype=bool)
            if o.wavelength_min is not None:mask&=wave>=o.wavelength_min
            if o.wavelength_max is not None:mask&=wave<=o.wavelength_max
            if not mask.any():raise ValueError('Выбранный диапазон не пересекает пробный спектр.')
            values=data['probe_transmission'][ti];values=100*(1-values) if metric=='absorption' else values
            add(metric,LABELS[metric]+f' · t={data["time_s"][ti]:.4g} с','Поглощение, %' if metric=='absorption' else 'I / I₀',[(LABELS[metric],values[mask])],xx=wave[mask],xlab='Длина волны в вакууме, нм',time_axis=False)
        elif metric=='probe_distance':
            if 'probe_wavelength_nm' not in data.arrays:continue
            k=int(np.argmin(abs(data['probe_wavelength_nm']-o.probe_wavelength_nm)));tr=np.exp(-data['probe_tau_position'][ti,:,k]);wave=data['probe_wavelength_nm'][k]
            add(metric,f'Пропускание · λ={wave:.7g} нм','I(x) / I₀',[(f'λ={wave:.7g} нм',tr)],xx=data['position_cm'],xlab='Расстояние, см',time_axis=False)
        elif metric=='emission':
            wave=2.99792458e17/data['photon_nu_Hz'];order=np.argsort(wave);wave=wave[order];mask=np.ones(len(wave),dtype=bool)
            if o.wavelength_min is not None:mask&=wave>=o.wavelength_min
            if o.wavelength_max is not None:mask&=wave<=o.wavelength_max
            if not mask.any():raise ValueError('Выбранный диапазон не пересекает расчётный спектр излучения.')
            # F_nu [erg cm^-2 s^-1 Hz^-1] -> F_lambda [W cm^-2 nm^-1].
            jac=2.99792458e10/(wave*1e-7)**2*1e-14
            curves=[('Выход влево',data['outgoing_left'][ti][order][mask]*jac[mask]),('Выход вправо',data['outgoing_right'][ti][order][mask]*jac[mask])]
            add(metric,LABELS[metric]+f' · t={data["time_s"][ti]:.4g} с','Спектральный поток, Вт см⁻² нм⁻¹',curves,xx=wave[mask],xlab='Длина волны в вакууме, нм',time_axis=False)
    return pages


def save_pages(pages,folder,export_csv=False):
    from pathlib import Path
    import csv
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);files=[]
    for page in pages:
        path=folder/(page.key+'.png');copy=2
        while path.exists():path=folder/(page.key+'_'+str(copy)+'.png');copy+=1
        previous=tuple(page.figure.get_size_inches())
        page.figure.set_size_inches(page.size_inches or previous,forward=False);page.figure.savefig(path,dpi=200,facecolor='white');page.figure.set_size_inches(previous,forward=False);files.append(path)
        if export_csv:
            fields=['X']+[label for label,_ in page.curves]
            with path.with_suffix('.csv').open('w',encoding='utf-8-sig',newline='') as stream:
                writer=csv.writer(stream,delimiter=';');writer.writerow(fields)
                for i,x in enumerate(page.x):writer.writerow([float(x)]+[float(values[i]) if np.isfinite(values[i]) else '' for _,values in page.curves])
    return files
