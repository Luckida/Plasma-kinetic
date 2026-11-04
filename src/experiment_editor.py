"""Gas, weak-probe and optional pump editor. No numerical dependencies."""
from copy import deepcopy
import json
import os
from pathlib import Path
import tkinter as tk
from tkinter import ttk,filedialog,messagebox
from experiment_config import editor_composition,finite,gas_initial,prepare_configuration
from spectrum_input import load_spectrum
from light_input import light_spectrum,light_enabled,integral

MODES={'Концентрации, см⁻³':'densities','Доли, % + полная концентрация':'fractions','Доли, % + давление заполнения':'pressure'}


class ExperimentEditor:
    def __init__(self,parent,path):
        self.path=Path(path).resolve();self.config=json.loads(self.path.read_text(encoding='utf-8-sig'))
        self.window=tk.Toplevel(parent);w=self.window;w.title('Состав, распределения и спектры · 2.8.1');w.geometry('1180x880');w.minsize(1000,760)
        self.vars={};self.widgets={};self.selected=None;self.last_mode=editor_composition(self.config)['mode']
        outer=ttk.Frame(w,padding=12);outer.pack(fill='both',expand=True)
        outer.columnconfigure(0,weight=1);outer.rowconfigure(1,weight=1)
        ttk.Label(outer,text='Состав, пробный свет и расчёт',font=('Segoe UI',16,'bold')).grid(row=0,column=0,sticky='w')
        book=ttk.Notebook(outer);book.grid(row=1,column=0,sticky='nsew',pady=10)
        gas=ttk.Frame(book,padding=12);probe=ttk.Frame(book,padding=12);pump=ttk.Frame(book,padding=12);run=ttk.Frame(book,padding=12);sources=ttk.Frame(book);distributions=ttk.Frame(book)
        self.book=book;self.distributions_tab=distributions
        for tab,label in [(gas,'Состав газа'),(distributions,'EEDF / ионы'),(probe,'Пробный спектр'),(pump,'Свет на границе'),(sources,'Источники в слое'),(run,'Расчёт')]:book.add(tab,text=label)
        composition=editor_composition(self.config)
        self.var('gas_mode',next(label for label,value in MODES.items() if value==composition['mode']))
        ttk.Label(gas,text='Как задать состав:').grid(row=0,column=0,sticky='w')
        choice=ttk.Combobox(gas,textvariable=self.vars['gas_mode'],values=list(MODES),state='readonly',width=42);choice.grid(row=0,column=1,columnspan=3,sticky='w');choice.bind('<<ComboboxSelected>>',self.change_mode)
        self.entry(gas,1,'Температура газа, K','temperature',composition.get('temperature_K',300.))
        self.entry(gas,2,'Полная концентрация, см⁻³','total',composition.get('total_density_cm3',sum(self.config.get('initial_cm3',{}).values()) or 1e10))
        self.entry(gas,3,'Давление заполнения, Па','pressure',composition.get('pressure_Pa',1.))
        self.tree=ttk.Treeview(gas,columns=('species','value'),show='headings',height=9,selectmode='browse')
        self.tree.heading('species',text='Компонент / уровень');self.tree.heading('value',text='Концентрация, см⁻³' if self.last_mode=='densities' else 'Доля, %')
        self.tree.column('species',width=370);self.tree.column('value',width=310);self.tree.grid(row=4,column=0,columnspan=4,sticky='nsew',pady=8)
        for row in composition['components']:
            value=row.get('density_cm3',0.) if self.last_mode=='densities' else 100*row.get('fraction',0.)
            self.tree.insert('','end',values=(row['species'],f'{value:.10g}'))
        self.tree.bind('<<TreeviewSelect>>',self.select_row)
        self.var('species','');self.var('amount','')
        from particle_choices import particle_choices
        self.available_particles=particle_choices(self.path.parent,self.config)
        self.species_box=ttk.Combobox(gas,textvariable=self.vars['species'],values=self.available_particles,width=25);self.species_box.grid(row=5,column=0,sticky='ew',padx=(0,6))
        ttk.Entry(gas,textvariable=self.vars['amount'],width=24).grid(row=5,column=1,sticky='ew')
        ttk.Button(gas,text='Добавить / применить',command=self.apply_row_dialog).grid(row=5,column=2,padx=6)
        ttk.Button(gas,text='Удалить',command=self.remove_row).grid(row=5,column=3)
        ttk.Button(gas,text='Нормировать до 100%',command=self.normalize_dialog).grid(row=6,column=1,columnspan=2,pady=8,sticky='w')
        ttk.Label(gas,text='Дважды вводить ионы не нужно: He+ — исходный ион, HeI_2_3S — метастабильный уровень.\nДоли относятся к числу тяжёлых частиц. Давление — давление заполнения, n = p/(kT).',wraplength=880).grid(row=7,column=0,columnspan=4,sticky='w',pady=6)
        for c in (0,1):gas.columnconfigure(c,weight=1)
        gas.rowconfigure(4,weight=1)
        e=self.config.get('electron_initial',{})
        self.entry(run,0,'Температура электронов, эВ','Te',e.get('Te_eV',.8))
        from distribution_editor import InitialDistributions
        self.distributions_editor=InitialDistributions(distributions,self.path.parent,self.config,self.vars['Te'],self.vars['temperature']);self.distributions_editor.pack(fill='both',expand=True);self.distributions_editor.set_targets([name.split(' · ',1)[0] for name in self.available_particles]+[row['species'] for row in composition['components']])
        self.var('ne_auto','density_cm3' not in e,boolean=True)
        ttk.Checkbutton(run,text='Исходная концентрация электронов — по заряду ионов',variable=self.vars['ne_auto']).grid(row=1,column=0,columnspan=3,sticky='w',pady=7)
        self.entry(run,2,'Или заданная nₑ, см⁻³','ne',e.get('density_cm3',1e7))
        solver=self.config.get('solver',{})
        for row,label,key,value in [(3,'Время расчёта, с','end',solver.get('t_end_s',2e-8)),(4,'Максимальный шаг, с','step',solver.get('max_step_s',1e-9)),(5,'Число сохраняемых моментов','samples',solver.get('samples',5)),(6,'Число ячеек по расстоянию','cells',self.config.get('cells',2)),(7,'Толщина слоя, см','length',self.config.get('thickness_cm',.01))]:self.entry(run,row,label,key,value)
        self.var('simple',True,boolean=True)
        ttk.Checkbutton(run,text='Обычный режим: без внешнего поля и неявных горячих добавок',variable=self.vars['simple']).grid(row=8,column=0,columnspan=3,sticky='w',pady=10)
        ttk.Label(run,text='Обычный режим также отключает оценочные сечения и BGK-релаксацию.\nСтолкновительная кинетика по атомным данным и перенос излучения сохраняются.\nНачальный газ и накачки задаются отдельно. Источники — на вкладке «Источники в слое».\nЕсли исходный газ полностью нейтрален, режим «по заряду ионов» даёт nₑ=0.\nТемпература EEDF не создаёт электроны: их число задаётся явно, ионами или источником.',wraplength=850).grid(row=9,column=0,columnspan=3,sticky='w',pady=10)
        pcfg=self.config.get('probe',{})
        self.var('probe_mode',pcfg.get('mode','flat'))
        ttk.Label(probe,text='Режим пробного света:').grid(row=0,column=0,sticky='w')
        ttk.Combobox(probe,textvariable=self.vars['probe_mode'],values=('flat','file','off'),state='readonly',width=22).grid(row=0,column=1,sticky='w')
        ttk.Label(probe,text='flat — ровный пробный спектр; file — файл; off — без зонда.\nСлабый зонд показывает пропускание и не меняет концентрации.',wraplength=850).grid(row=1,column=0,columnspan=3,sticky='w',pady=10)
        self.entry(probe,2,'Начало ровного спектра, нм','probe_min',pcfg.get('minimum_nm',1082.))
        self.entry(probe,3,'Конец ровного спектра, нм','probe_max',pcfg.get('maximum_nm',1084.))
        self.entry(probe,4,'Точек исходного спектра','probe_points',pcfg.get('points',501))
        self.spectrum_controls(probe,5,'probe',pcfg)
        ttk.Label(probe,text='Поддержаны TXT, CSV, TSV, DAT, ASC и XLSX. Заголовки произвольные или отсутствуют.\nВ выходных файлах: входной/прошедший свет, пропускание и поглощение по λ, времени и расстоянию.\nДля пробного света используются связанные линии; континуум и примесь спонтанного света не добавляются.\nДлины волн модели — в вакууме, из энергий LS; разрешённая структура линий ограничена банком данных.',wraplength=875).grid(row=12,column=0,columnspan=3,sticky='w',pady=12)
        rcfg=self.config.get('radiation',{})
        self.var('pump_enabled',light_enabled(rcfg),boolean=True)
        ttk.Checkbutton(pump,text='Включить световую накачку на входе слоя',variable=self.vars['pump_enabled']).grid(row=0,column=0,columnspan=3,sticky='w',pady=8)
        self.var('pump_use_probe',rcfg.get('spectrum_mode')=='probe' or not rcfg.get('spectrum_file'),boolean=True)
        ttk.Checkbutton(pump,text='Использовать форму из вкладки «Пробный спектр»',variable=self.vars['pump_use_probe']).grid(row=1,column=0,columnspan=3,sticky='w',pady=8)
        self.entry(pump,2,'Или отдельный файл формы','pump_file',rcfg.get('spectrum_file',''))
        button=ttk.Button(pump,text='Выбрать…',command=lambda:self.browse('pump'));button.grid(row=2,column=2,padx=8);self.widgets['pump_browse']=button
        opts=rcfg.get('spectrum_options',{})
        ttk.Label(pump,text='Единицы X файла').grid(row=3,column=0,sticky='w')
        widget=ttk.Combobox(pump,textvariable=self.var('pump_x_unit',opts.get('x_unit','auto')),values=('auto','nm','A','um','eV','Hz'),state='readonly',width=25);widget.grid(row=3,column=1,sticky='w',pady=6);self.widgets['pump_x_unit']=widget
        for key,value in [('y_kind','relative'),('energy_unit','cgs'),('x_column',int(opts.get('x_column',0))+1),('y_column',int(opts.get('y_column',1))+1)]:self.var('pump_'+key,value)
        self.var('pump_angular','normal_flux')
        power=rcfg.get('multiplier_W_cm2',1e-4)
        if light_enabled(rcfg) and rcfg.get('normalization')!='integral':
            old=light_spectrum(rcfg,self.path.parent);power=integral(old.wavelength_A,old.values)*float(rcfg.get('beam_mu',1.))*float(rcfg.get('scale',1.))*1e-7
        self.entry(pump,4,'Множитель: суммарная мощность, Вт/см²','pump_scale',power)
        self.entry(pump,5,'Начало накачки, с','pump_start',rcfg.get('start_s',0.))
        self.entry(pump,6,'Конец накачки, с','pump_stop',rcfg.get('stop_s',solver.get('t_end_s',2e-8)))
        button=ttk.Button(pump,text='Проверить форму и показать энергетику',command=self.pump_energy);button.grid(row=7,column=0,columnspan=3,sticky='w',pady=12);self.widgets['pump_preview']=button
        self.pump_info=tk.StringVar(value='Форма нормируется по интегралу до 1. Множитель задаёт мощность через плоскость входа.')
        ttk.Label(pump,textvariable=self.pump_info,wraplength=960,font=('Segoe UI',11)).grid(row=8,column=0,columnspan=3,sticky='w',pady=12)
        ttk.Label(pump,text='Можно использовать относительные отсчёты любого поддерживаемого файла.\nЭнергия каждого фотона определяется его длиной волны: E = hc/λ.\nДиапазон и средняя энергия фотонов, мощность и энергия импульса показываются после проверки.\nСвет воздействует только в заданном интервале. Пробная диагностика сама не меняет кинетику.',wraplength=960).grid(row=9,column=0,columnspan=3,sticky='w',pady=12)
        from source_editor import SourceEditor
        self.sources_editor=SourceEditor(sources,self.config,self.path.parent,self.current_probe);self.sources_editor.pack(fill='both',expand=True)
        self.sources_editor.geometry_provider=lambda:(finite(self.vars['length'].get(),'Толщина'),finite(self.vars['end'].get(),'Время'))
        footer=ttk.Frame(outer);footer.grid(row=2,column=0,sticky='ew')
        self.save_button=ttk.Button(footer,text='Сохранить',command=self.save_dialog);self.save_button.pack(side='right')
        ttk.Button(footer,text='Закрыть',command=w.destroy).pack(side='right',padx=8)
        self.status=tk.StringVar(value='Изменения вступают в силу после сохранения.');ttk.Label(footer,textvariable=self.status,wraplength=670).pack(side='left')
        for key in ('gas_mode','probe_mode','pump_enabled','pump_use_probe','ne_auto'):self.vars[key].trace_add('write',self.update_states)
        self.update_states()
        self.pump_energy(quiet=True)

    def var(self,key,value,boolean=False):
        self.vars[key]=(tk.BooleanVar if boolean else tk.StringVar)(master=self.window,value=value);return self.vars[key]

    def entry(self,parent,row,label,key,value):
        ttk.Label(parent,text=label).grid(row=row,column=0,sticky='w',pady=5,padx=(0,10))
        widget=ttk.Entry(parent,textvariable=self.var(key,value),width=28);widget.grid(row=row,column=1,sticky='ew',pady=5);self.widgets[key]=widget

    def spectrum_controls(self,parent,row,prefix,cfg):
        self.entry(parent,row,'Файл спектра',prefix+'_file',cfg.get('spectrum_file',''))
        widget=ttk.Button(parent,text='Выбрать…',command=lambda:self.browse(prefix));widget.grid(row=row,column=2,padx=6);self.widgets[prefix+'_browse']=widget
        options=cfg.get('spectrum_options',{})
        for offset,label,key,values,default in [(1,'Единицы X','x_unit',('auto','nm','A','um','eV','Hz'),'auto'),(2,'Интенсивность','y_kind',('relative','energy','photons'),'relative' if prefix=='probe' else 'energy'),(3,'Единицы мощности / геометрия','energy_unit',('cgs','W_m2'),'cgs')]:
            ttk.Label(parent,text=label).grid(row=row+offset,column=0,sticky='w',pady=4)
            widget=ttk.Combobox(parent,textvariable=self.var(prefix+'_'+key,options.get(key,default)),values=values,state='readonly',width=25);widget.grid(row=row+offset,column=1,sticky='w');self.widgets[prefix+'_'+key]=widget
        if prefix=='pump':
            widget=ttk.Combobox(parent,textvariable=self.var('pump_angular',options.get('angular_kind','four_pi_J')),values=('four_pi_J','normal_flux'),state='readonly',width=22);widget.grid(row=row+3,column=2,sticky='w');self.widgets['pump_angular']=widget
        cols=ttk.Frame(parent);cols.grid(row=row+4,column=0,columnspan=3,sticky='w',pady=6)
        for key,label,default in [('x_column','Столбец X (с 1)',0),('y_column','Столбец Y',1)]:
            ttk.Label(cols,text=label).pack(side='left',padx=(0,5));widget=ttk.Entry(cols,textvariable=self.var(prefix+'_'+key,int(options.get(key,default))+1),width=6);widget.pack(side='left',padx=(0,14));self.widgets[prefix+'_'+key]=widget
        widget=ttk.Button(parent,text='Проверить файл и показать точки',command=lambda:self.preview_dialog(prefix));widget.grid(row=row+5,column=0,columnspan=3,sticky='w',pady=8);self.widgets[prefix+'_preview']=widget

    def update_states(self,*_):
        gas_mode=MODES[self.vars['gas_mode'].get()];probe_mode=self.vars['probe_mode'].get();pump=bool(self.vars['pump_enabled'].get())
        states={'total':gas_mode=='fractions','pressure':gas_mode=='pressure','ne':not self.vars['ne_auto'].get()}
        for key in self.widgets:
            if key.startswith('pump_'):states[key]=pump and (not self.vars['pump_use_probe'].get() if key in ('pump_file','pump_browse','pump_x_unit') else True)
            elif key.startswith('probe_'):states[key]=(probe_mode=='flat' if key in ('probe_min','probe_max','probe_points') else probe_mode=='file')
        for key,enabled in states.items():
            widget=self.widgets[key];widget.configure(state=('readonly' if widget.winfo_class()=='TCombobox' else 'normal') if enabled else 'disabled')

    def browse(self,prefix):
        file=filedialog.askopenfilename(parent=self.window,filetypes=[('Спектры','*.txt *.csv *.tsv *.dat *.asc *.xlsx'),('Все файлы','*')])
        if file:
            try:file=os.path.relpath(file,self.path.parent)
            except ValueError:pass
            self.vars[prefix+'_file'].set(file)
            if prefix=='probe':self.vars['probe_mode'].set('file')
            else:self.vars['pump_use_probe'].set(False)

    def select_row(self,event=None):
        selected=self.tree.selection()
        if selected:
            self.selected=selected[0];name,value=self.tree.item(self.selected,'values');self.vars['species'].set(name);self.vars['amount'].set(value)

    def apply_row(self):
        from particle_choices import input_name
        name=input_name(self.vars['species'].get())
        if not name:raise ValueError('Укажите компонент, например He или H+.')
        value=finite(self.vars['amount'].get(),'Концентрация/доля')
        if value<0:raise ValueError('Количество не может быть отрицательным.')
        target=self.selected
        for item in self.tree.get_children():
            if self.tree.item(item,'values')[0]==name:target=item
        if target and self.tree.exists(target):self.tree.item(target,values=(name,f'{value:.12g}'))
        else:self.tree.insert('','end',values=(name,f'{value:.12g}'))
        self.selected=None;self.tree.selection_remove(*self.tree.selection());self.vars['species'].set('');self.vars['amount'].set('')

    def apply_row_dialog(self):
        try:self.apply_row()
        except ValueError as exc:messagebox.showerror('Состав',str(exc),parent=self.window)

    def remove_row(self):
        for item in self.tree.selection():self.tree.delete(item)
        self.selected=None;self.vars['species'].set('');self.vars['amount'].set('')

    def normalize(self):
        if MODES[self.vars['gas_mode'].get()]=='densities':raise ValueError('Выберите режим долей; концентрации нормировать не нужно.')
        items=self.tree.get_children();total=sum(finite(self.tree.item(item,'values')[1],'Доля') for item in items)
        if total<=0:raise ValueError('Нужна положительная сумма долей.')
        for item in items:
            name,value=self.tree.item(item,'values');self.tree.item(item,values=(name,f'{100*finite(value,"Доля")/total:.12g}'))

    def normalize_dialog(self):
        try:self.normalize()
        except ValueError as exc:messagebox.showerror('Доли',str(exc),parent=self.window)

    def change_mode(self,event=None):
        mode=MODES[self.vars['gas_mode'].get()]
        if self.vars['species'].get().strip():self.apply_row_dialog()
        items=self.tree.get_children();total=sum(finite(self.tree.item(item,'values')[1],'Количество') for item in items)
        if self.last_mode=='densities' and mode!='densities' and total>0:
            self.vars['total'].set(f'{total:.12g}')
            for item in items:
                name,value=self.tree.item(item,'values');self.tree.item(item,values=(name,f'{100*finite(value,"Количество")/total:.12g}'))
        elif self.last_mode!='densities' and mode=='densities':
            density=finite(self.vars['total'].get(),'Полная концентрация') if self.last_mode=='fractions' else finite(self.vars['pressure'].get(),'Давление')/(1.380649e-23*finite(self.vars['temperature'].get(),'Температура'))*1e-6
            for item in items:
                name,value=self.tree.item(item,'values');self.tree.item(item,values=(name,f'{finite(value,"Доля")/100*density:.12g}'))
        self.tree.heading('value',text='Концентрация, см⁻³' if mode=='densities' else 'Доля, %');self.last_mode=mode

    def options(self,prefix):
        out={key:self.vars[prefix+'_'+key].get() for key in ('x_unit','y_kind','energy_unit')}
        out.update({key:int(self.vars[prefix+'_'+key].get())-1 for key in ('x_column','y_column')})
        if prefix=='pump':out.update(angular_kind=self.vars['pump_angular'].get(),mu=float(self.config.get('radiation',{}).get('beam_mu',1.)))
        return out

    def preview(self,prefix):
        return load_spectrum(self.path.parent/self.vars[prefix+'_file'].get(),self.options(prefix),purpose='driving' if prefix=='pump' else 'probe')

    def preview_dialog(self,prefix):
        try:
            data=self.preview(prefix);panel=tk.Toplevel(self.window);panel.title('Предпросмотр импорта');panel.geometry('780x440')
            meta=data.metadata
            ttk.Label(panel,text=f"Точек: {meta['points']} · {meta['wavelength_min_nm']:.8g}…{meta['wavelength_max_nm']:.8g} нм\nВход: {meta['input_x_unit']}, {meta['input_y_kind']}\n"+'\n'.join(meta['warnings']),wraplength=740).pack(anchor='w',padx=12,pady=10)
            table=ttk.Treeview(panel,columns=('x','y'),show='headings',height=11);table.heading('x',text='Длина волны после перевода, нм');table.heading('y',text='Интенсивность после перевода');table.pack(fill='both',expand=True,padx=12,pady=10)
            for x,y in list(zip(data.wavelength_A,data.values))[:50]:table.insert('','end',values=(f'{x/10:.10g}',f'{y:.10g}'))
        except (OSError,ValueError,KeyError) as exc:messagebox.showerror('Спектр',str(exc),parent=self.window)

    def current_probe(self):
        mode=self.vars['probe_mode'].get();out={'mode':mode}
        if mode=='file':out.update(spectrum_file=self.vars['probe_file'].get(),spectrum_options=self.options('probe'))
        elif mode=='flat':out.update(minimum_nm=finite(self.vars['probe_min'].get(),'Начало спектра'),maximum_nm=finite(self.vars['probe_max'].get(),'Конец спектра'),points=int(self.vars['probe_points'].get()))
        return out

    def pump_config(self):
        return dict(normalization='integral',spectrum_mode='probe' if self.vars['pump_use_probe'].get() else 'file',spectrum_file=self.vars['pump_file'].get(),spectrum_options={'x_unit':self.vars['pump_x_unit'].get(),'y_kind':'relative'},multiplier_W_cm2=finite(self.vars['pump_scale'].get(),'Множитель света'),beam_mu=float(self.config.get('radiation',{}).get('beam_mu',1.)))

    def pump_energy(self,quiet=False):
        try:
            info=light_spectrum(self.pump_config(),self.path.parent,self.current_probe()).metadata
            duration=finite(self.vars['pump_stop'].get(),'Конец')-finite(self.vars['pump_start'].get(),'Начало')
            self.pump_info.set(f"Фотоны: {info['minimum_photon_eV']:.7g}…{info['maximum_photon_eV']:.7g} эВ. Средняя энергия: {info['mean_photon_eV']:.7g} эВ.\nМощность: {info['normal_power_W_cm2']:.7g} Вт/см²; энергия импульса: {duration*info['normal_power_W_cm2']:.7g} Дж/см².\nФотонный поток: {info['photon_flux_cm2_s']:.7g} см⁻² с⁻¹.")
        except (ValueError,OSError) as exc:
            if quiet:self.pump_info.set('Проверьте форму спектра: '+str(exc))
            else:messagebox.showerror('Энергетика спектра',str(exc),parent=self.window)

    def build_config(self):
        if self.vars['species'].get().strip():self.apply_row()
        cfg=json.loads(self.path.read_text(encoding='utf-8-sig'));mode=MODES[self.vars['gas_mode'].get()]
        gas=dict(mode=mode,temperature_K=finite(self.vars['temperature'].get(),'Температура газа'),components=[])
        for item in self.tree.get_children():
            name,value=self.tree.item(item,'values');gas['components'].append(dict(species=name,**{('density_cm3' if mode=='densities' else 'fraction'):finite(value,'Количество '+name)/(1 if mode=='densities' else 100)}))
        if mode=='fractions':gas['total_density_cm3']=finite(self.vars['total'].get(),'Полная концентрация')
        if mode=='pressure':gas['pressure_Pa']=finite(self.vars['pressure'].get(),'Давление заполнения')
        gas_initial(gas,self.path.parent,cfg.get('chemistry'));cfg['gas']=gas;cfg.pop('cell_overrides',None)
        electron,heavy=self.distributions_editor.build();cfg['gas']['distributions']=heavy
        cfg['phase_grid'],cfg['heavy_phase_grid']=self.distributions_editor.grids()
        e=electron;cfg['electron_initial']=e
        if self.vars['ne_auto'].get():e.pop('density_cm3',None)
        else:
            e['density_cm3']=finite(self.vars['ne'].get(),'Концентрация электронов')
            if e['density_cm3']<0:raise ValueError('Концентрация электронов >=0.')
        end=finite(self.vars['end'].get(),'Время расчёта');step=finite(self.vars['step'].get(),'Шаг');length=finite(self.vars['length'].get(),'Толщина')
        cells=int(self.vars['cells'].get());samples=int(self.vars['samples'].get())
        if end<=0 or not 0<step<=end or length<=0 or cells<1 or samples<2:raise ValueError('0<шаг<=время расчёта; толщина >0; ячеек >=1, моментов >=2.')
        cfg.update(cells=cells,thickness_cm=length);cfg.setdefault('solver',{}).update(t_end_s=end,max_step_s=step,samples=samples)
        if self.vars['simple'].get():
            cfg['spatial_sources']=[];cfg['electric_field']={'mode':'prescribed','Ez_V_cm':0.}
            cfg.setdefault('chemistry',{})['allow_estimated_rates']=False
            cfg.setdefault('collisions',{}).update(ee_bgk=False,ii_bgk=False)
        mode=self.vars['probe_mode'].get();cfg['probe']={'mode':mode}
        if mode=='file':
            self.preview('probe');cfg['probe'].update(spectrum_file=self.vars['probe_file'].get(),spectrum_options=self.options('probe'))
        elif mode=='flat':
            low=finite(self.vars['probe_min'].get(),'Начало спектра');high=finite(self.vars['probe_max'].get(),'Конец спектра');points=int(self.vars['probe_points'].get())
            if not 0<low<high or not 2<=points<=200000:raise ValueError('Пробный спектр: 0<начало<конец, точек 2…200000.')
            cfg['probe'].update(minimum_nm=low,maximum_nm=high,points=points)
        cfg['spatial_sources']=self.sources_editor.build(length,end)
        r=cfg.setdefault('radiation',{});r['enabled']=bool(self.vars['pump_enabled'].get());r['scale']=0.
        if r['enabled']:
            light=self.pump_config();start=finite(self.vars['pump_start'].get(),'Начало накачки');stop=finite(self.vars['pump_stop'].get(),'Конец накачки')
            if not 0<=start<stop<=end:raise ValueError('Накачка: 0<=начало<конец<=время расчёта.')
            light_spectrum(light,self.path.parent,cfg['probe']);r.update(light,start_s=start,stop_s=stop,scale=1.)
            r.pop('history_file',None);r.update(left_diffuse_scale=0.,right_diffuse_scale=0.)
        cfg.setdefault('diagnostics',{}).update(output_mode='compact',plots_on_finish=False)
        return prepare_configuration(cfg,self.path.parent)

    def save(self):
        cfg=self.build_config();temporary=self.path.with_suffix(self.path.suffix+'.tmp');temporary.write_text(json.dumps(cfg,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');temporary.replace(self.path)
        self.config=cfg;self.status.set('Сохранено: '+self.path.name);return cfg

    def save_dialog(self):
        try:self.save()
        except (ValueError,OSError,KeyError) as exc:messagebox.showerror('Параметры',str(exc),parent=self.window)


def open_editor(parent,path):return ExperimentEditor(parent,path)
