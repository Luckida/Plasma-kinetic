"""The same distribution form for electrons, ions and neutral particles."""
from copy import deepcopy
from pathlib import Path
import os
import tkinter as tk
from tkinter import ttk,filedialog,messagebox
from experiment_config import finite

KINDS={'Максвелловское (по умолчанию)':'maxwellian','Файл P(E) / EEPF · изотропный':'energy_table','Фазовый CSV · P(E, μ)':'phase_table'}
FORMS={'P(E), эВ⁻¹ · ∫P dE = 1':'probability','EEPF, эВ⁻³ᐟ² · ∫√E f dE = 1':'eepf'}


class DistributionEditor(ttk.LabelFrame):
    def __init__(self,parent,title,root,options,temperature_var=None,electron=True,allow_energy=True):
        super().__init__(parent,text=title,padding=10);self.root_path=Path(root);self.electron=electron;self.allow_energy=allow_energy;self.base=deepcopy(options)
        self.kind=tk.StringVar(self);self.file=tk.StringVar(self);self.form=tk.StringVar(self)
        self.temperature=temperature_var or tk.StringVar(self,value='1' if electron else '300')
        ttk.Label(self,text='Форма распределения').grid(row=0,column=0,columnspan=3,sticky='w',pady=3)
        self.kind_box=ttk.Combobox(self,textvariable=self.kind,values=[label for label,value in KINDS.items() if allow_energy or value!='energy_table'],state='readonly',width=28);self.kind_box.grid(row=1,column=0,columnspan=3,sticky='ew',pady=3)
        self.temp_label=ttk.Label(self,text='Te, эВ' if electron else 'Ti, K');self.temp_label.grid(row=2,column=0,sticky='w',pady=3)
        self.temp_entry=ttk.Entry(self,textvariable=self.temperature,width=12);self.temp_entry.grid(row=2,column=1,columnspan=2,sticky='ew',pady=3)
        ttk.Label(self,text='Файл распределения').grid(row=3,column=0,columnspan=3,sticky='w',pady=3)
        self.file_entry=ttk.Entry(self,textvariable=self.file,width=16);self.file_entry.grid(row=4,column=0,columnspan=2,sticky='ew',pady=3)
        self.browse_button=ttk.Button(self,text='Выбрать…',command=self.browse);self.browse_button.grid(row=4,column=2,padx=5)
        self.form_box=ttk.Combobox(self,textvariable=self.form,values=list(FORMS),state='readonly',width=28);self.form_box.grid(row=5,column=0,columnspan=3,sticky='ew',pady=3)
        self.hint=tk.StringVar(self);ttk.Label(self,textvariable=self.hint,wraplength=410,justify='left').grid(row=6,column=0,columnspan=3,sticky='w',pady=5)
        self.columnconfigure(1,weight=1);self.kind.trace_add('write',self.update_state);self.load(options,keep_temperature=temperature_var is not None)

    def load(self,options,keep_temperature=False):
        self.base=deepcopy(options);kind=options.get('kind') or ('energy_table' if options.get('energy_distribution_file') else 'phase_table' if options.get('distribution_file') else 'maxwellian')
        self.kind.set(next(label for label,value in KINDS.items() if value==kind))
        self.file.set(options.get('energy_distribution_file',options.get('distribution_file','')))
        self.form.set(next(label for label,value in FORMS.items() if value==options.get('energy_form','probability')))
        if not keep_temperature:self.temperature.set(str(options.get('Te_eV' if self.electron else 'Ti_K',1. if self.electron else 300.)))
        self.update_state()

    def set_particle(self,electron):
        self.electron=electron;self.temp_label.configure(text='Te, эВ' if electron else 'Ti, K');self.update_state()

    def update_state(self,*_):
        kind=KINDS.get(self.kind.get(),'maxwellian');table=kind!='maxwellian'
        self.temp_entry.configure(state='disabled' if table else 'normal');self.file_entry.configure(state='normal' if table else 'disabled');self.browse_button.configure(state='normal' if table else 'disabled');self.form_box.configure(state='readonly' if kind=='energy_table' else 'disabled')
        note={'maxwellian':'P(E) ∝ √E exp(−E/θ); ⟨E⟩ = 3θ/2 в системе покоя. θ = Te [эВ] для электронов; θ = kBTi [эВ] для ионов. Температура задаёт начальную форму, затем распределение меняется в расчёте.',
            'energy_table':'Два столбца: E [эВ] и P(E) или EEPF. TXT, CSV, TSV, DAT, ASC, XLSX; произвольный заголовок. Форма нормируется автоматически. Углы изотропны, средняя скорость равна нулю.',
            'phase_table':'CSV: energy_eV;mu;probability_per_eV_mu. Прямоугольная сетка E, μ, μ = cos угла к z. Интеграл по E и μ должен быть 1; файл покрывает расчётную сетку. Это позволяет задать направленный пучок.'}[kind]
        self.hint.set(note)

    def browse(self):
        path=filedialog.askopenfilename(parent=self.winfo_toplevel(),filetypes=[('Распределения','*.txt *.csv *.tsv *.dat *.asc *.xlsx'),('Все файлы','*')])
        if path:
            try:path=os.path.relpath(path,self.root_path)
            except ValueError:pass
            self.file.set(path)

    def build(self):
        kind=KINDS[self.kind.get()]
        if kind=='maxwellian':
            temperature=finite(self.temperature.get(),'Температура частиц')
            if temperature<=0:raise ValueError('Температура распределения должна быть >0.')
            return dict(kind=kind,**{'Te_eV' if self.electron else 'Ti_K':temperature},match_moments=True)
        path=self.file.get().strip()
        if not path:raise ValueError('Выберите файл распределения или максвелловскую форму.')
        if not (self.root_path/path).is_file():raise ValueError('Файл распределения не найден: '+path)
        if kind=='phase_table':return dict(kind=kind,distribution_file=path)
        if not self.allow_energy:raise ValueError('Изотропный P(E) имеет u_z=0. Для направленного источника задайте максвелловское распределение или фазовый CSV.')
        from energy_input import read_energy_table
        form=FORMS[self.form.get()];read_energy_table(self.root_path/path,form)
        return dict(kind=kind,energy_distribution_file=path,energy_form=form)


class InitialDistributions(ttk.Frame):
    def __init__(self,parent,root,config,Te_var,Tgas_var):
        super().__init__(parent,padding=10);self.root_path=Path(root);self.config=config;self.overrides=deepcopy(config.get('gas',{}).get('distributions',{}));self.current='default';self.Tgas_var=Tgas_var
        left=ttk.Frame(self);left.grid(row=0,column=0,sticky='nsew',padx=(0,10));right=ttk.Frame(self);right.grid(row=0,column=1,sticky='nsew');self.columnconfigure((0,1),weight=1,uniform='distribution_columns');self.rowconfigure(0,weight=1)
        self.electron=DistributionEditor(left,'Электроны · EEDF',root,config.get('electron_initial',{}),Te_var,True);self.electron.pack(fill='x')
        ttk.Label(left,text='Начальные электроны: концентрация задаётся на вкладке «Расчёт». Для ионов энергетическое распределение называется IEDF. Оба ввода используют ту же нормировку.',wraplength=440).pack(anchor='w',pady=8)
        self.heavy=DistributionEditor(right,'Ионы и нейтральные частицы · IEDF',root,self.overrides.get('default',dict(Ti_K=Tgas_var.get())),electron=False);self.heavy.pack(fill='x')
        self.heavy.temperature=Tgas_var;self.heavy.temp_entry.configure(textvariable=Tgas_var)
        self.target=tk.StringVar(self,value='default');ttk.Label(right,text='Для кого задаётся распределение:').pack(anchor='w',pady=(8,3))
        self.target_box=ttk.Combobox(right,textvariable=self.target,values=['default'],state='normal');self.target_box.pack(fill='x');self.target_box.bind('<<ComboboxSelected>>',self.switch)
        self.use_default=tk.BooleanVar(self,value=False);ttk.Checkbutton(right,text='Общее распределение для этой частицы',variable=self.use_default).pack(anchor='w',pady=6)
        ttk.Label(right,text='default — все тяжёлые частицы. Можно задать отдельную форму для He+, H+, Li+ или точного ID уровня. Виды и уровни предлагаются из банка данных. Источники накачки имеют собственное распределение на вкладке «Источники в слое».',wraplength=450).pack(anchor='w',pady=5)
        self.override_temperature=tk.StringVar(self,value='300')
        self.grid_vars={};grid=ttk.LabelFrame(left,text='Сетки кинетической энергии',padding=8);grid.pack(fill='x',pady=10)
        for row,label,key,value in [(0,'Электроны · Emin, эВ','e_min',config.get('phase_grid',{}).get('energy_min_eV',1e-5)),(1,'Электроны · Emax, эВ / auto','e_max',config.get('phase_grid',{}).get('energy_max_eV','auto')),(2,'Электроны · число энергий (≥16)','e_points',config.get('phase_grid',{}).get('energy_points',32)),(3,'Электроны · число углов (≥4)','e_angles',config.get('phase_grid',{}).get('angular_points',4)),(4,'Тяжёлые частицы · Emin, эВ','h_min',config.get('heavy_phase_grid',{}).get('energy_min_eV',1e-6)),(5,'Тяжёлые частицы · Emax / auto','h_max',config.get('heavy_phase_grid',{}).get('energy_max_eV','auto')),(6,'Тяжёлые частицы · число энергий (≥16)','h_points',config.get('heavy_phase_grid',{}).get('energy_points',16)),(7,'Тяжёлые частицы · число углов (≥4)','h_angles',config.get('heavy_phase_grid',{}).get('angular_points',4))]:
            ttk.Label(grid,text=label).grid(row=row,column=0,sticky='w',pady=3);var=tk.StringVar(self,value=str(value));self.grid_vars[key]=var;ttk.Entry(grid,textvariable=var,width=14).grid(row=row,column=1,sticky='ew',pady=3)
        grid.columnconfigure(1,weight=1)

    def set_targets(self,names):self.target_box.configure(values=['default']+list(dict.fromkeys(names)))

    def save_current(self):
        if self.current!='default' and self.use_default.get():self.overrides.pop(self.current,None)
        else:self.overrides[self.current]=self.heavy.build()

    def switch(self,event=None):
        wanted=self.target.get().strip()
        if not wanted:return
        try:self.save_current()
        except (ValueError,OSError) as exc:self.target.set(self.current);messagebox.showerror('Распределение',str(exc),parent=self.winfo_toplevel());return
        self.current=wanted;self.use_default.set(wanted!='default' and wanted not in self.overrides)
        var=self.Tgas_var if wanted=='default' else self.override_temperature;self.heavy.temperature=var;self.heavy.temp_entry.configure(textvariable=var)
        self.heavy.load(self.overrides.get(wanted,self.overrides.get('default',dict(Ti_K=self.Tgas_var.get()))),keep_temperature=wanted=='default')

    def build(self):
        if self.target.get()!=self.current:self.switch()
        self.save_current();return self.electron.build(),deepcopy(self.overrides)

    def grids(self):
        from grid_validation import validate_phase_grid
        def count(prefix, field, label):
            value=finite(self.grid_vars[prefix+'_'+field].get(),label)
            if not value.is_integer():raise ValueError(label+': требуется целое число, получено '+str(value))
            return int(value)
        out=[]
        for prefix,key in [('e','phase_grid'),('h','heavy_phase_grid')]:
            grid=deepcopy(self.config.get(key,{}));low=finite(self.grid_vars[prefix+'_min'].get(),key+'.energy_min_eV');text=self.grid_vars[prefix+'_max'].get().strip();high='auto' if text=='auto' else finite(text,key+'.energy_max_eV')
            grid.update(energy_min_eV=low,energy_max_eV=high,energy_points=count(prefix,'points',key+'.energy_points'),angular_points=count(prefix,'angles',key+'.angular_points'))
            validate_phase_grid(grid,key,allow_auto=True);out.append(grid)
        return tuple(out)
