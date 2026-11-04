"""Explicit beam density/speed and normalised light source controls."""
from copy import deepcopy
from pathlib import Path
import os
import tkinter as tk
from tkinter import ttk,filedialog,messagebox
from experiment_config import finite
from light_input import particle_inputs,light_spectrum

SOURCE_MODES={'Ионы':'ions','Электроны':'electrons','Свет':'radiation','Нейтральные частицы':'particles'}


class SourceEditor(ttk.Frame):
    def __init__(self,parent,config,root,probe_provider=None):
        super().__init__(parent,padding=10);self.root_path=Path(root);self.config=config;self.probe_provider=probe_provider or (lambda:config.get('probe',{}));self.rows=deepcopy(config.get('spatial_sources',[]));self.selected=None;self.vars={};self.fields={};self.dirty=False
        ttk.Label(self,text='Независимые источники в заданной точке слоя',font=('Segoe UI',13,'bold')).pack(anchor='w')
        self.table=ttk.Treeview(self,columns=('mode','x','start','stop','enabled'),show='tree headings',height=5,selectmode='browse');self.table.heading('#0',text='Имя');self.table.column('#0',width=190)
        for key,label,width in [('mode','Тип',120),('x','x, см',100),('start','Начало, с',125),('stop','Конец, с',125),('enabled','Включён',80)]:self.table.heading(key,text=label);self.table.column(key,width=width)
        self.table.pack(fill='x',pady=8);self.table.bind('<<TreeviewSelect>>',self.load)
        form=ttk.Frame(self);form.pack(fill='x');geo=ttk.LabelFrame(form,text='Положение и время',padding=8);geo.grid(row=0,column=0,sticky='nsew',padx=(0,10));params=ttk.LabelFrame(form,text='Параметры источника',padding=8);params.grid(row=0,column=1,sticky='nsew');form.columnconfigure(0,weight=1);form.columnconfigure(1,weight=1)
        end=config.get('solver',{}).get('t_end_s',1e-7);length=config.get('thickness_cm',.01)
        self.vars['mode']=tk.StringVar(self,value='Ионы');ttk.Label(geo,text='Тип').grid(row=0,column=0,sticky='w');ttk.Combobox(geo,textvariable=self.vars['mode'],values=list(SOURCE_MODES),state='readonly',width=25).grid(row=0,column=1,sticky='ew',pady=4)
        for row,label,key,value in [(1,'Имя','id','Источник 1'),(2,'Положение x, см','x',length/2),(3,'Ширина, см (0 = точка)','width',0.),(4,'Начало импульса, с','start',0.),(5,'Конец импульса, с','stop',end)]:self.field(geo,row,label,key,value)
        self.enabled=tk.BooleanVar(self,value=True);ttk.Checkbutton(geo,text='Включить этот источник',variable=self.enabled).grid(row=6,column=0,columnspan=2,sticky='w',pady=5)
        self.particle=ttk.Frame(params);self.light=ttk.Frame(params)
        for row,label,key,value in [(0,'Частица / уровень','species','H+'),(1,'Концентрация пучка, см⁻³','density',1e7),(2,'Скорость u_z, см/с','velocity',1e6)]:self.field(self.particle,row,label,key,value)
        from distribution_editor import DistributionEditor
        self.vars['temperature']=tk.StringVar(self,value='300.0');self.distribution_editor=DistributionEditor(self.particle,'Распределение частиц источника',root,dict(Ti_K=300.),self.vars['temperature'],electron=False,allow_energy=False);self.distribution_editor.grid(row=3,column=0,columnspan=2,sticky='ew')
        self.fields['temperature']=(self.distribution_editor.temp_label,self.distribution_editor.temp_entry)
        ttk.Label(self.particle,text='Знак u_z задаёт направление; инжекция через плоскость = n × |u_z|. При u_z = 0 направленная накачка равна нулю. Для фазового файла средняя скорость должна совпадать с u_z.',wraplength=430).grid(row=4,column=0,columnspan=2,sticky='w',pady=6)
        self.use_probe=tk.BooleanVar(self,value=True);ttk.Checkbutton(self.light,text='Форма из вкладки «Пробный спектр»',variable=self.use_probe).grid(row=0,column=0,columnspan=2,sticky='w',pady=5)
        self.field(self.light,1,'Или файл формы','file','');ttk.Button(self.light,text='Выбрать файл…',command=self.browse).grid(row=2,column=1,sticky='w')
        self.vars['unit']=tk.StringVar(self,value='auto');ttk.Label(self.light,text='Единицы X файла').grid(row=3,column=0,sticky='w');ttk.Combobox(self.light,textvariable=self.vars['unit'],values=('auto','nm','A','um','eV','Hz'),state='readonly',width=22).grid(row=3,column=1,pady=4)
        for row,label,key,value in [(4,'Множитель, Вт/см²','power',1e-4),(5,'Косинус угла |μ|','mu',1.)]:self.field(self.light,row,label,key,value)
        self.vars['direction']=tk.StringVar(self,value='right');ttk.Label(self.light,text='Направление света').grid(row=6,column=0,sticky='w');ttk.Combobox(self.light,textvariable=self.vars['direction'],values=('right','left','both'),state='readonly',width=22).grid(row=6,column=1,pady=4)
        ttk.Button(self.light,text='Показать энергетику спектра',command=self.energy_dialog).grid(row=7,column=0,columnspan=2,sticky='w',pady=8)
        self.info=tk.StringVar(self,value='');ttk.Label(self,textvariable=self.info,wraplength=900).pack(anchor='w',pady=8)
        actions=ttk.Frame(self);actions.pack(fill='x',pady=4)
        for label,action in [('Добавить / применить',self.apply_dialog),('Новый',self.new),('Удалить выбранный',self.remove)]:ttk.Button(actions,text=label,command=action).pack(side='left',padx=(0,8))
        self.vars['mode'].trace_add('write',self.change_mode);self.refresh();self.change_mode();self.loading=False
        self.geometry_provider=lambda:(float(self.config.get('thickness_cm',.01)),float(self.config.get('solver',{}).get('t_end_s',1e-7)))
        for var in [*self.vars.values(),self.enabled,self.use_probe,self.distribution_editor.kind,self.distribution_editor.file,self.distribution_editor.form]:var.trace_add('write',self.mark_dirty)

    def mark_dirty(self,*_):
        if not self.loading:self.dirty=True

    def field(self,parent,row,label,key,value):
        label_widget=ttk.Label(parent,text=label,wraplength=225);label_widget.grid(row=row,column=0,sticky='w',pady=4,padx=(0,8));var=tk.StringVar(self,value=str(value));self.vars[key]=var
        if key=='species':
            from particle_choices import particle_choices
            entry=ttk.Combobox(parent,textvariable=var,values=particle_choices(self.root_path,self.config),width=22)
        else:entry=ttk.Entry(parent,textvariable=var,width=22)
        entry.grid(row=row,column=1,sticky='ew',pady=4);self.fields[key]=(label_widget,entry)

    def change_mode(self,*_):
        mode=SOURCE_MODES[self.vars['mode'].get()]
        self.particle.grid_remove();self.light.grid_remove();(self.light if mode=='radiation' else self.particle).grid(row=0,column=0,sticky='nsew')
        self.distribution_editor.set_particle(mode=='electrons')
        for widget in self.fields['species']:
            if mode=='electrons':widget.grid_remove()
            else:widget.grid()
        if mode=='electrons' and self.vars['temperature'].get()=='300.0':self.vars['temperature'].set('1.0')
        self.info.set('Форма нормируется по интегралу. Множитель задаёт суммарную мощность через плоскость; при both она делится между направлениями.' if mode=='radiation' else 'Источник добавляет частицы в указанной области. Электроны и ионы включаются независимо; автоматической нейтрализации нет.')

    def refresh(self):
        self.table.delete(*self.table.get_children())
        for i,row in enumerate(self.rows):self.table.insert('','end',iid=str(i),text=row['id'],values=(next(k for k,v in SOURCE_MODES.items() if v==row['mode']),row.get('position_cm',0.),row.get('start_s',0.),row.get('stop_s',self.config.get('solver',{}).get('t_end_s',1e-7)),'да' if row.get('enabled',True) else 'нет'))

    def load(self,event=None):
        if not self.table.selection():return
        self.loading=True
        self.selected=int(self.table.selection()[0]);row=self.rows[self.selected];dist=row.get('distribution',{});mode=row['mode']
        self.vars['mode'].set(next(k for k,v in SOURCE_MODES.items() if v==mode));self.enabled.set(row.get('enabled',True));self.use_probe.set(row.get('spectrum_mode')=='probe')
        speed=row.get('velocity_cm_s',dist.get('uz_cm_s',0.));density=row.get('density_cm3',float(row.get('rate_cm2_s',0.))/abs(speed) if speed else '')
        mapping=dict(id=row['id'],x=row.get('position_cm',0.),width=row.get('width_cm',0.),start=row.get('start_s',0.),stop=row.get('stop_s',self.config.get('solver',{}).get('t_end_s',1e-7)),species=row.get('species','H+'),density=density,velocity=speed,temperature=dist.get('Te_eV',1.) if mode=='electrons' else dist.get('Ti_K',300.),file=row.get('spectrum_file',''),unit=row.get('spectrum_options',{}).get('x_unit','auto'),power=row.get('multiplier_W_cm2',1e-4),mu=row.get('mu',1.),direction=row.get('direction','right'))
        for key,value in mapping.items():self.vars[key].set(str(value))
        self.distribution_editor.load(dist,keep_temperature=True)
        if mode!='radiation' and 'density_cm3' not in row:self.info.set('Старый источник задан потоком. Для изменения задайте концентрацию и скорость. Неизменённый источник сохраняется в прежнем виде.')
        self.loading=False;self.dirty=False

    def new(self):
        self.selected=None;self.table.selection_remove(*self.table.selection());self.vars['id'].set('Источник '+str(len(self.rows)+1))
        self.dirty=False

    def remove(self):
        if self.table.selection():self.rows.pop(int(self.table.selection()[0]));self.selected=None;self.refresh();self.dirty=False

    def browse(self):
        file=filedialog.askopenfilename(parent=self.winfo_toplevel(),filetypes=[('Спектры','*.txt *.csv *.tsv *.dat *.asc *.xlsx'),('Все файлы','*')])
        if file:
            try:file=os.path.relpath(file,self.root_path)
            except ValueError:pass
            self.vars['file'].set(file);self.use_probe.set(False)

    def light_config(self):
        return dict(normalization='integral',spectrum_mode='probe' if self.use_probe.get() else 'file',spectrum_file=self.vars['file'].get(),spectrum_options={'x_unit':self.vars['unit'].get(),'y_kind':'relative'},multiplier_W_cm2=finite(self.vars['power'].get(),'Множитель света'),mu=finite(self.vars['mu'].get(),'μ'))

    def energy_dialog(self):
        try:
            row=self.light_config();data=light_spectrum(row,self.root_path,self.probe_provider());info=data.metadata
            self.info.set(f"Фотоны: {info['minimum_photon_eV']:.7g}…{info['maximum_photon_eV']:.7g} эВ; средняя {info['mean_photon_eV']:.7g} эВ. Мощность {info['normal_power_W_cm2']:.7g} Вт/см²; фотонный поток {info['photon_flux_cm2_s']:.7g} см⁻² с⁻¹.")
        except (ValueError,OSError) as exc:messagebox.showerror('Энергетика спектра',str(exc),parent=self.winfo_toplevel())

    def build_row(self):
        mode=SOURCE_MODES[self.vars['mode'].get()];row=dict(id=self.vars['id'].get().strip(),mode=mode,enabled=self.enabled.get(),position_cm=finite(self.vars['x'].get(),'Положение'),width_cm=finite(self.vars['width'].get(),'Ширина'),start_s=finite(self.vars['start'].get(),'Начало'),stop_s=finite(self.vars['stop'].get(),'Конец'))
        if not row['id']:raise ValueError('Укажите имя источника.')
        if mode=='radiation':
            row.update(self.light_config(),direction=self.vars['direction'].get())
            if row['enabled']:light_spectrum(row,self.root_path,self.probe_provider())
        else:
            row.update(density_cm3=finite(self.vars['density'].get(),'Концентрация пучка'),velocity_cm_s=finite(self.vars['velocity'].get(),'Скорость пучка'))
            row['distribution']=self.distribution_editor.build()
            if mode!='electrons':
                from particle_choices import input_name
                row['species']=input_name(self.vars['species'].get())
            particle_inputs(row)
        return row

    def apply(self):
        row=self.build_row();candidate=deepcopy(self.rows)
        if self.selected is None:candidate.append(row)
        else:candidate[self.selected]=row
        length,end=self.geometry_provider();self.validate(candidate,length,end)
        self.rows=candidate;self.selected=None;self.refresh();self.dirty=False;return row

    def apply_dialog(self):
        try:self.apply()
        except (ValueError,OSError) as exc:messagebox.showerror('Источник',str(exc),parent=self.winfo_toplevel())

    @staticmethod
    def validate(rows,length,end):
        names=set()
        for row in rows:
            if not row.get('enabled',True):continue
            if row['id'] in names:raise ValueError('Повтор имени источника: '+row['id'])
            names.add(row['id']);x=float(row.get('position_cm',0.));width=float(row.get('width_cm',0.));start=float(row.get('start_s',0.));stop=float(row.get('stop_s',end))
            if width<0 or not 0<=x-width/2<=x+width/2<=length:raise ValueError('Источник должен целиком лежать в слое.')
            if not 0<=start<stop<=end:raise ValueError('Источник: 0<=начало<конец<=время расчёта.')
            if row['mode']=='radiation' and width!=0:raise ValueError('Источник света — плоскость, ширина должна быть 0.')

    def build(self,length,end):
        if self.dirty:self.apply()
        self.validate(self.rows,length,end);return deepcopy(self.rows)
