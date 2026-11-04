"""Standard-library Tk editor for spatial sources and the two slab walls."""
from copy import deepcopy
import json
import math
import os
from pathlib import Path


def number(text):
    value=float(str(text).strip().replace(',','.'))
    if not math.isfinite(value):raise ValueError('Нужно конечное число')
    return value


def validate(config):
    length=number(config['thickness_cm']);cells=config['cells'];end=number(config['solver']['t_end_s'])
    if length<=0 or not isinstance(cells,int) or cells<1 or end<=0:raise ValueError('L>0, cells>=1, t_end>0')
    identifiers=set()
    for item in config.get('spatial_sources',[]):
        if not item.get('enabled',True):continue
        if not item['id'] or item['id'] in identifiers:raise ValueError('У каждого источника должен быть уникальный id')
        identifiers.add(item['id']);position=number(item['position_cm']);width=number(item.get('width_cm',0.))
        if width<0 or not 0<=position-width/2<=position+width/2<=length:raise ValueError('Источник должен лежать внутри слоя')
        if not 0<=number(item.get('start_s',0.))<number(item.get('stop_s',end))<=end:raise ValueError('0<=start_s<stop_s<=t_end_s')
        if item['mode']=='radiation':
            if width!=0 or number(item.get('scale',1.))<0 or not 0<number(item.get('mu',1.))<=1:raise ValueError('Излучение: width=0, scale>=0, 0<mu<=1')
            if item.get('spectrum_mode') not in ('probe','flat') and not item.get('spectrum_file'):raise ValueError('Выберите форму спектра источника')
        elif item['mode'] in ('ions','electrons','particles'):
            from light_input import particle_inputs
            row=particle_inputs(item)
            if number(row.get('rate_cm2_s',0.))<0:raise ValueError('Концентрация пучка должна быть >=0')
        else:raise ValueError('Неизвестный режим источника')
    for surface in config.get('particle_boundary',{}).get('walls',{}).values():
        for key,value in surface.items():
            if key=='overrides':continue
            if not 0<=number(value)<=1:raise ValueError('Вероятности на стенке должны лежать в [0,1]')
    return config


def open_editor(parent,path):
    import tkinter as tk
    from tkinter import ttk,filedialog,messagebox
    path=Path(path);cfg=json.loads(path.read_text(encoding='utf-8-sig'));items=deepcopy(cfg.get('spatial_sources',[]))
    panel=tk.Toplevel(parent);panel.title('Slab — источники, стенки и пространственный вывод');panel.geometry('1040x790')
    panel.minsize(900,700);body=ttk.Frame(panel,padding=12);body.pack(fill='both',expand=True)
    tabs=ttk.Notebook(body);tabs.pack(fill='both',expand=True)
    geo=ttk.Frame(tabs,padding=12);src=ttk.Frame(tabs,padding=12)
    tabs.add(geo,text='Слой и стенки');tabs.add(src,text='Пространственные источники')
    fields={}
    def entry(frame,label,key,value,row,column=0,width=18):
        ttk.Label(frame,text=label).grid(row=row,column=column,sticky='w',padx=4,pady=4)
        var=tk.StringVar(value=str(value));ttk.Entry(frame,textvariable=var,width=width).grid(row=row,column=column+1,sticky='we',padx=4,pady=4);fields[key]=var
        return var
    entry(geo,'Длина слоя, см','L',cfg.get('thickness_cm',1.),0)
    entry(geo,'Число ячеек','cells',cfg.get('cells',8),1)
    entry(geo,'Конец расчёта, с','end',cfg.get('solver',{}).get('t_end_s',1e-6),2)
    ttk.Label(geo,text='Границы частиц').grid(row=3,column=0,sticky='w',padx=4,pady=4)
    boundary=tk.StringVar(value=cfg.get('particle_boundary',{}).get('kind','periodic'))
    ttk.Combobox(geo,textvariable=boundary,values=('walls','vacuum','reservoir','periodic'),state='readonly').grid(row=3,column=1,sticky='we')
    wallvars={};walls=cfg.get('particle_boundary',{}).get('walls',{})
    for col,side,label in [(0,'left','Начало слоя · x=0'),(2,'right','Конец слоя · x=L')]:
        box=ttk.LabelFrame(geo,text=label,padding=10);box.grid(row=4,column=col,columnspan=2,sticky='nsew',padx=4,pady=15)
        defaults={'quenching_probability':1.,'neutralization_probability':1.,'sticking_probability':0.,'electron_absorption_probability':1.}
        labels=['Тушение возбуждённых','Нейтрализация ионов','Прилипание тяжёлых','Поглощение электронов']
        for row,((key,value),text) in enumerate(zip(defaults.items(),labels)):
            var=tk.StringVar(value=str(walls.get(side,{}).get(key,value)));wallvars[(side,key)]=var
            ttk.Label(box,text=text+' · P').grid(row=row,column=0,sticky='w',pady=4)
            ttk.Entry(box,textvariable=var,width=10).grid(row=row,column=1,padx=8)
    ttk.Label(geo,text='Вероятности зависят от материала и задаются пользователем. Возвращённые частицы отражаются зеркально; тушение переводит нейтральный атом в основное состояние. Нейтрализация ионов, ток на стенку и выделенная энергия входят в балансы.',wraplength=900).grid(row=5,column=0,columnspan=4,sticky='w',pady=8)
    chem=cfg.get('chemistry',{});estimated=tk.BooleanVar(value=chem.get('allow_estimated_rates',False))
    ttk.Checkbutton(geo,text='Разрешить явно обозначенные оценочные сечения He*–Li (модель захвата, см. отчёт)',variable=estimated).grid(row=6,column=0,columnspan=4,sticky='w',pady=8)
    dataset=tk.StringVar(value=chem.get('datasets',{}).get('He_H_chemiionization','RD82'))
    ttk.Label(geo,text='He*–H: набор коэффициентов').grid(row=7,column=0,sticky='w')
    ttk.Combobox(geo,textvariable=dataset,values=('GM25','RD82'),state='readonly').grid(row=7,column=1,sticky='we')
    ttk.Label(geo,text='GM25: Garcia Muñoz 2025, 500–10 000 K; RD82: прежний вариант. Вне диапазона коэффициент не экстраполируется.',wraplength=900).grid(row=8,column=0,columnspan=4,sticky='w',pady=8)
    from source_editor import SourceEditor
    source_editor=SourceEditor(src,cfg,path.parent);source_editor.pack(fill='both',expand=True)
    source_editor.geometry_provider=lambda:(number(fields['L'].get()),number(fields['end'].get()))
    ttk.Label(body,text='Результаты: четыре понятные таблицы; выбранные графики — через кнопку «Результаты».',wraplength=990).pack(anchor='w',pady=8)
    actions=ttk.Frame(body);actions.pack(fill='x')
    def save():
        try:
            # Read current content to preserve changes from the other editor.
            updated=json.loads(path.read_text(encoding='utf-8-sig'));updated.update(cells=int(fields['cells'].get()),thickness_cm=number(fields['L'].get()),spatial_sources=source_editor.build(number(fields['L'].get()),number(fields['end'].get())))
            updated.setdefault('solver',{})['t_end_s']=number(fields['end'].get())
            updated['solver']['max_step_s']=min(updated['solver'].get('max_step_s',updated['solver']['t_end_s']/20),updated['solver']['t_end_s'])
            bound=updated.setdefault('particle_boundary',{});bound['kind']=boundary.get()
            if boundary.get()=='walls':
                bound['walls']={side:dict(walls.get(side,{}),**{key:number(var.get()) for (name,key),var in wallvars.items() if name==side}) for side in ('left','right')}
            else:bound.pop('walls',None)
            chem=updated.setdefault('chemistry',{});chem['allow_estimated_rates']=estimated.get();chem.setdefault('datasets',{})['He_H_chemiionization']=dataset.get()
            updated.setdefault('diagnostics',{}).update(output_mode='compact',plots_on_finish=False)
            validate(updated);path.write_text(json.dumps(updated,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
            messagebox.showinfo('Сохранено','Параметры слоя и источников сохранены.',parent=panel)
        except (ValueError,OSError,KeyError) as exc:messagebox.showerror('Параметры Slab',str(exc),parent=panel)
    ttk.Button(actions,text='Сохранить параметры',command=save).pack(side='left');ttk.Button(actions,text='Закрыть',command=panel.destroy).pack(side='right')
    return panel
