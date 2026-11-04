"""User-selected plots, readable labels, source starts and independent scales."""
from pathlib import Path
import queue
import threading
import tkinter as tk
from tkinter import ttk,filedialog,messagebox
from experiment_config import finite


class ResultBrowser:
    def __init__(self,parent,folder):
        from result_store import ResultData
        from result_plots import LABELS
        self.data=ResultData(folder);self.folder=Path(folder);self.pages=[];self.canvases=[];self.busy=False;self.closed=False;self.messages=queue.Queue();self.vars={}
        self.window=tk.Toplevel(parent);self.window.title('Результаты · выбор графиков');self.window.geometry('1320x860');self.window.minsize(1050,720)
        outer=ttk.Frame(self.window,padding=12);outer.pack(fill='both',expand=True)
        ttk.Label(outer,text='Результаты: выберите нужные данные',font=('Segoe UI',16,'bold')).pack(anchor='w')
        ttk.Label(outer,text='Графики строятся по отмеченным пунктам. Времена начала накачек показаны на временных графиках.').pack(anchor='w',pady=(5,10))
        quick=ttk.Frame(outer);quick.pack(fill='x',pady=(0,6));ttk.Label(quick,text='Быстрый выбор:').pack(side='left',padx=(0,8));self.quick_buttons={}
        for metric,label in [('components','Концентрации компонентов'),('levels','Населённости уровней'),('rates','Скорости реакций')]:
            button=ttk.Button(quick,text=label,command=lambda key=metric:self.select_plot(key));button.pack(side='left',padx=(0,6));self.quick_buttons[metric]=button
        scale_bar=ttk.Frame(outer);scale_bar.pack(fill='x',pady=(0,8));self.scale_widgets={}
        scale_values=('Линейный','Логарифмический','Симметричный лог')
        for key,label in [('x_scale','Масштаб X'),('y_scale','Масштаб Y')]:
            ttk.Label(scale_bar,text=label).pack(side='left',padx=(0,6));var=self.var(key,'Линейный');box=ttk.Combobox(scale_bar,textvariable=var,values=scale_values,state='readonly',width=23);box.pack(side='left',padx=(0,16));self.scale_widgets[key]=box
        body=ttk.Panedwindow(outer,orient='horizontal');self.body=body;body.pack(fill='both',expand=True);left=ttk.Frame(body,width=460);right=ttk.Frame(body);body.add(left,weight=0);body.add(right,weight=1)
        controls=ttk.Notebook(left);self.controls=controls;controls.pack(fill='both',expand=True,padx=(0,10));graphs=ttk.Frame(controls,padding=10);species=ttk.Frame(controls,padding=10);reactions=ttk.Frame(controls,padding=10);scales=ttk.Frame(controls,padding=10)
        for frame,label in [(graphs,'Графики'),(species,'Компоненты'),(reactions,'Реакции'),(scales,'Масштабы')]:controls.add(frame,text=label)
        self.components_tab=species;self.graphs_tab=graphs
        self.metric_vars={};self.metric_widgets={}
        for key,label in LABELS.items():
            var=tk.BooleanVar(self.window,value=key=='components');self.metric_vars[key]=var;widget=ttk.Checkbutton(graphs,text=label,variable=var);widget.pack(anchor='w',pady=3);self.metric_widgets[key]=widget
            if key in ('absorption','transmission','probe_distance') and 'probe_wavelength_nm' not in self.data.arrays:widget.configure(state='disabled')
        self.var('axis','По времени');ttk.Label(species,text='Ось динамики').pack(anchor='w');ttk.Combobox(species,textvariable=self.vars['axis'],values=('По времени','По расстоянию'),state='readonly').pack(fill='x',pady=4)
        elements=sorted({e for s in self.data.states for e in s.get('nuclei',{s['element']:1})})
        self.var('element','Все элементы');ttk.Combobox(species,textvariable=self.vars['element'],values=['Все элементы']+elements,state='readonly').pack(fill='x',pady=4)
        self.nonzero=tk.BooleanVar(self.window,value=False);ttk.Checkbutton(species,text='Скрыть уровни с нулевой населённостью',variable=self.nonzero,command=self.fill_components).pack(anchor='w',pady=5)
        self.inventory_info=tk.StringVar(self.window);ttk.Label(species,textvariable=self.inventory_info,wraplength=420).pack(anchor='w')
        self.component_list=self.listbox(species,'Компоненты / ионы / молекулы',6);self.level_list=self.listbox(species,'Уровни и реагенты',9)
        actions=ttk.Frame(species);actions.pack(fill='x',pady=6)
        ttk.Button(actions,text='Все в списках',command=self.select_all_components).pack(side='left');ttk.Button(actions,text='Снять выбор',command=lambda:(self.component_list.selection_clear(0,'end'),self.level_list.selection_clear(0,'end'))).pack(side='left',padx=6)
        self.fractions=tk.BooleanVar(self.window,value=False);ttk.Checkbutton(species,text='Показывать доли от числа тяжёлых частиц',variable=self.fractions).pack(anchor='w',pady=8)
        self.vars['element'].trace_add('write',lambda *_:self.fill_components())
        from reaction_selector import ReactionSelector
        self.reactions_tab=reactions
        self.reaction_panel=ReactionSelector(reactions,self.data,self.reaction_context,self.build_reactions);self.reaction_panel.pack(fill='both',expand=True)
        self.reaction_list=self.reaction_panel.tree;self.target_map=self.reaction_panel.target_labels
        self.vars['target']=self.reaction_panel.vars['target'];self.vars['direction']=self.reaction_panel.vars['direction'];self.vars['reaction_search']=self.reaction_panel.vars['search']
        for row,label,key,default,values in [(0,'Ось X','x_scale','Линейный',('Линейный','Логарифмический','Симметричный лог')),(1,'Ось Y','y_scale','Линейный',('Линейный','Логарифмический','Симметричный лог'))]:
            ttk.Label(scales,text=label).grid(row=row,column=0,sticky='w',pady=5);ttk.Combobox(scales,textvariable=self.vars[key],values=values,state='readonly',width=24).grid(row=row,column=1,sticky='ew',pady=5)
        for row,label,key,value in [(2,'Размер текста, pt','fontsize',16),(3,'Кривых на изображении','curves',5),(4,'Начало времени, с','time_min',''),(5,'Конец времени, с','time_max',''),(6,'Начало спектра, нм','wave_min',''),(7,'Конец спектра, нм','wave_max',''),(8,'λ для профиля по x, нм','probe_wave',1083.3)]:
            ttk.Label(scales,text=label).grid(row=row,column=0,sticky='w',pady=5);ttk.Entry(scales,textvariable=self.var(key,value),width=24).grid(row=row,column=1,sticky='ew',pady=5)
        self.mark=tk.BooleanVar(self.window,value=True);ttk.Checkbutton(scales,text='Отмечать начало каждой накачки',variable=self.mark).grid(row=9,column=0,columnspan=2,sticky='w',pady=8)
        self.export_csv=tk.BooleanVar(self.window,value=False);ttk.Checkbutton(scales,text='С PNG сохранять данные графика в CSV',variable=self.export_csv).grid(row=10,column=0,columnspan=2,sticky='w',pady=8)
        ttk.Label(scales,text='Пустые границы — весь диапазон.\nLog показывает только положительные значения.\nДля знаковых скоростей можно выбрать симметричный log.\nКривые делятся на отдельные изображения без уменьшения текста.',wraplength=420).grid(row=11,column=0,columnspan=2,sticky='w',pady=10)
        choices=ttk.Frame(right);choices.pack(fill='x',pady=(0,8))
        self.cell_values=['Среднее по слою']+[f'Ячейка {i+1} · x={x:.5g} см' for i,x in enumerate(self.data['position_cm'])];self.var('cell',self.cell_values[0]);ttk.Label(choices,text='Ячейка:').pack(side='left');ttk.Combobox(choices,textvariable=self.vars['cell'],values=self.cell_values,state='readonly',width=27).pack(side='left',padx=6)
        self.time_values=[f'{i}: {t:.7g} с' for i,t in enumerate(self.data['time_s'])];self.var('time',self.time_values[-1]);ttk.Label(choices,text='Момент:').pack(side='left',padx=(8,0));ttk.Combobox(choices,textvariable=self.vars['time'],values=self.time_values,state='readonly',width=22).pack(side='left',padx=6)
        self.preview=ttk.Notebook(right);self.preview.pack(fill='both',expand=True);self.placeholder=ttk.Frame(self.preview,padding=25);self.preview.add(self.placeholder,text='Предпросмотр')
        ttk.Label(self.placeholder,text='1. Отметьте нужные графики.\n2. Выберите компоненты, уровни или реакции.\n3. Нажмите «Построить выбранное».\n4. Сохраните нужные изображения в PNG.',font=('Segoe UI',13),justify='left').pack(anchor='w',pady=20)
        source_lines=[]
        for source in self.data.meta.get('sources',[]):
            line=f"{source['id']}: старт {source['start_s']:.5g} с"
            if source['mode']=='radiation':line+=f", {source['normal_power_W_cm2']:.5g} Вт/см²; фотоны {source['minimum_photon_eV']:.6g}…{source['maximum_photon_eV']:.6g} эВ; средняя {source['mean_photon_eV']:.6g} эВ"
            else:line+=f", n={source.get('density_cm3','старый формат')} см⁻³, u={source.get('velocity_cm_s','старый формат')} см/с"
            source_lines.append(line)
        if source_lines:ttk.Label(self.placeholder,text='\n'.join(source_lines),wraplength=650,justify='left').pack(anchor='w',pady=15)
        footer=ttk.Frame(outer);footer.pack(fill='x',pady=(10,0));self.build_button=ttk.Button(footer,text='Построить выбранное',command=self.build);self.build_button.pack(side='left');self.save_button=ttk.Button(footer,text='Сохранить PNG',command=self.save_dialog,state='disabled');self.save_button.pack(side='left',padx=8)
        ttk.Button(footer,text='Таблицы и файлы',command=self.open_folder).pack(side='left');ttk.Button(footer,text='Закрыть',command=self.close).pack(side='right')
        self.status=tk.StringVar(value='Данные загружены. Ничего лишнего автоматически не строится.');ttk.Label(outer,textvariable=self.status,wraplength=1250).pack(anchor='w',pady=(8,0))
        self.fill_components();self.fill_reactions()
        self.controls.bind('<<NotebookTabChanged>>',self.on_controls_tab)
        for key in ('cell','time','axis','time_min','time_max'):self.vars[key].trace_add('write',self.reaction_panel.schedule_refresh)
        self.window.protocol('WM_DELETE_WINDOW',self.close);self.poll_id=self.window.after(100,self.poll)

    def var(self,key,value):self.vars[key]=tk.StringVar(self.window,value=str(value));return self.vars[key]

    def listbox(self,parent,label,height):
        ttk.Label(parent,text=label).pack(anchor='w',pady=(7,3));frame=ttk.Frame(parent);frame.pack(fill='both',expand=True);box=tk.Listbox(frame,selectmode='extended',exportselection=False,height=height,font=('Segoe UI',11));scroll=ttk.Scrollbar(frame,orient='vertical',command=box.yview);box.configure(yscrollcommand=scroll.set);box.pack(side='left',fill='both',expand=True);scroll.pack(side='right',fill='y');return box

    def fill_components(self):
        import numpy as np
        element=self.vars['element'].get();self.component_keys=[];self.level_keys=[];n=self.data['n_cm3'];self.component_list.delete(0,'end');self.level_list.delete(0,'end')
        eligible=[]
        for j,s in enumerate(self.data.states):
            nuclei=s.get('nuclei',{s['element']:1})
            if element!='Все элементы' and element not in nuclei:continue
            if self.nonzero.get() and not np.any(n[:,:,j]>0):continue
            eligible.append(j);self.level_keys.append(s['id']);self.level_list.insert('end',s.get('display_label',s['label']))
        for name,ids in self.data.meta['components'].items():
            if not set(ids).intersection(eligible):continue
            self.component_keys.append(name);self.component_list.insert('end',name)
        self.component_keys.append('e⁻');self.component_list.insert('end','e⁻ — электроны');self.component_list.selection_set(0,'end')
        best=sorted(eligible,key=lambda j:-float(n[:,:,j].max()))[:5]
        for j in best:self.level_list.selection_set(self.level_keys.index(self.data.states[j]['id']))
        occupied=sum(bool(np.any(n[:,:,j]>0)) for j in eligible)
        self.inventory_info.set(f'Показано уровней: {len(eligible)} из {len(self.data.states)}; ненулевых: {occupied}. Нулевая населённость не исключает уровень из выбора.')

    def select_all_components(self):self.component_list.selection_set(0,'end');self.level_list.selection_set(0,'end')

    @property
    def reaction_keys(self):return self.reaction_panel.matching_ids

    def fill_reactions(self):self.reaction_panel.refresh()

    def show_reaction(self,event=None):self.reaction_panel.show_details(event)

    def reaction_context(self):
        def number(key):
            var=self.vars.get(key);value=var.get().strip() if var else ''
            return finite(value,key) if value else None
        cell=self.cell_values.index(self.vars['cell'].get())-1 if 'cell' in self.vars else -1
        time_index=self.time_values.index(self.vars['time'].get()) if 'time' in self.vars else -1
        return dict(cell=cell,axis='space' if self.vars['axis'].get()=='По расстоянию' else 'time',time_index=time_index,time_min=number('time_min'),time_max=number('time_max'))

    def on_controls_tab(self,event=None):
        if not hasattr(self,'build_button'):return
        reaction_active=self.controls.select()==str(self.reactions_tab)
        if self.body.winfo_width()>800:self.body.sashpos(0,600 if reaction_active else 460)
        if reaction_active:self.reaction_panel.refresh()

    def select_plot(self,metric):
        """An explicit plot choice; browsing tabs never changes that choice."""
        for key,var in self.metric_vars.items():var.set(key==metric)
        if metric in ('components','levels','electrons'):self.controls.select(self.components_tab)
        elif metric in ('rates','integrals','top5'):
            from reaction_selector import MODES
            self.reaction_panel.vars['mode'].set(next(label for label,key in MODES.items() if key==metric));self.controls.select(self.reactions_tab)
        else:self.controls.select(self.graphs_tab)
        self.status.set('Выбрано: '+self.metric_widgets[metric].cget('text')+'. Уточните строки и нажмите «Построить выбранное».')

    def build_reactions(self):
        self.select_plot(self.reaction_panel.metric);self.build()

    def options(self):
        from result_plots import PlotOptions
        scales={'Линейный':'linear','Логарифмический':'log','Симметричный лог':'symlog'}
        def optional(key):return finite(self.vars[key].get(),key) if self.vars[key].get().strip() else None
        metrics=[key for key,var in self.metric_vars.items() if var.get()]
        if not metrics:raise ValueError('Отметьте хотя бы один график.')
        reaction_metrics=[key for key in metrics if key in ('rates','integrals','top5')]
        request=self.reaction_panel.request('top5' if 'top5' in reaction_metrics else reaction_metrics[0]) if reaction_metrics else None
        return PlotOptions(metrics=metrics,components=[self.component_keys[i] for i in self.component_list.curselection()],levels=[self.level_keys[i] for i in self.level_list.curselection()],reactions=request.selected_ids if request else None,reaction_filter=request.matching_ids if request else None,
            axis='time' if self.vars['axis'].get()=='По времени' else 'space',cell=self.cell_values.index(self.vars['cell'].get())-1,time_index=self.time_values.index(self.vars['time'].get()),target=request.target if request else 'all',direction=request.direction if request else 'source',
            x_scale=scales[self.vars['x_scale'].get()],y_scale=scales[self.vars['y_scale'].get()],fontsize=int(self.vars['fontsize'].get()),curves_per_figure=int(self.vars['curves'].get()),mark_starts=self.mark.get(),fractions=self.fractions.get(),time_min=optional('time_min'),time_max=optional('time_max'),wavelength_min=optional('wave_min'),wavelength_max=optional('wave_max'),probe_wavelength_nm=finite(self.vars['probe_wave'].get(),'λ профиля'))

    def build_sync(self):
        from result_plots import build_plots
        self.show_pages(build_plots(self.data,self.options()));return self.pages

    def build(self):
        if self.busy:return
        try:options=self.options()
        except ValueError as exc:messagebox.showerror('Графики',str(exc),parent=self.window);return
        self.busy=True;self.build_button.configure(state='disabled');self.save_button.configure(state='disabled');self.status.set('Построение выбранных графиков…')
        def job():
            try:
                from result_plots import build_plots
                self.messages.put(('pages',build_plots(self.data,options)))
            except Exception as exc:self.messages.put(('error',str(exc)))
        threading.Thread(target=job,daemon=True).start()

    def show_pages(self,pages):
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
        for canvas in self.canvases:canvas.get_tk_widget().destroy()
        self.canvases=[]
        for tab in self.preview.tabs():self.preview.nametowidget(tab).destroy()
        self.pages=pages
        for i,page in enumerate(pages):
            tab=ttk.Frame(self.preview);self.preview.add(tab,text=f'{i+1}. {page.title[:24]}')
            preview_dpi=min(100.,max(60.,(self.preview.winfo_width()-38)/page.size_inches[0]))
            page.figure.set_size_inches(page.size_inches,forward=False);page.figure.set_dpi(preview_dpi)
            view=ttk.Frame(tab);view.pack(fill='both',expand=True);viewport=tk.Canvas(view,highlightthickness=0,bg='white');hs=ttk.Scrollbar(view,orient='horizontal',command=viewport.xview);vs=ttk.Scrollbar(view,orient='vertical',command=viewport.yview)
            viewport.configure(xscrollcommand=hs.set,yscrollcommand=vs.set);viewport.grid(row=0,column=0,sticky='nsew');vs.grid(row=0,column=1,sticky='ns');hs.grid(row=1,column=0,sticky='ew');view.rowconfigure(0,weight=1);view.columnconfigure(0,weight=1)
            width,height=page.figure.get_size_inches()*preview_dpi;fixed=ttk.Frame(viewport,width=int(width),height=int(height));fixed.pack_propagate(False);viewport.create_window(0,0,window=fixed,anchor='nw',width=int(width),height=int(height));viewport.configure(scrollregion=(0,0,int(width),int(height)))
            canvas=FigureCanvasTkAgg(page.figure,master=fixed);canvas.draw();canvas.get_tk_widget().pack();toolbar=NavigationToolbar2Tk(canvas,tab,pack_toolbar=False);toolbar.update();toolbar.pack(side='bottom',fill='x');self.canvases.append(canvas)
        self.save_button.configure(state='normal' if pages else 'disabled');self.status.set(f'Построено изображений: {len(pages)}. PNG сохраняются с крупным текстом, независимо от размера окна.')

    def poll(self):
        if self.closed:return
        try:
            while True:
                kind,value=self.messages.get_nowait();self.busy=False;self.build_button.configure(state='normal')
                if kind=='pages':self.show_pages(value)
                else:self.status.set(value);messagebox.showerror('Графики',value,parent=self.window)
        except queue.Empty:pass
        self.poll_id=self.window.after(100,self.poll)

    def save_to(self,folder):
        from result_plots import save_pages
        files=save_pages(self.pages,folder,self.export_csv.get());self.status.set(f'Сохранено PNG: {len(files)} · {folder}');return files

    def save_dialog(self):
        if not self.pages:return
        path=filedialog.askdirectory(parent=self.window,initialdir=self.folder,title='Куда сохранить выбранные PNG')
        if path:
            try:self.save_to(Path(path)/'рисунки')
            except OSError as exc:messagebox.showerror('Сохранение',str(exc),parent=self.window)

    def open_folder(self):
        from launcher import open_path
        open_path(self.folder)

    def close(self):
        self.closed=True
        if getattr(self,'poll_id',None) is not None:
            try:self.window.after_cancel(self.poll_id)
            except tk.TclError:pass
            self.poll_id=None
        self.window.destroy()


def open_results(parent,folder):return ResultBrowser(parent,folder)
