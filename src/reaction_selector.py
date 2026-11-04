"""Automatic reaction/reagent chooser shared by result plots and top-five."""
import tkinter as tk
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText
from reaction_catalog import ReactionCatalog,ReactionRequest,KIND_LABELS

MODES={'Top-5 по времени':'top5','Скорости выбранных реакций':'rates','Интегралы выбранных реакций':'integrals'}
DIRECTIONS={'Скорости событий':'events','Заселение и расселение':'both','Заселение':'source','Расселение':'sink','Участие в реакциях':'participation'}


class ReactionSelector(ttk.Frame):
    def __init__(self,parent,data,context_provider,on_build):
        super().__init__(parent);self.catalog=ReactionCatalog(data);self.context_provider=context_provider;self.on_build=on_build;self.loading=True;self.pending=None;self.matching_indices=[];self.last_target='all';self.error=''
        scroll=tk.Canvas(self,highlightthickness=0);bar=ttk.Scrollbar(self,orient='vertical',command=scroll.yview);scroll.configure(yscrollcommand=bar.set);bar.pack(side='right',fill='y');scroll.pack(side='left',fill='both',expand=True)
        body=ttk.Frame(scroll,padding=10);window=scroll.create_window(0,0,window=body,anchor='nw');body.bind('<Configure>',lambda _:scroll.configure(scrollregion=scroll.bbox('all')));scroll.bind('<Configure>',lambda e:scroll.itemconfigure(window,width=e.width))
        self.vars={};self.combos={};self.choice_labels={c.label:c.key for c in self.catalog.choices};self.target_labels={c.label:c.key for c in self.catalog.targets}
        def combo(row,label,key,values,default):
            ttk.Label(body,text=label).grid(row=row,column=0,sticky='w',pady=3,padx=(0,8));var=tk.StringVar(self,value=default);self.vars[key]=var;widget=ttk.Combobox(body,textvariable=var,values=values,state='readonly',width=35);widget.grid(row=row,column=1,sticky='ew',pady=3);self.combos[key]=widget
        combo(0,'Что изучаем','target',list(self.target_labels),self.catalog.targets[0].label)
        combo(1,'Что строить','mode',list(MODES),next(iter(MODES)))
        combo(2,'Вклад в цель','direction',['Скорости событий'],'Скорости событий')
        combo(3,'Входной реагент А','partner_a',list(self.choice_labels),self.catalog.choices[0].label)
        combo(4,'Входной реагент Б','partner_b',list(self.choice_labels),self.catalog.choices[0].label)
        self.kind_labels={'Любой процесс':'all'}
        for kind in sorted({r['kind'] for r in data.reactions}):self.kind_labels[KIND_LABELS.get(kind,kind)]=kind
        combo(5,'Тип процесса','kind',list(self.kind_labels),'Любой процесс')
        ttk.Label(body,text='Поиск (необязательно)').grid(row=6,column=0,sticky='w',pady=3);self.vars['search']=tk.StringVar(self,value='');ttk.Entry(body,textvariable=self.vars['search']).grid(row=6,column=1,sticky='ew',pady=3)
        self.nonzero=tk.BooleanVar(self,value=False);ttk.Checkbutton(body,text='Только ненулевые скорости в выбранной области',variable=self.nonzero).grid(row=7,column=0,columnspan=2,sticky='w',pady=4)
        self.hint=tk.StringVar(self,value='');ttk.Label(body,textvariable=self.hint,wraplength=550).grid(row=8,column=0,columnspan=2,sticky='w',pady=4)
        frame=ttk.Frame(body);frame.grid(row=9,column=0,columnspan=2,sticky='nsew',pady=4)
        self.tree=ttk.Treeview(frame,columns=('formula','peak','status'),show='headings',height=7,selectmode='extended')
        for key,label,width in [('formula','Реакция',360),('peak','max, см⁻³ с⁻¹',120),('status','Скорость',95)]:self.tree.heading(key,text=label);self.tree.column(key,width=width,minwidth=70,stretch=key=='formula')
        y=ttk.Scrollbar(frame,orient='vertical',command=self.tree.yview);x=ttk.Scrollbar(frame,orient='horizontal',command=self.tree.xview);self.tree.configure(yscrollcommand=y.set,xscrollcommand=x.set);self.tree.grid(row=0,column=0,sticky='nsew');y.grid(row=0,column=1,sticky='ns');x.grid(row=1,column=0,sticky='ew');frame.columnconfigure(0,weight=1);frame.rowconfigure(0,weight=1)
        self.tree.tag_configure('zero',foreground='#666666');self.tree.bind('<<TreeviewSelect>>',self.show_details);self.tree.bind('<Double-1>',self.plot_selected)
        self.summary=tk.StringVar(self,value='');ttk.Label(body,textvariable=self.summary,wraplength=550).grid(row=10,column=0,columnspan=2,sticky='w',pady=4)
        actions=ttk.Frame(body);actions.grid(row=11,column=0,columnspan=2,sticky='ew',pady=5)
        self.build_button=ttk.Button(actions,text='Построить',command=on_build);self.build_button.pack(side='left')
        ttk.Button(actions,text='Автовыбор 5',command=self.suggest).pack(side='left',padx=5);ttk.Button(actions,text='Все найденные',command=self.select_all).pack(side='left')
        more=ttk.Frame(body);more.grid(row=12,column=0,columnspan=2,sticky='w',pady=3)
        ttk.Button(more,text='Справочник всех реакций',command=self.show_catalogue).pack(side='left')
        ttk.Button(more,text='Сбросить фильтры',command=self.reset_filters).pack(side='left',padx=6)
        self.detail=ScrolledText(body,height=6,wrap='word',font=('Segoe UI',11),state='disabled');self.detail.grid(row=13,column=0,columnspan=2,sticky='nsew',pady=5)
        self.code_button=ttk.Button(actions,text='Код процесса',command=self.show_code);self.code_button.pack(side='left',padx=5)
        body.columnconfigure(1,weight=1);body.rowconfigure(9,weight=1)
        self.loading=False
        for var in [*self.vars.values(),self.nonzero]:var.trace_add('write',self.schedule_refresh)
        self.refresh()

    @property
    def metric(self):return MODES[self.vars['mode'].get()]

    @property
    def target(self):return self.target_labels[self.vars['target'].get()]

    @property
    def direction(self):return DIRECTIONS[self.vars['direction'].get()]

    @property
    def selected_ids(self):return [self.catalog.reactions[int(i)]['id'] for i in self.tree.selection() if self.tree.exists(i)]

    @property
    def matching_ids(self):return [self.catalog.reactions[k]['id'] for k in self.matching_indices]

    def context(self):
        context=self.context_provider()
        if getattr(self,'request_metric',None)=='top5' or (getattr(self,'request_metric',None) is None and self.metric=='top5'):context=dict(context,axis='time')
        return context

    def schedule_refresh(self,*_):
        if self.loading:return
        if self.pending is not None:self.after_cancel(self.pending)
        self.pending=self.after(120,self.refresh)

    def _values(self,key,choices):
        values=[c.label for c in choices];self.combos[key].configure(values=values)
        if self.vars[key].get() not in values:self.vars[key].set(values[0])

    def refresh(self):
        if self.pending is not None:self.after_cancel(self.pending);self.pending=None
        if self.loading:return
        old=set(self.selected_ids);self.loading=True;self.error=''
        try:
            target=self.target
            if target=='all':self.combos['direction'].configure(values=['Скорости событий'],state='disabled');self.vars['direction'].set('Скорости событий')
            else:
                allowed=['Заселение и расселение','Заселение','Расселение','Участие в реакциях'];self.combos['direction'].configure(values=allowed,state='readonly')
                if target!=self.last_target or self.vars['direction'].get() not in allowed:self.vars['direction'].set('Участие в реакциях' if target.startswith('component:') else 'Заселение и расселение')
            self.last_target=target;direction=self.direction
            self._values('partner_a',self.catalog.choices)
            a=self.choice_labels[self.vars['partner_a'].get()];compatible=self.catalog.choices_for_partners(a,target,direction);keys={c.key for c in compatible}
            self._values('partner_b',compatible+[c for c in self.catalog.choices if c.key not in keys]);b=self.choice_labels[self.vars['partner_b'].get()]
            candidates=self.catalog.matching(target=target,direction=direction,partner_a=a,partner_b=b)
            kinds={'all'}|{self.catalog.reactions[k]['kind'] for k in candidates};kind_values=[label for label,key in self.kind_labels.items() if key in kinds];self.combos['kind'].configure(values=kind_values)
            if self.vars['kind'].get() not in kind_values:self.vars['kind'].set('Любой процесс')
            context=self.context();indices=self.catalog.matching(target=target,direction=direction,partner_a=a,partner_b=b,kind=self.kind_labels[self.vars['kind'].get()],search=self.vars['search'].get(),nonzero=self.nonzero.get(),context=context)
            self.matching_indices=self.catalog.ranked(indices,target,direction,context);peaks=self.catalog.context_rates(**context).max(axis=0)
            self.tree.delete(*self.tree.get_children())
            for k in self.matching_indices:
                r=self.catalog.reactions[k];positive=peaks[k]>0;self.tree.insert('','end',iid=str(k),values=(r['reaction'],f'{peaks[k]:.5g}','ненулевая' if positive else '0'),tags=() if positive else ('zero',))
            preserved=[str(k) for k in self.matching_indices if self.catalog.reactions[k]['id'] in old]
            chosen=preserved or [str(k) for k in self.matching_indices[:5]]
            if chosen:self.tree.selection_set(chosen);self.tree.focus(chosen[0]);self.tree.see(chosen[0])
            count=sum(peaks[k]>0 for k in self.matching_indices);self.summary.set(f'Найдено: {len(indices)}; ненулевых: {count}; выбрано: {len(chosen)}. '+('Предложены наиболее быстрые реакции.' if not preserved and chosen else ''))
            if not indices:self._detail(self.catalog.absence_message())
            elif count==0:self._detail('Все найденные скорости равны нулю в этой области. Реакции остаются доступными; положительный top-5 отсутствует. Для просмотра нулей выберите «Скорости выбранных реакций».\n\n'+self.catalog.describe(self.matching_indices[0],context))
            else:self.show_details()
            self.build_button.configure(state='normal' if indices else 'disabled')
            self.code_button.configure(state='normal' if indices else 'disabled')
            if target.startswith('component:') and direction=='participation':note='Суммируются уровни компонента. «Участие» включает возбуждение внутри него; оно не меняет общее число частиц компонента.'
            elif self.metric=='top5':note='Top-5 выбирается в каждый момент среди всех найденных реакций. Таблица предлагает лидеров; фильтры задают область поиска.'
            else:note='Строятся выбранные строки. Скорость — число событий реакции, а не коэффициент k(T). Реагенты и источники данных показаны ниже.'
            self.hint.set(note)
        except (ValueError,KeyError,IndexError) as exc:
            self.error=str(exc);self.matching_indices=[];self.tree.delete(*self.tree.get_children());self.summary.set(self.error);self._detail(self.error);self.build_button.configure(state='disabled');self.code_button.configure(state='disabled')
        finally:self.loading=False

    def _detail(self,text):self.detail.configure(state='normal');self.detail.delete('1.0','end');self.detail.insert('1.0',text);self.detail.configure(state='disabled')

    def show_details(self,event=None):
        selection=self.tree.selection()
        if not selection:return
        item=self.tree.focus() if self.tree.focus() in selection else selection[0]
        try:
            from process_map import process_locations
            kind=self.catalog.reactions[int(item)]['kind'];location=process_locations(kind)[0]
            ordered=[int(item)]+[int(value) for value in selection if value!=item]
            self._detail(f'Описание выбранных реакций: {len(ordered)}.\n'+self.catalog.describe_many(ordered,self.context())+'\nРеализация активного процесса: '+location['file']+':'+str(location['line'])+' · '+location['function'])
        except (ValueError,KeyError):return
        if event is not None:self.summary.set(f'Найдено: {len(self.matching_indices)}; выбрано: {len(selection)}. Двойной щелчок: скорости.')

    def suggest(self):
        self.refresh();self.tree.selection_set([str(k) for k in self.matching_indices[:5]]);self.show_details()

    def select_all(self):self.tree.selection_set(self.tree.get_children());self.show_details()

    def reset_filters(self):
        self.loading=True
        for key in ('target','partner_a','partner_b'):self.vars[key].set(self.catalog.targets[0].label if key=='target' else self.catalog.choices[0].label)
        self.vars['kind'].set('Любой процесс');self.vars['search'].set('');self.nonzero.set(False);self.loading=False;self.refresh()

    def show_catalogue(self):
        from reaction_reference import ReactionReference
        return ReactionReference(self.winfo_toplevel(),self.catalog,self.context()).window

    def show_code(self):
        if not self.tree.selection():return None
        from process_map import show_process
        item=self.tree.focus() if self.tree.focus() in self.tree.selection() else self.tree.selection()[0]
        return show_process(self.winfo_toplevel(),self.catalog.reactions[int(item)]['kind'])

    def plot_selected(self,event=None):self.vars['mode'].set('Скорости выбранных реакций');self.refresh();self.on_build()

    def request(self,metric=None):
        self.request_metric=metric
        try:
            self.refresh()
            if self.error:raise ValueError(self.error)
            if not self.matching_indices:raise ValueError(self.catalog.absence_message())
            if not self.selected_ids:self.suggest()
            return ReactionRequest(metric or self.metric,self.target,self.direction,self.matching_ids,self.selected_ids)
        finally:self.request_metric=None

    def destroy(self):
        if self.pending is not None:
            try:self.after_cancel(self.pending)
            except tk.TclError:pass
            self.pending=None
        super().destroy()
