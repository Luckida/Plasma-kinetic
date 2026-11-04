"""Searchable reference for every saved channel and data-bank alternative."""
import tkinter as tk
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText
from reaction_catalog import folded,KIND_LABELS
from reaction_data import kinetics_summary


class ReactionReference:
    def __init__(self,parent,catalog,context):
        self.catalog=catalog;self.context=context;self.rows=catalog.data.meta.get('reaction_bank',catalog.reactions);self.pending=None
        self.window=tk.Toplevel(parent);w=self.window;w.title('Справочник · все доступные реакции');w.geometry('1200x820');w.minsize(900,600)
        outer=ttk.Frame(w,padding=12);outer.pack(fill='both',expand=True)
        ttk.Label(outer,text='Все реакции для состава этого расчёта',font=('Segoe UI',15,'bold')).pack(anchor='w')
        ttk.Label(outer,text='Нулевые скорости не скрываются. Альтернативы и выключенные каналы описаны отдельно; их скорости не подменяются нулями.',wraplength=1100).pack(anchor='w',pady=6)
        self.search=tk.StringVar(w);ttk.Entry(outer,textvariable=self.search).pack(fill='x',pady=6)
        self.count=tk.StringVar(w);ttk.Label(outer,textvariable=self.count).pack(anchor='w')
        frame=ttk.Frame(outer);frame.pack(fill='both',expand=True)
        self.tree=ttk.Treeview(frame,columns=('reaction','kind','status'),show='headings',height=11,selectmode='browse')
        for key,label,width in [('reaction','Реакция',630),('kind','Процесс',200),('status','Наличие в сети',290)]:self.tree.heading(key,text=label);self.tree.column(key,width=width,minwidth=100)
        self.tree.grid(row=0,column=0,sticky='nsew');ys=ttk.Scrollbar(frame,orient='vertical',command=self.tree.yview);ys.grid(row=0,column=1,sticky='ns');xs=ttk.Scrollbar(frame,orient='horizontal',command=self.tree.xview);xs.grid(row=1,column=0,sticky='ew');self.tree.configure(yscrollcommand=ys.set,xscrollcommand=xs.set);frame.rowconfigure(0,weight=1);frame.columnconfigure(0,weight=1)
        self.detail=ScrolledText(outer,height=14,wrap='word',font=('Segoe UI',12),state='disabled');self.detail.pack(fill='both',expand=True,pady=8)
        self.tree.bind('<<TreeviewSelect>>',self.describe);self.search.trace_add('write',self.schedule);w.bind('<Destroy>',self.destroyed);self.refresh()

    def schedule(self,*_):
        if self.pending is not None:self.window.after_cancel(self.pending)
        self.pending=self.window.after(160,self.refresh)

    def refresh(self):
        self.pending=None;query=folded(self.search.get());self.tree.delete(*self.tree.get_children());count=0
        for i,row in enumerate(self.rows):
            haystack=folded(row['reaction']+' '+row['id']+' '+row.get('source','')+' '+KIND_LABELS.get(row['kind'],row['kind'])+' '+kinetics_summary(row))
            if query and query not in haystack:continue
            self.tree.insert('','end',iid=str(i),values=(row['reaction'],KIND_LABELS.get(row['kind'],row['kind']),row.get('status','В расчётной сети')));count+=1
        self.count.set(f'Показано: {count} из {len(self.rows)}. Каналов расчётной сети: {len(self.catalog.reactions)}.')
        children=self.tree.get_children()
        if children:self.tree.selection_set(children[0]);self.describe()
        else:self.set_detail('Для этого поиска записей нет. Очистите поле поиска, чтобы увидеть все реакции.')

    def set_detail(self,text):self.detail.configure(state='normal');self.detail.delete('1.0','end');self.detail.insert('1.0',text);self.detail.configure(state='disabled')

    def describe(self,event=None):
        selected=self.tree.selection()
        if not selected:return
        row=self.rows[int(selected[0])]
        if row['id'] in self.catalog.by_id:text=self.catalog.describe(self.catalog.by_id[row['id']],self.context)
        else:text=row['reaction']+'\nРеагенты: '+row.get('left','')+'\nПродукты: '+row.get('right','')+'\nID: '+row['id']+'\nСтатус: '+row.get('status','Вне расчётной сети')+'\nИсточник данных: '+str(row.get('source',''))+'\n'+kinetics_summary(row)+'\nСкорость в этом результате не рассчитывалась. Выбор записи в справочнике не включает процесс в решатель.'
        self.set_detail(text)

    def destroyed(self,event):
        if event.widget==self.window and self.pending is not None:
            try:self.window.after_cancel(self.pending)
            except tk.TclError:pass
            self.pending=None
