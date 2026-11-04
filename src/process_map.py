"""Read-only, AST-derived locations for the actual implemented operators."""
from functools import lru_cache
from pathlib import Path
import ast

ROOT=Path(__file__).resolve().parent.parent


@lru_cache(maxsize=64)
def definition(module,qualified_name):
    path=ROOT/'src'/(module+'.py');source=path.read_text(encoding='utf-8');node=ast.parse(source)
    for part in qualified_name.split('.'):
        node=next(child for child in node.body if isinstance(child,(ast.ClassDef,ast.FunctionDef,ast.AsyncFunctionDef)) and child.name==part)
    return dict(file='src/'+path.name,path=str(path),function=qualified_name,line=node.lineno,end_line=node.end_lineno,
        code='\n'.join(f'{i:5d}  {line}' for i,line in enumerate(source.splitlines(),1) if node.lineno<=i<=node.end_lineno))


def process_locations(kind):
    common=[('Скорость события','mixture_kinetics','MixtureSlab.rates'),('Стехиометрия','reaction_data','equation'),('Сборка сети','mixture_kinetics','MixtureSlab._network')]
    if kind in ('excitation','deexcitation','ionization'):
        entries=[('Регистрация канала','mixture_kinetics','MixtureSlab._atomic_channels'),('Сумма по EEDF','kinetic_kernels','impact_rates'),('Изменение EEDF','kinetic_kernels','impact_events')]
    elif kind in ('spontaneous','two_photon','photoexcitation','stimulated_emission'):
        entries=[('Атомные фотоканалы','mixture_kinetics','MixtureSlab._atomic_channels'),('Перенос света','slab_transport','solve_transfer')]
    elif kind in ('photoionization','radiative_recombination','stimulated_recombination'):
        entries=[('PI и рекомбинация','mixture_kinetics','MixtureSlab._atomic_channels'),('Фото-сечения','mixture_photo','PhotoBank.__init__'),('Парная рекомбинация','kinetic_kernels','build_rr')]
    else:
        entries=[('Чтение сети','mixture_chemistry','Chemistry.__init__'),('Формулы k(T)','mixture_chemistry','rate_fit'),('Парное ядро','mixture_chemistry','Chemistry.kernel'),('Продукты реакции','chemistry_kernels','products')]
    return [dict(definition(module,name),label=label) for label,module,name in entries+common]


def show_process(parent,kind):
    import tkinter as tk
    from tkinter import ttk
    from tkinter.scrolledtext import ScrolledText
    window=tk.Toplevel(parent);window.title('Физический процесс · исходный код');window.geometry('1100x760');window.minsize(850,550)
    ttk.Label(window,text='Это текущий код программы. Числа слева — строки файла.',padding=10).pack(anchor='w')
    book=ttk.Notebook(window);book.pack(fill='both',expand=True,padx=10,pady=(0,10))
    for location in process_locations(kind):
        tab=ttk.Frame(book,padding=8);book.add(tab,text=location['label']);ttk.Label(tab,text=location['file']+':'+str(location['line'])+' · '+location['function'],font=('Segoe UI',12,'bold')).pack(anchor='w',pady=(0,7))
        text=ScrolledText(tab,wrap='none',font=('Consolas',11));text.pack(fill='both',expand=True);horizontal=ttk.Scrollbar(tab,orient='horizontal',command=text.xview);text.configure(xscrollcommand=horizontal.set);horizontal.pack(fill='x');text.insert('1.0',location['code']);text.configure(state='disabled')
    return window
