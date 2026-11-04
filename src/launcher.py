"""Окно запуска многокомпонентного решателя и импорта физических данных.

Windows: двойной щелчок по START_WINDOWS.bat. Только стандартная библиотека
нужна самому окну; numpy/scipy/matplotlib устанавливаются при запуске расчёта.
"""
from project_version import VERSION
import os
import json
from pathlib import Path
import queue
import subprocess
import sys
import threading
from progress import Progress,parse_event,describe,duration
from grid_validation import validate_configuration_grids

ROOT = Path(__file__).resolve().parent.parent
CODE=Path(__file__).resolve().parent
PRESETS = tuple('configs/'+name for name in ('gas_probe.json','study_01_decay.json','gas_probe_user_fixed.json','pulsed_sources.json','gas_probe_H_He_Li.json','mixture_H_He.json','mixture_H_He_chemistry.json','mixture_He_chemistry.json','mixture_H_He_Li.json','mixture_H_He_Na.json','mixture_He_ion_decay.json','slab_H_He_Li_sources.json'))


def configuration(path,validate=False):
    path=Path(path)
    if path.suffix.lower()=='.json':
        config=json.loads(path.read_text(encoding='utf-8-sig'))
        engines={'AtomicData-1.0':'prepare_species.py','MixtureKineticSlab-1.8':'mixture_kinetics.py','MixtureKineticSlab-2.1':'mixture_kinetics.py'}
        if config.get('schema') not in engines:raise ValueError('Неизвестная схема JSON')
        if validate:validate_configuration_grids(config)
        return config,engines[config['schema']]
    raise ValueError('Текущий решатель использует JSON: MixtureKineticSlab-2.1.')


def settings(path):
    config,engine=configuration(path)
    if engine=='prepare_species.py':
        output=config.get('atomic_data',{}).get('prepared_directory','atomic_data')
        return ROOT/'data',(Path(path).resolve().parent/output).resolve()
    return ROOT/'data',(Path(path).resolve().parent/config.get('output_directory','output_point')).resolve()


def python_executable():
    executable = Path(sys.executable)
    # У pythonw нет стандартной консоли. Дочерний python пишет в PIPE,
    # а CREATE_NO_WINDOW запрещает создание отдельного консольного окна.
    if os.name == "nt" and executable.name.lower() == "pythonw.exe":
        console = executable.with_name("python.exe")
        if console.exists():
            executable = console
    return str(executable)


def command(arguments, emit):
    environment = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
    options = {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {}
    with subprocess.Popen([python_executable(), "-u", *map(str, arguments)], cwd=ROOT,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            encoding="utf-8", errors="replace", env=environment, **options) as process:
        for line in process.stdout:
            emit(line)
        return process.wait()


def run_model(config_path, emit, on_progress=None):
    """Подготовить окружение, выполнить обычный CLI и вернуть его код/папку.

    Параметры и скорости реакций не изменяются. Код 2 означает нарушение
    диагностик и передаётся окну отдельно от ошибки запуска.
    """
    config_path = Path(config_path).resolve()
    with (ROOT / "launcher.log").open("w", encoding="utf-8") as log:
        def report(line):
            log.write(line)
            log.flush()
            event=parse_event(line)
            if event is not None:
                if on_progress is not None:on_progress(event)
                emit(describe(event)+'\n')
            else:emit(line)

        progress=Progress('json',sink=lambda line:report(line+'\n'))

        try:
            data, output = settings(config_path)
            config, engine = configuration(config_path,validate=True)
            progress.start('Подготовка окружения','Проверка библиотек')
            report("Проверка библиотек…\n")
            if command(["-c", "import numpy, scipy, matplotlib, threadpoolctl, numba, mpmath"], lambda line: None):
                progress.module('Установка библиотек Python',force=True)
                report("Установка библиотек. Для этого нужен интернет…\n")
                # Готовые wheels: не запускаем компиляторы в пользовательской системе.
                if command(["-m", "pip", "install", "--only-binary=:all:", "-r", ROOT / "requirements.txt"], report):
                    raise RuntimeError("Не удалось установить библиотеки. Подробности в журнале окна и launcher.log.")
            report(f"\nРасчёт: {config_path.name}\nОжидайте; окно остаётся доступным.\n\n")
            progress.module('Загрузка численных модулей Python',force=True)
            arguments=[CODE / engine, '--config', config_path]
            if engine in ('kinetic_slab.py','mixture_kinetics.py','prepare_species.py'):arguments+=['--progress','json']
            code = command(arguments, report)
            if code not in (0, 2):
                raise RuntimeError("Расчёт остановлен с ошибкой. Подробности выше и в launcher.log.")
            if engine=='mixture_kinetics.py' and not (output/'данные_для_графиков.npz').exists() and config.get('diagnostics',{}).get('plots_on_finish',True):
                progress.start('Графики','Населённости, EEDF и скорости ионов')
                if command([CODE/'plot_mixture.py',output],report):raise RuntimeError('Расчёт записан, но построение графиков не удалось')
                progress.done('Графики сохранены')
            return code, output
        except Exception as error:
            report("\nОшибка: " + str(error) + "\n")
            raise


def open_path(path, edit=False):
    path = Path(path)
    if os.name == "nt":
        if edit:
            subprocess.Popen(["notepad.exe", str(path)], creationflags=subprocess.CREATE_NO_WINDOW)
        else:
            os.startfile(str(path))
    else:
        subprocess.Popen(["open" if sys.platform == "darwin" else "xdg-open", str(path)])


def main():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    from tkinter.scrolledtext import ScrolledText

    window = tk.Tk()
    window.title(f"Plasma kinetics {VERSION} — состав, источники и выбранные графики")
    window.geometry("1020x650")
    window.minsize(680, 430)
    frame = ttk.Frame(window, padding=14)
    frame.pack(fill="both", expand=True)
    ttk.Label(frame, text="Кинетика и атомные данные", font=("Segoe UI", 16, "bold")).pack(anchor="w")
    ttk.Label(frame, text="Выберите параметры и нажмите «Запустить». Результаты сохраняются в CSV и PNG.").pack(anchor="w", pady=(4, 12))
    selection = ttk.Frame(frame)
    selection.pack(fill="x")
    ttk.Label(selection, text="Параметры:").pack(side="left", padx=(0, 8))
    chosen = tk.StringVar(value=PRESETS[0])
    dropdown = ttk.Combobox(selection, textvariable=chosen, values=PRESETS, state="readonly")
    dropdown.pack(side="left", fill="x", expand=True)
    busy = False
    messages = queue.Queue()

    def config_path():
        return (ROOT / chosen.get()).resolve()

    def browse():
        path = filedialog.askopenfilename(parent=window, initialdir=ROOT, filetypes=[("Параметры JSON/TOML", "*.json *.toml"),("Все файлы", "*")])
        if path:
            chosen.set(path)

    browse_button = ttk.Button(selection, text="Выбрать файл…", command=browse)
    browse_button.pack(side="left", padx=(8, 0))
    ttk.Label(frame, text="Смеси: атомные уровни, ионы, общая EEDF и межкомпонентные реакции. «Элементы и данные» открывает импорт и проверку источников.",
              wraplength=820).pack(anchor="w", pady=(6, 12))
    controls = ttk.Frame(frame)
    controls.pack(fill="x")
    status = tk.StringVar(value="Готово к запуску")
    log = ScrolledText(frame, wrap="word", height=19, font=("Consolas", 10), state="disabled")

    def worker(path):
        try:
            code, output = run_model(path, lambda line: messages.put(("text", line)),lambda event:messages.put(('progress',event)))
            messages.put(("done", (code, output)))
        except Exception as error:
            messages.put(("error", str(error)))

    def run():
        nonlocal busy
        path = config_path()
        if not path.exists():
            messagebox.showerror("Нет файла", "Не найден файл параметров: " + str(path), parent=window)
            return
        busy = True
        for button in (run_button, parameters_button, browse_button):
            button.configure(state="disabled")
        dropdown.configure(state="disabled")
        log.configure(state="normal")
        log.delete("1.0", "end")
        log.configure(state="disabled")
        status.set("Выполняется подготовка и расчёт…")
        progress.configure(mode='indeterminate',value=0);progress.start(12)
        eta.set('Оставшееся время: оценивается после начала измеряемого этапа')
        threading.Thread(target=worker, args=(path,), daemon=True).start()

    def edit():
        try:
            _,engine=configuration(config_path())
            if engine!='mixture_kinetics.py':
                open_path(config_path(),edit=True);return
            from experiment_editor import open_editor
            open_editor(window,config_path())
        except (OSError,ValueError,KeyError) as error:
            messagebox.showerror("Параметры", str(error), parent=window)

    def results():
        try:
            _, output = settings(config_path())
            if (output/'данные_для_графиков.npz').exists():
                from result_browser import open_results
                open_results(window,output)
            elif output.exists():
                messagebox.showinfo('Прежние результаты','Для окна выбранных графиков выполните расчёт в версии 2.6. Прежние таблицы можно открыть в папке результатов.',parent=window)
                open_path(output)
            else:
                messagebox.showinfo("Результаты", "Папка появится после расчёта. Нажмите «Запустить».", parent=window)
        except (OSError, ValueError, KeyError,ImportError) as error:
            messagebox.showerror("Результаты", str(error), parent=window)

    def data_panel():
        """A separate, usable import block; the solver never guesses sources."""
        if busy:
            messagebox.showinfo('Выполняется расчёт','Дождитесь завершения перед изменением физических данных.',parent=window)
            return
        from datetime import datetime
        panel=tk.Toplevel(window);panel.title('Элементы и данные из сети');panel.geometry('920x650')
        body=ttk.Frame(panel,padding=14);body.pack(fill='both',expand=True)
        ttk.Label(body,text='Элементы, сечения и коэффициенты',font=('Segoe UI',15,'bold')).pack(anchor='w')
        ttk.Label(body,text='Атомные уровни/A/σ(E) и тепловые k(T) хранятся раздельно. Импорт UMIST создаёт кандидатов; неизвестные уровни и законы реакций перечисляются в отчёте.',wraplength=860).pack(anchor='w',pady=(5,12))
        cfg,engine=configuration(config_path());base=config_path().parent
        element_dir=(base/cfg.get('atomic_data',{}).get('element_directory','data/elements')).resolve()
        table=ttk.Treeview(body,columns=('stages','provider'),show='tree headings',height=6)
        table.heading('#0',text='Элемент');table.heading('stages',text='Импортируемые стадии');table.heading('provider',text='Электронные сечения')
        table.column('#0',width=90);table.column('stages',width=200);table.column('provider',width=510);table.pack(fill='x')
        def refresh():
            table.delete(*table.get_children())
            for file in sorted(element_dir.glob('*/atom.json')):
                document=json.loads(file.read_text(encoding='utf-8-sig'))
                entries=document['species'];providers=', '.join(str(v.get('electron_cross_sections') or ('native CSV' if v.get('native_csv') else 'не подключены')) for v in entries.values())
                table.insert('','end',text=document['element'],values=(', '.join(entries),providers))
        refresh()
        names=tk.StringVar(value=', '.join(cfg.get('atomic_data',{}).get('species',['He'])))
        selection=ttk.Frame(body);selection.pack(fill='x',pady=10)
        ttk.Label(selection,text='Виды в конфигурации:').pack(side='left');ttk.Entry(selection,textvariable=names,state='readonly').pack(side='left',fill='x',expand=True,padx=8)
        def save_species():
            try:
                from experiment_editor import open_editor
                open_editor(panel,config_path())
            except (ValueError,OSError) as exc:messagebox.showerror('Виды',str(exc),parent=panel)
        ttk.Button(selection,text='Изменить состав',command=save_species).pack(side='left')
        actions=ttk.Frame(body);actions.pack(fill='x',pady=(0,10))
        data_log=ScrolledText(body,wrap='word',height=15,font=('Consolas',10));data_log.pack(fill='both',expand=True)
        data_messages=queue.Queue();data_busy=False;buttons=[]
        def execute(arguments):
            nonlocal data_busy,busy
            if data_busy:return
            data_busy=True;busy=True
            for b in buttons:b.configure(state='disabled')
            for b in (run_button,parameters_button,browse_button):b.configure(state='disabled')
            dropdown.configure(state='disabled')
            data_log.delete('1.0','end')
            def job():
                try:
                    if command(['-c','import numpy, scipy, matplotlib, threadpoolctl, numba, mpmath'],lambda line:None):
                        data_messages.put(('text','Установка библиотек для импорта данных…\n'))
                        if command(['-m','pip','install','--only-binary=:all:','-r',ROOT/'requirements.txt'],lambda line:data_messages.put(('text',line))):raise RuntimeError('Не удалось установить зависимости')
                    code=command([CODE/'data_import.py',*arguments],lambda line:data_messages.put(('text',line)))
                    data_messages.put(('done',code))
                except Exception as exc:data_messages.put(('text',str(exc)+'\n'));data_messages.put(('done',1))
            threading.Thread(target=job,daemon=True).start()
        def import_network():
            from atomic_core import Species
            try:
                elements=sorted({Species.parse(v.strip()).element for v in names.get().split(',') if v.strip()})
                destination=ROOT/'cache'/'network_candidates'/datetime.now().strftime('%Y%m%d_%H%M%S')
                execute(['network','--elements',*elements,'--out',destination])
            except ValueError as exc:messagebox.showerror('UMIST',str(exc),parent=panel)
        def register():
            path=filedialog.askopenfilename(parent=panel,title='ElementData-1: дескриптор нового элемента',filetypes=[('JSON','*.json')])
            if path:execute(['register','--descriptor',path,'--element-directory',element_dir])
        def import_cross_sections():
            source=filedialog.askopenfilename(parent=panel,title='Экспорт LXCat/BOLSIG',filetypes=[('Text','*.txt'),('Все','*')])
            if not source:return
            mapping=filedialog.askopenfilename(parent=panel,title='Явное соответствие уровней и единиц',filetypes=[('JSON','*.json')])
            if not mapping:return
            destination=filedialog.asksaveasfilename(parent=panel,title='Куда записать electron_cross_sections.csv',defaultextension='.csv')
            if destination:execute(['lxcat','--file',source,'--mapping',mapping,'--out',destination])
        for label,callback in [('Подготовить атомные данные',lambda:execute(['prepare','--config',config_path()])),('Загрузить сеть UMIST',import_network),('Добавить элемент',register),('Импорт σ(E)',import_cross_sections),('Проверить смесь',lambda:execute(['validate','--config',config_path()]))]:
            button=ttk.Button(actions,text=label,command=callback);button.pack(side='left',padx=(0,5));buttons.append(button)
        ttk.Label(body,text='Порядок добавления элемента и форматы файлов: docs/ADDING_ELEMENTS_RU.md. Новое имя без физических данных не становится готовой моделью.',wraplength=860).pack(anchor='w',pady=(10,0))
        def poll_data():
            nonlocal data_busy,busy
            if not panel.winfo_exists():return
            while True:
                try:kind,value=data_messages.get_nowait()
                except queue.Empty:break
                if kind=='text':data_log.insert('end',value);data_log.see('end')
                else:
                    data_busy=False;busy=False
                    for b in buttons:b.configure(state='normal')
                    for b in (run_button,parameters_button,browse_button):b.configure(state='normal')
                    dropdown.configure(state='readonly')
                    data_log.insert('end','\nЗавершено.\n' if value==0 else '\nИмпорт/проверка остановлены; причина выше.\n');refresh()
            panel.after(150,poll_data)
        panel.after(150,poll_data)
        def close_panel():
            if data_busy:messagebox.showinfo('Импорт выполняется','Дождитесь завершения операции с данными.',parent=panel)
            else:panel.destroy()
        panel.protocol('WM_DELETE_WINDOW',close_panel)

    run_button = ttk.Button(controls, text="Запустить", command=run)
    run_button.pack(side="left")
    parameters_button = ttk.Button(controls, text="Состав и спектр", command=edit)
    parameters_button.pack(side="left", padx=8)
    ttk.Button(controls, text="Результаты", command=results).pack(side="left")
    ttk.Button(controls,text='Элементы и данные',command=data_panel).pack(side='left',padx=8)
    def slab_panel():
        if busy:
            messagebox.showinfo('Выполняется расчёт','Дождитесь завершения перед изменением параметров.',parent=window);return
        try:
            from spatial_editor import open_editor
            open_editor(window,config_path())
        except (OSError,ValueError,KeyError) as exc:messagebox.showerror('Slab',str(exc),parent=window)
    ttk.Button(controls,text='Дополнительно: Slab',command=slab_panel).pack(side='left',padx=4)
    ttk.Button(controls,text='JSON',command=lambda:open_path(config_path(),edit=True)).pack(side='left',padx=4)
    ttk.Label(frame, textvariable=status, wraplength=820).pack(anchor="w", pady=(12, 5))
    progress = ttk.Progressbar(frame, mode="indeterminate")
    progress.pack(fill="x", pady=(0, 10))
    eta=tk.StringVar(value='Процент и оставшееся время относятся к текущему этапу.')
    ttk.Label(frame,textvariable=eta,wraplength=820).pack(anchor='w',pady=(0,6))
    log.pack(fill="both", expand=True)
    ttk.Label(frame, text="Журнал также записывается в launcher.log. При первом запуске загрузки могут занять несколько минут.",
              wraplength=820).pack(anchor="w", pady=(8, 0))

    def poll():
        nonlocal busy
        for _ in range(100):
            try:
                kind, value = messages.get_nowait()
            except queue.Empty:
                break
            if kind == "text":
                log.configure(state="normal")
                log.insert("end", value)
                log.see("end")
                log.configure(state="disabled")
            elif kind=='progress':
                fraction=value.get('fraction')
                if fraction is None:
                    if str(progress['mode'])!='indeterminate':
                        progress.stop();progress.configure(mode='indeterminate');progress.start(12)
                else:
                    progress.stop();progress.configure(mode='determinate',maximum=100,value=100*fraction)
                status.set(value['stage']+' · '+value.get('module',''))
                remaining=value.get('eta_s')
                detail=value.get('detail','')
                label='Осталось до конца этапа ≈ '+duration(remaining) if remaining is not None else 'Оставшееся время оценивается'
                if value.get('status')=='done':label='Этап завершён'
                eta.set(label+(' · '+detail if detail else ''))
            else:
                busy = False
                progress.stop()
                for button in (run_button, parameters_button, browse_button):
                    button.configure(state="normal")
                dropdown.configure(state="readonly")
                if kind == "done":
                    code, output = value
                    status.set("Готово. Откройте «Результаты»." if code == 0 else
                               "Расчёт завершён с замечаниями к данным. См. О_результатах.txt; графики доступны в «Результатах».")
                    progress.configure(mode='determinate',value=100)
                    eta.set('Выберите нужные графики в «Результатах».' if code==0 else 'Код 2: расчёт сохранён; физические ограничения перечислены в О_результатах.txt.')
                else:
                    status.set("Ошибка. Подробности в журнале ниже.")
                    log.configure(state="normal")
                    log.insert("end", "\n" + value + "\n")
                    log.see("end")
                    log.configure(state="disabled")
                    messagebox.showerror("Не удалось завершить расчёт", value, parent=window)
        window.after(100, poll)

    def close():
        if busy:
            messagebox.showinfo("Расчёт выполняется", "Дождитесь завершения расчёта перед закрытием окна.", parent=window)
        else:
            window.destroy()

    window.protocol("WM_DELETE_WINDOW", close)
    window.after(100, poll)
    window.mainloop()
