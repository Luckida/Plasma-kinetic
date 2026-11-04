"""Generate code line links and an exhaustive registered chemistry inventory."""
from pathlib import Path
import argparse
import ast
import csv
import json

ROOT=Path(__file__).resolve().parent.parent
OPERATORS='''mixture_kinetics|MixtureSlab.__init__|Состав и начальные распределения
mixture_kinetics|MixtureSlab.enabled|Флаги процессов
mixture_kinetics|MixtureSlab._atomic_channels|Атомные электронные и фотонные каналы
mixture_kinetics|MixtureSlab._spectrum|Общая фотонная сетка и внешняя накачка
mixture_kinetics|MixtureSlab._network|Стехиометрия общей сети
mixture_kinetics|MixtureSlab._pairs|Парная геометрия и рекомбинационные ядра
mixture_kinetics|MixtureSlab.rates|Частоты атомных процессов и химические события
mixture_kinetics|MixtureSlab._check_data|Недостающие атомные процессы и диапазоны
mixture_kinetics|MixtureSlab.field_E|Предписанное поле и экспериментальный Poisson
mixture_kinetics|MixtureSlab.rr_values|Парные RR/SR скорости
mixture_kinetics|MixtureSlab.field|Общее поле излучения
mixture_kinetics|MixtureSlab.optical_matrices|Поглощение и излучение активного базиса
mixture_kinetics|MixtureSlab._binary|e–neutral, e–ion и BGK столкновения
mixture_kinetics|MixtureSlab._local|Связанный шаг населённостей, EEDF и химии
mixture_kinetics|MixtureSlab._transport|Ускорение и пространственный перенос
mixture_kinetics|MixtureSlab.solve|Расщепление и контроль шага
mixture_kinetics|MixtureSlab.total_energy|Полная энергия частиц и уровней
mixture_kinetics|MixtureSlab.snapshot|Балансы, события и сохранённый срез
mixture_kinetics|conservative_profile|Согласование спектрального профиля
mixture_kinetics|isolated_line_profile|Конечное окно Voigt-профиля кинетических линий
kinetic_kernels|impact_rates|Сумма σv по EEDF
kinetic_kernels|impact_events|Электронные потери и продукты неупругих событий
kinetic_kernels|pair_geometry|Относительные электрон–тяжёлые скорости
kinetic_kernels|build_rr|Обратный PI оператор Милна
kinetic_kernels|rr_emission|Спонтанная рекомбинация
kinetic_kernels|rr_stimulated|Стимулированная рекомбинация
kinetic_kernels|rr_heavy_products|Тяжёлые продукты рекомбинации
kinetic_kernels|momentum_transfer_rate|Экранированное кулоновское транспортное сечение
kinetic_kernels|binary_generator_rates|Генератор парных упругих столкновений
mixture_photo|PhotoBank.__init__|Выбор фото-сечений
mixture_photo|hydrogen_sigma|Кулоновские фото-интегралы
mixture_photo|verner_sigma|Фото-фиты Verner
mixture_chemistry|append_components|Молекулярные носители и энергетические якоря
mixture_chemistry|rate_fit|Все используемые законы k(T)
mixture_chemistry|Chemistry.__init__|Чтение сети и выбор альтернатив
mixture_chemistry|Chemistry.kernel|Парное σv и моментное замыкание
mixture_chemistry|Chemistry.event_rates|Скорости химических событий
mixture_chemistry|Chemistry.evolve|Стоки и продукты химии
mixture_chemistry|Chemistry.photon_emission|Фотоны радиационной перезарядки
chemistry_kernels|geometry|Относительная энергия тяжёлых пар
chemistry_kernels|products|Кинематика продуктов в координатах Якоби
particle_kinetics|PhaseGrid.__init__|Энергетическая и угловая сетки
particle_kinetics|PhaseGrid.initial|Начальная EEDF
particle_kinetics|PhaseGrid.ion_initial|Начальные тяжёлые распределения и пучки
particle_kinetics|PhaseGrid.from_csv|Фазовый CSV
particle_kinetics|PhaseGrid.moments|Концентрация, энергия, скорость и температуры моментов
particle_kinetics|PhaseGrid.deposit|Проекция частиц и энергии
particle_kinetics|PhaseGrid.kick|Ускорение заданным Δv
particle_kinetics|PhaseGrid.equilibrium_moments|Начальное согласование дискретных моментов
particle_kinetics|PhaseGrid.ee_bgk|Приближённая e–e релаксация
particle_kinetics|stream_periodic|Периодический перенос
particle_kinetics|stream_open|Открытый перенос
collision_evolution|evolve_frozen|Положительный столкновительный подшаг
collision_evolution|uniformize|Пуассоновская сумма марковского генератора
level_selection|frozen_evolution|Экспонента населённостей и интегралы событий
level_selection|LevelSelection.__init__|Выбор активных уровней
atomic_core|Network.__init__|Ядра элементов и заряд стехиометрии
slab_transport|solve_transfer|Стационарный перенос света S_N
slab_transport|add_internal_beams|Внутренние листовые источники света
wall_boundaries|WallBoundaries.apply|Потоки и конверсия на z=0,L
spatial_sources|SpatialSources.__init__|Параметры и формы источников
spatial_sources|SpatialSources.apply|Инжекция частиц
spatial_sources|SpatialSources.photons|Включённые внутренние источники света
light_input|light_spectrum|Нормировка спектра и энергетика
light_input|particle_inputs|Концентрация/скорость источника → колонная инжекция
energy_input|energy_phase|Энергетический файл P(E)/EEPF → узлы
probe_spectrum|Probe.evaluate|Пропускание слабого зонда
reaction_data|equation|Единые полные реагенты и продукты
grid_validation|validate_phase_grid|Структурная проверка сеток
phase_grid_setup|resolve_heavy_grid|Автоматическая тяжёлая Emax
result_browser|ResultBrowser.options|Выбранные величины и масштабы осей
result_browser|ResultBrowser.select_plot|Явный выбор концентраций/уровней/реакций'''
OPERATORS+='\natomic_database|SourceCache.get|Проверка URL/хешей и чтение поставленных источников'


def code_map():
    rows=['| Процесс / роль | Файл и строка начала | Функция |','|---|---|---|']
    for entry in OPERATORS.splitlines():
        module,name,label=entry.split('|');path=ROOT/'src'/(module+'.py');node=ast.parse(path.read_text(encoding='utf-8'))
        for part in name.split('.'):
            node=next(c for c in node.body if isinstance(c,(ast.ClassDef,ast.FunctionDef,ast.AsyncFunctionDef)) and c.name==part)
        file='src/'+path.name
        rows.append(f'| {label} | [{file}:{node.lineno}](../{file}#L{node.lineno}) | `{name}` |')
    return '\n'.join(rows)+'\n'


def chemical_registry():
    header='# Зарегистрированные химические данные\n\nЭто все декларативные записи банка, включая альтернативы и оценки. Рабочая сеть конкретного запуска зависит от состава, флагов, доступных состояний и политики коэффициентов. Это не список всех возможных реакций природы.\n\n'
    sections=[];total=0
    for file in sorted((ROOT/'data/reactions').glob('*/reactions.json')):
        text=file.read_text(encoding='utf-8-sig');doc=json.loads(text);rows=doc.get('reactions',[]);total+=len(rows)
        section=[f'## {file.parent.name}',f'Записей: {len(rows)}. Выбор по умолчанию: `{json.dumps(doc.get("default_datasets",{}),ensure_ascii=False)}`. Участники-селекторы раскрыты в `selectors` этого JSON.','| ID и строка JSON | Реагенты → продукты | Тип | Коэффициенты | Статус / набор |','|---|---|---|---|---|']
        coefficients={}
        for name in ['cross_sections.csv','rate_coefficients.csv']:
            table=file.parent/name
            if table.exists():
                with table.open(encoding='utf-8-sig',newline='') as f:
                    reader=csv.DictReader(f,delimiter=';')
                    for row in reader:coefficients.setdefault(row['id'],{}).setdefault(name,[]).append(reader.line_num)
        for r in rows:
            line=next(i for i,t in enumerate(text.splitlines(),1) if ('"id": "'+r['id']+'"') in t)
            def side(values):return ' + '.join((str(v)+' ' if v!=1 else '')+k for k,v in values.items())
            formula=side(r['reactants'])+' → '+side(r['products']);links=[]
            for name,lines in coefficients.get(r.get('rate_id',r['id']),{}).items():
                p=(file.parent/name).relative_to(ROOT).as_posix();kind='σ(E)' if name=='cross_sections.csv' else 'k(T)'
                links.append(f'[{kind}, строка {min(lines)}](../{p}#L{min(lines)})')
            status=('оценка, по умолчанию выключена' if r.get('estimated') else 'запись банка')
            if r.get('alternative'):status+='; '+r['alternative']
            if r.get('radiative'):status+='; с фотоном'
            target=file.relative_to(ROOT).as_posix()
            section.append(f'| [`{r["id"]}`](../{target}#L{line}) | {formula} | `{r["kind"]}` | '+('; '.join(links) or 'см. rate_id и исходный JSON')+f' | {status} |')
        sections.append('\n\n'.join(section[:2])+'\n\n'+'\n'.join(section[2:]))
    return header+f'Всего зарегистрировано **{total} записей** в {len(sections)} группах. Источники и Q/области — в полных записях и CSV; таблица не назначает отсутствующие коэффициенты.\n\n'+'\n\n'.join(sections)+'\n'


def outputs():
    version=(ROOT/'VERSION.txt').read_text().strip();table=code_map()
    intro=f'# Где находится физика и химия — {version}\n\nСтроки начала функций получены из AST текущих файлов. После изменения кода выполните `python tools/build_reference.py`. Ссылки #L работают после загрузки проекта на GitHub; локально найдите указанную строку в редакторе.\n\n'
    physical=ROOT/'docs/PHYSICS_RU.md';base=physical.read_text(encoding='utf-8').split('<!-- GENERATED_CODE_MAP -->')[0].rstrip()
    appendix='\n\n<!-- GENERATED_CODE_MAP -->\n\n## Код со строками\n\n'+table
    return {ROOT/'docs/PHYSICS_CODE_MAP_RU.md':intro+table,physical:base+appendix,ROOT/'docs/REACTION_DATA_RU.md':chemical_registry()}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args();different=[]
    for file,value in outputs().items():
        if args.check:
            if not file.exists() or file.read_text(encoding='utf-8')!=value:different.append(str(file.relative_to(ROOT)))
        else:file.write_text(value,encoding='utf-8')
    if different:raise SystemExit('Обновите документацию: python tools/build_reference.py; '+', '.join(different))
    print('Ссылки кода и полный реестр химии актуальны.' if args.check else 'Обновлены физическая карта и реестр химии.')


if __name__=='__main__':main()
