# Где находится физика и химия — 2.8.2

Строки начала функций получены из AST текущих файлов. После изменения кода выполните `python tools/build_reference.py`. Ссылки #L работают после загрузки проекта на GitHub; локально найдите указанную строку в редакторе.

| Процесс / роль | Файл и строка начала | Функция |
|---|---|---|
| Состав и начальные распределения | [src/mixture_kinetics.py:75](../src/mixture_kinetics.py#L75) | `MixtureSlab.__init__` |
| Флаги процессов | [src/mixture_kinetics.py:222](../src/mixture_kinetics.py#L222) | `MixtureSlab.enabled` |
| Атомные электронные и фотонные каналы | [src/mixture_kinetics.py:286](../src/mixture_kinetics.py#L286) | `MixtureSlab._atomic_channels` |
| Стехиометрия общей сети | [src/mixture_kinetics.py:343](../src/mixture_kinetics.py#L343) | `MixtureSlab._network` |
| Парная геометрия и рекомбинационные ядра | [src/mixture_kinetics.py:352](../src/mixture_kinetics.py#L352) | `MixtureSlab._pairs` |
| Частоты атомных процессов и химические события | [src/mixture_kinetics.py:463](../src/mixture_kinetics.py#L463) | `MixtureSlab.rates` |
| Недостающие атомные процессы и диапазоны | [src/mixture_kinetics.py:484](../src/mixture_kinetics.py#L484) | `MixtureSlab._check_data` |
| Предписанное поле и экспериментальный Poisson | [src/mixture_kinetics.py:392](../src/mixture_kinetics.py#L392) | `MixtureSlab.field_E` |
| Парные RR/SR скорости | [src/mixture_kinetics.py:402](../src/mixture_kinetics.py#L402) | `MixtureSlab.rr_values` |
| Общее поле излучения | [src/mixture_kinetics.py:419](../src/mixture_kinetics.py#L419) | `MixtureSlab.field` |
| Поглощение и излучение активного базиса | [src/mixture_kinetics.py:449](../src/mixture_kinetics.py#L449) | `MixtureSlab.optical_matrices` |
| e–neutral, e–ion и BGK столкновения | [src/mixture_kinetics.py:541](../src/mixture_kinetics.py#L541) | `MixtureSlab._binary` |
| Связанный шаг населённостей, EEDF и химии | [src/mixture_kinetics.py:579](../src/mixture_kinetics.py#L579) | `MixtureSlab._local` |
| Ускорение и пространственный перенос | [src/mixture_kinetics.py:643](../src/mixture_kinetics.py#L643) | `MixtureSlab._transport` |
| Расщепление и контроль шага | [src/mixture_kinetics.py:683](../src/mixture_kinetics.py#L683) | `MixtureSlab.solve` |
| Полная энергия частиц и уровней | [src/mixture_kinetics.py:389](../src/mixture_kinetics.py#L389) | `MixtureSlab.total_energy` |
| Балансы, события и сохранённый срез | [src/mixture_kinetics.py:667](../src/mixture_kinetics.py#L667) | `MixtureSlab.snapshot` |
| Согласование спектрального профиля | [src/mixture_kinetics.py:46](../src/mixture_kinetics.py#L46) | `conservative_profile` |
| Сумма σv по EEDF | [src/kinetic_kernels.py:102](../src/kinetic_kernels.py#L102) | `impact_rates` |
| Электронные потери и продукты неупругих событий | [src/kinetic_kernels.py:112](../src/kinetic_kernels.py#L112) | `impact_events` |
| Относительные электрон–тяжёлые скорости | [src/kinetic_kernels.py:143](../src/kinetic_kernels.py#L143) | `pair_geometry` |
| Обратный PI оператор Милна | [src/kinetic_kernels.py:159](../src/kinetic_kernels.py#L159) | `build_rr` |
| Спонтанная рекомбинация | [src/kinetic_kernels.py:179](../src/kinetic_kernels.py#L179) | `rr_emission` |
| Стимулированная рекомбинация | [src/kinetic_kernels.py:208](../src/kinetic_kernels.py#L208) | `rr_stimulated` |
| Тяжёлые продукты рекомбинации | [src/kinetic_kernels.py:231](../src/kinetic_kernels.py#L231) | `rr_heavy_products` |
| Экранированное кулоновское транспортное сечение | [src/kinetic_kernels.py:257](../src/kinetic_kernels.py#L257) | `momentum_transfer_rate` |
| Генератор парных упругих столкновений | [src/kinetic_kernels.py:304](../src/kinetic_kernels.py#L304) | `binary_generator_rates` |
| Выбор фото-сечений | [src/mixture_photo.py:64](../src/mixture_photo.py#L64) | `PhotoBank.__init__` |
| Кулоновские фото-интегралы | [src/mixture_photo.py:27](../src/mixture_photo.py#L27) | `hydrogen_sigma` |
| Фото-фиты Verner | [src/mixture_photo.py:54](../src/mixture_photo.py#L54) | `verner_sigma` |
| Молекулярные носители и энергетические якоря | [src/mixture_chemistry.py:29](../src/mixture_chemistry.py#L29) | `append_components` |
| Все используемые законы k(T) | [src/mixture_chemistry.py:77](../src/mixture_chemistry.py#L77) | `rate_fit` |
| Чтение сети и выбор альтернатив | [src/mixture_chemistry.py:107](../src/mixture_chemistry.py#L107) | `Chemistry.__init__` |
| Парное σv и моментное замыкание | [src/mixture_chemistry.py:234](../src/mixture_chemistry.py#L234) | `Chemistry.kernel` |
| Скорости химических событий | [src/mixture_chemistry.py:275](../src/mixture_chemistry.py#L275) | `Chemistry.event_rates` |
| Стоки и продукты химии | [src/mixture_chemistry.py:284](../src/mixture_chemistry.py#L284) | `Chemistry.evolve` |
| Фотоны радиационной перезарядки | [src/mixture_chemistry.py:325](../src/mixture_chemistry.py#L325) | `Chemistry.photon_emission` |
| Относительная энергия тяжёлых пар | [src/chemistry_kernels.py:15](../src/chemistry_kernels.py#L15) | `geometry` |
| Кинематика продуктов в координатах Якоби | [src/chemistry_kernels.py:27](../src/chemistry_kernels.py#L27) | `products` |
| Энергетическая и угловая сетки | [src/particle_kinetics.py:27](../src/particle_kinetics.py#L27) | `PhaseGrid.__init__` |
| Начальная EEDF | [src/particle_kinetics.py:44](../src/particle_kinetics.py#L44) | `PhaseGrid.initial` |
| Начальные тяжёлые распределения и пучки | [src/particle_kinetics.py:69](../src/particle_kinetics.py#L69) | `PhaseGrid.ion_initial` |
| Фазовый CSV | [src/particle_kinetics.py:95](../src/particle_kinetics.py#L95) | `PhaseGrid.from_csv` |
| Концентрация, энергия, скорость и температуры моментов | [src/particle_kinetics.py:116](../src/particle_kinetics.py#L116) | `PhaseGrid.moments` |
| Проекция частиц и энергии | [src/particle_kinetics.py:127](../src/particle_kinetics.py#L127) | `PhaseGrid.deposit` |
| Ускорение заданным Δv | [src/particle_kinetics.py:161](../src/particle_kinetics.py#L161) | `PhaseGrid.kick` |
| Начальное согласование дискретных моментов | [src/particle_kinetics.py:177](../src/particle_kinetics.py#L177) | `PhaseGrid.equilibrium_moments` |
| Приближённая e–e релаксация | [src/particle_kinetics.py:212](../src/particle_kinetics.py#L212) | `PhaseGrid.ee_bgk` |
| Периодический перенос | [src/particle_kinetics.py:250](../src/particle_kinetics.py#L250) | `stream_periodic` |
| Открытый перенос | [src/particle_kinetics.py:265](../src/particle_kinetics.py#L265) | `stream_open` |
| Положительный столкновительный подшаг | [src/collision_evolution.py:14](../src/collision_evolution.py#L14) | `evolve_frozen` |
| Пуассоновская сумма марковского генератора | [src/collision_evolution.py:33](../src/collision_evolution.py#L33) | `uniformize` |
| Экспонента населённостей и интегралы событий | [src/level_selection.py:31](../src/level_selection.py#L31) | `frozen_evolution` |
| Выбор активных уровней | [src/level_selection.py:51](../src/level_selection.py#L51) | `LevelSelection.__init__` |
| Ядра элементов и заряд стехиометрии | [src/atomic_core.py:130](../src/atomic_core.py#L130) | `Network.__init__` |
| Стационарный перенос света S_N | [src/slab_transport.py:40](../src/slab_transport.py#L40) | `solve_transfer` |
| Внутренние листовые источники света | [src/slab_transport.py:106](../src/slab_transport.py#L106) | `add_internal_beams` |
| Потоки и конверсия на z=0,L | [src/wall_boundaries.py:61](../src/wall_boundaries.py#L61) | `WallBoundaries.apply` |
| Параметры и формы источников | [src/spatial_sources.py:43](../src/spatial_sources.py#L43) | `SpatialSources.__init__` |
| Инжекция частиц | [src/spatial_sources.py:119](../src/spatial_sources.py#L119) | `SpatialSources.apply` |
| Включённые внутренние источники света | [src/spatial_sources.py:129](../src/spatial_sources.py#L129) | `SpatialSources.photons` |
| Нормировка спектра и энергетика | [src/light_input.py:28](../src/light_input.py#L28) | `light_spectrum` |
| Концентрация/скорость источника → колонная инжекция | [src/light_input.py:53](../src/light_input.py#L53) | `particle_inputs` |
| Энергетический файл P(E)/EEPF → узлы | [src/energy_input.py:40](../src/energy_input.py#L40) | `energy_phase` |
| Пропускание слабого зонда | [src/probe_spectrum.py:61](../src/probe_spectrum.py#L61) | `Probe.evaluate` |
| Единые полные реагенты и продукты | [src/reaction_data.py:23](../src/reaction_data.py#L23) | `equation` |
| Структурная проверка сеток | [src/grid_validation.py:10](../src/grid_validation.py#L10) | `validate_phase_grid` |
| Автоматическая тяжёлая Emax | [src/phase_grid_setup.py:15](../src/phase_grid_setup.py#L15) | `resolve_heavy_grid` |
| Выбранные величины и масштабы осей | [src/result_browser.py:135](../src/result_browser.py#L135) | `ResultBrowser.options` |
| Явный выбор концентраций/уровней/реакций | [src/result_browser.py:122](../src/result_browser.py#L122) | `ResultBrowser.select_plot` |
| Проверка URL/хешей и чтение поставленных источников | [src/atomic_database.py:98](../src/atomic_database.py#L98) | `SourceCache.get` |
