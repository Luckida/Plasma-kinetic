# Источники версии 2.4

Числа в новых пакетах привязаны к конкретному процессу и диапазону. Поиск Google Scholar в среде подготовки не открылся; использованы первичные публикации, авторские копии и официальные базы. Ниже явно различены подключённые данные и найденные, но не подключённые процессы.

## Межкомпонентные процессы, использованные в расчёте

1. **Belyaev, A. K., Rodionov, D. S., Augustovičová, L., Soldán, P., Kraemer, W. P. (2015).** Full quantum study of non-radiative inelastic processes in lithium–helium ion–atom collisions. MNRAS 449, 3323–3332. [DOI](https://doi.org/10.1093/mnras/stv391). Использованы рисунки 5–7 и соответствие состояний таблицы 1. [Авторское воспроизведение статьи](https://www.mff.cuni.cz/cs/vedecka-rada/habilitacni-prace-dr-augustovicova.pdf). [Каталог исходных численных приложений CDS](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/MNRAS/449/3323) — найден, но в эту версию числа оттуда не получены. В пакете — извлечённые векторные кривые, не подмена их машинными таблицами CDS.

2. **García Muñoz, A. (2025).** Modeling helium in exoplanet atmospheres. A revised network with photoelectron-driven processes. A&A 698, A199. [DOI](https://doi.org/10.1051/0004-6361/202555145), [авторская версия](https://arxiv.org/abs/2505.12148). Таблицы A.5–A.7: тепловые коэффициенты He*–H, HeH⁺ и перезарядки. Разбиение Пеннинг/ассоциация 90:10 — допущение автора. Указанные там первичные работы Movre & Meyer (1997), Movre et al. (1994), Loreau et al. (2018), De Fazio (2014), Black (1978) не выдаются за независимо оцифрованные здесь таблицы.

3. **Sixt, T., Chung, T., Stienkemeier, F., Dulitz, K. (2023).** Symmetry Dependence of the Continuum Coupling in the Chemi-ionization of Li(2²S₁/₂) by He(2³S₁,2³P_J). J. Phys. Chem. A 127, 4407–4414. [DOI](https://doi.org/10.1021/acs.jpca.3c00431), [полный текст](https://pmc.ncbi.nlm.nih.gov/articles/PMC10226102/). Таблица 2 даёт якорные полу-классические коэффициенты при 44 мэВ. Продолжение к другим энергиям σ∝E^(-1/3) — **модель версии 2.4**, явно помеченная как оценочная.

4. **Dulitz, K., Sixt, T., Guan, J., Grzesiak, J., Debatin, M., Stienkemeier, F. (2020).** Suppression of Penning ionization by orbital angular momentum conservation. Phys. Rev. A 102, 022818. [DOI](https://doi.org/10.1103/PhysRevA.102.022818), [авторская версия](https://arxiv.org/abs/1909.01694). Обоснование необходимости различать спиновые/орбитальные каналы; модель скалярных невыстроенных состояний не воспроизводит квантовую подготовку из этой работы.

5. **Grzesiak, J., Momose, T., Stienkemeier, F., Mudrich, M., Dulitz, K. (2019).** Penning collisions between supersonically expanded metastable helium atoms and laser-cooled lithium atoms. J. Chem. Phys. 150, 034201. [DOI](https://doi.org/10.1063/1.5063709), [авторская версия](https://arxiv.org/abs/1901.02482). Использовано при оценке неполноты: обнаруженная малая молекулярная ветвь не заменена произвольной константой для HeLi⁺.

6. **Galli, D., Palla, F. (1998).** The Chemistry of the Early Universe. A&A 335, 403–420. [Авторская версия](https://arxiv.org/abs/astro-ph/9803315). Таблица 4, Li7/Li8: нерадиационная и радиационная перезарядка Li+H⁺. Для Li7 первичный источник: Kimura, Dutta & Shimakura (1994), [DOI](https://doi.org/10.1086/174418). Для Li8 — Stancil & Zygelman (1996), как процитировано в таблице. Введённый в коде защитный диапазон 100–20 000 K не является новой оценкой погрешности исходных фит-функций.

7. **Roberge, W., Dalgarno, A. (1982).** ApJ 255, 489. [DOI](https://doi.org/10.1086/159849). Сохранённый альтернативный набор He(2³S)+H. Интерполяция ассоциативной доли 90/T — прежнее приближение этого кода.

## Атомные и фотонные данные общего ядра

8. **OPEN-ADAS.** [Каталог Li](https://open.adas.ac.uk/element/li), [H](https://open.adas.ac.uk/element/h), [He](https://open.adas.ac.uk/element/he). Уровни, A и термические Υ; используются закреплённые дескрипторами ADF04 LS-наборы. Υ не превращаются в σ(E). Первичные ссылки и контрольные суммы указаны в `data/elements/*/atom.json` и `manifest.json` подготовленных наборов.

9. **CCC, Curtin University.** [Официальная база](https://atom.curtin.edu.au/CCC-WWW/index.html). Электронные σ(E) H I, He I/II, Li I; конкретные URL и хеши — в каталогах каждого атома. Для Li/Na исходная работа Bray (1994), [DOI](https://doi.org/10.1103/PhysRevA.49.1066); для H — Bray & Stelbovics (1992), [DOI](https://doi.org/10.1103/PhysRevA.46.6995); для He — [DOI](https://doi.org/10.1103/PhysRevA.52.1279).

10. **Verner, D. A., Ferland, G. J., Korista, K. T., Yakovlev, D. G. (1996).** Atomic Data for Astrophysics. II. New Analytic Fits for Photoionization Cross Sections of Atoms and Ions. ApJ 465, 487. [DOI](https://doi.org/10.1086/177435). Фото-ионизация основных состояний; используется существующим PhotoBank.

11. **Hummer, D. G., Storey, P. J. (1998).** Recombination of helium-like ions. I. Photoionization cross-sections and total recombination and cooling coefficients for atomic helium. MNRAS 297, 1073. [DOI](https://doi.org/10.1046/j.1365-8711.1998.2970041073.x). Фото-ионизация He I и обратная RR по Милну. Таблицы получаются из [официальных данных Cloudy](https://gitlab.nublado.org/cloudy/cloudy/-/tree/master/data/sh98_he1/pi); лицензия сохранена в `data/radiation/Cloudy_license.txt`.

12. **UMIST Database for Astrochemistry, Rate22.** [База](https://umistdatabase.net/), [публикация Rate22](https://doi.org/10.1051/0004-6361/202346908). Сохранённые исходные молекулярные реакции H/He. DOI первичных измерений и статус коэффициентов остаются в каждой строке `rate_coefficients.csv`.

13. **He, Brinkmann, Kemaneci & Gibson (2025).** Plasma Sources Sci. Technol. 34, 095018. [DOI](https://doi.org/10.1088/1361-6595/adffe6). Сохранённые исходные реакции He*+He*, He₂⁺ и диссоциативная рекомбинация. Это прежние данные архива; они не представляют новое полное состояние-разрешённое описание.

## Найдены, но не подключены как готовые численные каналы

14. **Zámečníková, M., Kraemer, W. P., Soldán, P. (2018).** Radiative Charge Transfer between Metastable Helium and Lithium Cations. ApJ 867, 157. [DOI](https://doi.org/10.3847/1538-4357/aae64f). Нужна температурная таблица/фит. Одного диапазона величин из аннотации недостаточно для k(T).

15. **Augustovičová, L., Kraemer, W. P., Soldán, P. (2014).** Depopulation of metastable helium by radiative association with hydrogen and lithium ions. ApJ 782, 46. [DOI](https://doi.org/10.1088/0004-637X/782/1/46). Найдены суммарные фиты и энергии молекулярных состояний. Для корректного подключения нужны отдельные триплетные молекулярные носители, распределение по v,J и разрушение/фотопроцессы; не отождествлены с основным HeH⁺.

16. **Loreau, J., Ryabchenko, S., Burgos, J. M. M., Vaeck, N. (2018).** Charge transfer in low-energy collisions of H with He⁺ and H⁺ with He in excited states. J. Phys. B. [DOI](https://doi.org/10.1088/1361-6455/aab425). [Авторская версия](https://arxiv.org/abs/1801.06821). Используются суммарные коэффициенты García Muñoz (2025), с явным указанием происхождения, экстраполяции и приближённого основного продукта H. Обратные реакции из суммарного коэффициента не строятся.

17. **Plowman et al. (2025).** Convergent close-coupling approach to excitation in proton collisions with helium at intermediate energies. Phys. Rev. A 112. [DOI](https://doi.org/10.1103/y3hw-kyc4). Обнаруженный новый набор CCC для возбуждения He(основное)+p относится к кэВ-энергиям. Не экстраполирован к тепловым столкновениям эВ-диапазона данного примера.

Для поиска работ в Google Scholar можно использовать их точные названия или DOI. Прямой доступ к Scholar был недоступен; факт поиска первоисточников не означает, что прочитаны или подключены все когда-либо опубликованные данные.


## Определения распределений и ввод в версии 2.8

- G. J. M. Hagelaar. *Brief documentation of BOLSIG+ version 03/2016*. [Официальное руководство](https://www.bolsig.laplace.univ-tlse.fr/wp-content/uploads/2016/03/bolsigdoc0316.pdf). Разделы Rate coefficients и List of symbols: конвенция EEPF с размерностью эВ⁻³ᐟ² и интегрирование электронных сечений. Использовано для формата входа; внешний BOLSIG+ не запускается.
- OpenMC. *openmc.stats.Maxwell*. [Официальная документация](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.stats.Maxwell.html). Определение максвелловской энергетической плотности через √E exp(−E/θ). OpenMC не является зависимостью модели.

Численные атомные и химические коэффициенты в версии 2.8 не добавлялись и не менялись. Общий справочник показывает все доступные для состава записи существующего банка, их первичные ссылки, диапазоны и ограничения. Предлагаемые примеры P(E) являются аналитическими входными формами, а не новыми измерениями.
