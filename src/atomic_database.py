"""Named-species data preparation, independent of the helium time integrator.

Reviewed source choices live in atomic_sources.json. New species without a
reviewed importer get a candidate list and explicit gaps, never helium data.
Downloaded numerical tables are personal caches, not part of the release ZIP.
"""
import csv
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import re
import tempfile
from urllib.parse import urlencode, urljoin, urlsplit
from urllib.request import Request, urlopen
import zipfile
import numpy as np
from scipy import constants as C
from atomic_core import Species, selected_species
from progress import Progress

ROOT = Path(__file__).resolve().parent.parent
HC_EV_CM = C.h * C.c / C.e * 100
BOHR_CM2 = (C.physical_constants['Bohr radius'][0] * 100) ** 2
LEVEL_FIELDS = ['id', 'species', 'configuration', 'n', 'multiplicity', 'L', 'g', 'energy_cm1', 'energy_reference', 'mass_u', 'source_index', 'label', 'source']
LINE_FIELDS = ['id', 'lower', 'upper', 'A_s', 'wavelength_nm', 'wavelength_medium', 'source']
RATE_FIELDS = ['id', 'kind', 'initial', 'final', 'T_K', 'value', 'units', 'source']
XS_FIELDS = ['id', 'kind', 'initial', 'final', 'energy_eV', 'sigma_cm2', 'source']
CATALOG_FIELDS = ['id', 'kind', 'initial', 'final', 'energy_min_eV', 'energy_max_eV', 'model_gap_eV', 'source_first_energy_eV', 'source', 'source_sha256']


def write_csv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter=';')
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def number(text):
    return float(re.sub(r'(?<=[0-9.])([+-]\d+)$', r'e\1', text.replace('D', 'E')))


def science_hash(raw, mode):
    if mode == 'adas_science':
        text = raw.decode('utf-8')
        raw = re.split(r'(?m)^C\s+OPEN-ADAS\s*$', text)[0].rstrip().encode()
    elif mode != 'bytes':
        raise ValueError('Неизвестный hash_mode')
    return hashlib.sha256(raw).hexdigest()


def catalogue_selection_hash(registry,species_key):
    """Unrelated new elements must not invalidate a verified existing pack."""
    entry=registry['species'].get(species_key)
    keys=[entry.get(k) for k in ('adf04','electron_cross_sections')] if entry else []
    if entry:keys.extend(entry.get('native_csv',{}).values())
    selected=dict(species=entry,sources={key:registry['sources'][key] for key in keys if key in registry.get('sources',{})})
    return hashlib.sha256(json.dumps(selected,sort_keys=True,ensure_ascii=False).encode()).hexdigest()


def load_catalogue(root,options):
    """One atom.json per element; legacy single-file registries still read."""
    if not options.get('element_directory'):
        return json.loads((Path(root)/options.get('catalogue_file','atomic_sources.json')).read_text(encoding='utf-8-sig'))
    registry=dict(schema='AtomicSources-1',sources={},species={},ground_stages={},descriptor_files_sha256={})
    folder=Path(root)/options['element_directory']
    for path in sorted(folder.glob('*/atom.json')):
        item=json.loads(path.read_text(encoding='utf-8-sig'))
        if item.get('schema')!='ElementData-1':raise ValueError('Нужен ElementData-1: '+str(path))
        for key in ('sources','species'):
            if set(registry[key])&set(item[key]):raise ValueError('Повтор '+key+' в '+str(path))
            registry[key].update(item[key])
        registry['ground_stages'][item['element']]=item['ground_stages']
        registry['descriptor_files_sha256'][item['element']]=hashlib.sha256(path.read_bytes()).hexdigest()
    if not registry['species']:raise ValueError('Нет дескрипторов элементов в '+str(folder))
    return registry


def pack_directory(output,item,nmax,options):
    return Path(output)/(item.element if options.get('element_directory') else '')/f'{item.key}_n{nmax}'


class SourceCache:
    def __init__(self, directory, offline=False, progress=None, bundled_directory=None):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.progress = progress or Progress()
        self.receipts = {}
        self.bundled_directory = Path(bundled_directory) if bundled_directory is not None else None

    def get(self, name, info):
        url = info['url']
        if urlsplit(url).scheme != 'https':
            raise ValueError('Для атомных источников нужен https')
        safe = re.sub(r'[^a-zA-Z0-9_.-]', '_', name)
        path = self.directory / (hashlib.sha256(url.encode()).hexdigest()[:16] + '_' + safe)
        receipt_path = path.with_suffix(path.suffix + '.json')
        if not path.exists() and self.bundled_directory is not None:
            bundled = self.bundled_directory / path.name
            if bundled.is_file():
                path = bundled
                receipt_path = path.with_suffix(path.suffix + '.json')
                if not receipt_path.is_file():
                    raise ValueError('Нет квитанции происхождения поставленного источника ' + name)
        if path.exists():
            raw = path.read_bytes()
            if receipt_path.exists():
                receipt = json.loads(receipt_path.read_text(encoding='utf-8-sig'))
                if receipt['url'] != url or receipt['sha256_bytes'] != hashlib.sha256(raw).hexdigest():
                    raise ValueError('Повреждён кэш источника ' + name)
        else:
            if self.offline:
                raise OSError('Нет локального источника ' + name + ' (offline=true)')
            self.progress.module('Загрузка атомных данных', name, force=True)
            request = Request(url, headers={'User-Agent': 'AtomicKinetics/1.7 research data loader'})
            chunks, size = [], 0
            with urlopen(request, timeout=30) as response:
                while True:
                    block = response.read(131072)
                    if not block:
                        break
                    chunks.append(block)
                    size += len(block)
                    if size > 40_000_000:
                        raise ValueError('Атомный источник больше 40 MB: ' + name)
                    self.progress.module('Загрузка атомных данных', f'{name}: {size / 1e6:.2f} MB')
            raw = b''.join(chunks)
            if not raw:
                raise ValueError('Пустой атомный источник: ' + name)
        digest = science_hash(raw, info.get('hash_mode', 'bytes'))
        if info.get('sha256') and info['sha256'] != digest:
            raise ValueError('Изменился проверенный источник ' + name + '; новая версия не подключена автоматически')
        receipt = dict(url=url, sha256_bytes=hashlib.sha256(raw).hexdigest(), sha256_science=digest,
                       source_version='pinned' if info.get('sha256') else 'recorded_first_use_unreviewed',
                       reference=info.get('reference', ''), retrieved_utc=datetime.now(timezone.utc).isoformat())
        if not path.exists():
            temp = path.with_suffix(path.suffix + '.part')
            temp.write_bytes(raw)
            temp.replace(path)
            write_json(receipt_path, receipt)
        elif receipt_path.exists():
            receipt = json.loads(receipt_path.read_text(encoding='utf-8-sig'))
        else:
            write_json(receipt_path, receipt)
        self.receipts[name] = receipt
        return raw


def parse_adf04_ls3(raw, species, mass_u, source, nmax):
    """Strict type-3 LS importer. Extra infinite-T value is NOT another T node.

    Unsupported transition codes are counted explicitly. Parent energies are
    separate, avoiding the false identification of different ionic grounds.
    """
    text = raw.decode('utf-8')
    lines = text.splitlines()
    header = re.match(r'^\s*([A-Z][a-z]?)\+\s*(\d+)\s+(\d+)\s+(\d+)\s+(.+)$', lines[0], re.I)
    if not header or header.group(1).casefold() != species.element.casefold() or int(header.group(2)) != species.charge:
        raise ValueError('ADF04 не соответствует запрошенному атому/иону')
    parents = [dict(parent=k+1, threshold_cm1=float(m.group(1)), term=m.group(2))
               for k, m in enumerate(re.finditer(r'([\d.]+)\(([^)]+)\)', header.group(5)))]
    if not parents:
        raise ValueError('Нет родительского порога ADF04')
    stop = next(i for i, line in enumerate(lines[1:], 1) if line.strip() == '-1')
    pattern = re.compile(r'^\s*(\d+)\s+(.+?)\s+\((\d+)\)(\d+)\(\s*([\d.]+)\)\s+([\d.]+)')
    levels, by_index = [], {}
    for line in lines[1:stop]:
        m = pattern.match(line)
        if not m:
            raise ValueError('Неизвестная строка LS-уровня: ' + line)
        index, configuration, mult, angular, pseudo_j, energy = m.groups()
        configuration_for_n=re.sub(r'\([^)]*\)','',configuration)
        orbitals = re.findall(r'(\d+)[SPDFGHIJK]\d*', configuration_for_n,re.I)
        if not orbitals:
            raise ValueError('Неизвестная конфигурация: ' + configuration)
        n, mult, angular = max(map(int, orbitals)), int(mult), int(angular)
        weight = mult * (2 * angular + 1)
        if abs(2 * float(pseudo_j) + 1 - weight) > 1e-9:
            raise ValueError('Набор не является LS-термами; нужен отдельный IC/J-импортёр')
        if n > nmax:
            continue
        ident = f'{species.key}:L{int(index):03d}'
        row = dict(id=ident, species=species.label, configuration=configuration.strip(), n=n,
                   multiplicity=mult, L=angular, g=weight, energy_cm1=float(energy),
                   energy_reference='ground_of_this_ion', mass_u=mass_u, source_index=int(index),
                   label=f'{species.label} {configuration.strip()} {mult}{"SPDFGHIJK"[angular]}', source=source)
        levels.append(row)
        by_index[int(index)] = row
    temp_line = lines[stop + 1].split()
    if int(temp_line[1]) != 3:
        raise ValueError('Импортёр поддерживает ADF04 ITYP=3 (upsilon)')
    temperatures = [number(value) for value in temp_line[2:]]
    if not temperatures or min(temperatures) <= 0 or np.any(np.diff(temperatures) <= 0):
        raise ValueError('Неверная температурная сетка ADF04')
    rates, radiative, omitted = [], [], {}
    for line in lines[stop+2:]:
        if line.strip().startswith('-1'):
            break
        cols = line.split()
        if not cols:
            continue
        if not cols[0].isdigit():
            omitted[cols[0]] = omitted.get(cols[0], 0) + 1
            continue
        upper, lower = int(cols[0]), int(cols[1])
        if upper not in by_index or lower not in by_index:
            continue
        up, low = by_index[upper], by_index[lower]
        aval = number(cols[2])
        values = [number(v) for v in cols[3:]]
        if len(values) not in (len(temperatures), len(temperatures)+1) or min(values) < 0:
            raise ValueError('Неверное число/знак upsilon ADF04')
        ident = f'CE:{low["id"]}:{up["id"]}'
        for temp, value in zip(temperatures, values):
            rates.append(dict(id=ident, kind='upsilon', initial=low['id'], final=up['id'],
                              T_K=temp, value=value, units='dimensionless_Maxwellian_effective_collision_strength', source=source))
        if aval > 1e-20:
            gap = up['energy_cm1'] - low['energy_cm1']
            if gap <= 0:
                raise ValueError('Положительное A для неположительной разности энергий')
            radiative.append(dict(id=f'A:{up["id"]}:{low["id"]}', lower=low['id'], upper=up['id'],
                                  A_s=aval, wavelength_nm=1e7/gap, wavelength_medium='vacuum_LS_centroid', source=source))
    return levels, rates, radiative, dict(parents=parents, unsupported_transition_codes=omitted,
        infinite_temperature_upsilon='extra final value, when present, not exported as a finite-temperature node',
        temperature_min_K=min(temperatures), temperature_max_K=max(temperatures))


def ccc_label(level, style):
    if style == 'helium':
        return ('s' if level['multiplicity'] == 1 else 't') + str(level['n']) + 'SPDFGHI'[level['L']]
    return str(level['n']) + 'SPDFGHI'[level['L']]


def ccc_zip_tables(raw, levels, style, source):
    labels = {ccc_label(row, style): row for row in levels}
    rows, catalogue = [], []
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        if sum(info.file_size for info in archive.infolist()) > 200_000_000:
            raise ValueError('Неожиданный размер распакованного CCC')
        for info in archive.infolist():
            if info.is_dir():
                continue
            with archive.open(info) as stream:
                header = stream.readline().decode('utf-8')
            match = re.search(r'#\s*(\S+)\s*->\s*([^;\s]+)\s*;\s*Units:\s*(\S+)', header)
            if not match:
                continue
            initial, final, units = match.groups()
            if initial not in labels or (final not in labels and final not in ('TICS', 'MTCS')):
                continue
            if units not in ('a0^2','cm^2'):
                raise ValueError('Неизвестные единицы CCC: ' + units)
            scale=BOHR_CM2 if units=='a0^2' else 1.
            values = []
            for line in archive.read(info).decode().splitlines():
                if not line.strip() or line.lstrip().startswith('#'):
                    continue
                cols = line.split()
                values.append((float(cols[0]), float(cols[1]) * scale))
            values = np.array(values)
            if len(values) < 2 or not np.all(np.isfinite(values)) or np.any(np.diff(values[:, 0]) <= 0) or np.any(values < 0):
                raise ValueError('Неверная CCC-таблица ' + info.filename)
            lo = labels[initial]
            if final == 'TICS':
                kind, target, gap = 'electron_ionization', f'{lo["species"]}:parent_unresolved', ''
            elif final == 'MTCS':
                kind, target, gap = 'electron_momentum_transfer', lo['id'], 0.
            else:
                up = labels[final]
                gap = (up['energy_cm1'] - lo['energy_cm1']) * HC_EV_CM
                kind = 'electron_elastic' if final == initial else ('electron_excitation' if gap >= 0 else 'electron_deexcitation')
                target = up['id']
            ident = f'{kind}:{lo["id"]}:{target}'
            catalogue.append(dict(id=ident, kind=kind, initial=lo['id'], final=target,
                                  energy_min_eV=float(values[0, 0]), energy_max_eV=float(values[-1, 0]),
                                  model_gap_eV=gap, source_first_energy_eV=float(values[0, 0]),
                                  source=source + ' :: ' + info.filename, source_sha256=hashlib.sha256(raw).hexdigest()))
            rows.extend(dict(id=ident, kind=kind, initial=lo['id'], final=target, energy_eV=float(e), sigma_cm2=float(s), source=source + ' :: ' + info.filename) for e, s in values)
    if len({row['id'] for row in catalogue}) != len(catalogue):
        raise ValueError('Повтор канала CCC: нельзя автоматически выбрать набор')
    return rows, catalogue


class FormInputs(HTMLParser):
    def __init__(self):
        super().__init__()
        self.inputs = {'low': {}, 'upp': {}}

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag.lower() == 'input' and attrs.get('name') in self.inputs:
            value = attrs.get('value', '')
            self.inputs[attrs['name']][value.split('_')[0].upper()] = value


def alkali_tables(cache, levels, element='Li'):
    if element not in ('Li','Na'):raise ValueError('Этот CCC-адаптер проверен только для нейтральных Li/Na')
    query = urlencode(dict(element=element, charge='I', process='Electron Impact'))
    raw = cache.get('ccc_'+element+'_form.html', dict(url='https://atom.curtin.edu.au/cgi-bin/c4nick.1?' + query))
    form = FormInputs()
    form.feed(raw.decode('utf-8'))
    if not form.inputs['low'] or not form.inputs['upp']:
        raise ValueError('Изменился опубликованный интерфейс CCC '+element)
    rows, catalogue, gaps = [], [], []
    for lo in levels:
        lab = ccc_label(lo, 'lithium')
        if lab not in form.inputs['low']:
            gaps.append(f'CCC: нет начального {lo["label"]}')
            continue
        targets = [(ccc_label(up, 'lithium'), up) for up in levels if up['energy_cm1'] >= lo['energy_cm1']]
        targets += [('TICS', None), ('MTCS', lo)]
        for target, up in targets:
            if target not in form.inputs['upp']:
                gaps.append(f'CCC: нет конечного {target} для {lo["label"]}')
                continue
            q = urlencode(dict(element=element, charge='I', process='ELECTRON_IMPACT', low=form.inputs['low'][lab],
                               upp=form.inputs['upp'][target], data_type='cross', cs_unit='cm2'))
            url = 'https://atom.curtin.edu.au/cgi-bin/c4nick.in?' + q
            name = f'ccc_{element}_{lab}_{target}.html'
            try:
                response = cache.get(name, dict(url=url, reference='I. Bray, Phys. Rev. A 49, 1066 (1994)'))
            except OSError as exc:
                gaps.append(f'CCC download error {lab} → {target}: {type(exc).__name__}: {exc}')
                continue
            text = response.decode('utf-8')
            pairs = re.findall(r'<TD\s+ALIGN=RIGHT>\s*([\d.+Ee-]+)<TD><TD><TD\s+ALIGN=RIGHT>\s*([\d.+Ee-]+)', text, re.I)
            if not pairs:
                gaps.append(f'CCC: пустая таблица {lab} → {target}')
                continue
            if not re.search(r'Cross section in units of cm<SUP>2</SUP>', text, re.I):
                raise ValueError('Не подтверждены единицы cm2 CCC Li')
            values = np.array(pairs, dtype=float)
            if len(values) < 2 or not np.all(np.isfinite(values)) or np.any(np.diff(values[:, 0]) <= 0) or np.any(values < 0):
                gaps.append(f'CCC rejected nonphysical/nonmonotone table {element} {lab} -> {target}; no clipping or invented values')
                continue
            if target == 'TICS':
                kind, final, gap = 'electron_ionization', element+' I:parent_unresolved', ''
            elif target == 'MTCS':
                kind, final, gap = 'electron_momentum_transfer', lo['id'], 0.
            else:
                kind = 'electron_elastic' if up['id'] == lo['id'] else 'electron_excitation'
                final, gap = up['id'], (up['energy_cm1'] - lo['energy_cm1']) * HC_EV_CM
            ident = f'{kind}:{lo["id"]}:{final}'
            catalogue.append(dict(id=ident, kind=kind, initial=lo['id'], final=final,
                                  energy_min_eV=float(values[0, 0]), energy_max_eV=float(values[-1, 0]), model_gap_eV=gap,
                                  source_first_energy_eV=float(values[0, 0]), source=url,
                                  source_sha256=cache.receipts[name]['sha256_bytes']))
            rows.extend(dict(id=ident, kind=kind, initial=lo['id'], final=final, energy_eV=float(e), sigma_cm2=float(s), source=url) for e, s in values)
    return rows, catalogue, gaps


# Compatibility for third-party callers of the original Li data importer.
lithium_tables=alkali_tables


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        href = dict(attrs).get('href', '')
        if tag.lower() == 'a' and '/detail/adf' in href:
            self.links.append(urljoin('https://open.adas.ac.uk/', href))


def discover_adas(cache, species):
    url = 'https://open.adas.ac.uk/element/' + species.element.lower()
    raw = cache.get('adas_index_' + species.element + '.html', dict(url=url))
    parser = Links()
    parser.feed(raw.decode('utf-8'))
    # Filter by explicit target-charge suffix; an unlabelled file is not
    # automatically assigned to one charge state.
    pattern = re.compile(r'(?:%5D%5B|#|_)' + species.element.lower() + str(species.charge) + r'(?:ls|ic)?(?:_|\.)', re.I)
    return sorted({link for link in parser.links if pattern.search(link)})


def coverage(levels, rates, lines, catalogue):
    result = []
    for kind in ('levels', 'radiative_lines', 'maxwellian_upsilon', 'electron_excitation', 'electron_deexcitation',
                 'electron_ionization', 'electron_elastic', 'electron_momentum_transfer', 'photoionization',
                 'radiative_recombination', 'differential_scattering', 'charge_exchange', 'penning'):
        count = len(levels) if kind == 'levels' else len(lines) if kind == 'radiative_lines' else len({r['id'] for r in rates}) if kind == 'maxwellian_upsilon' else sum(r['kind'] == kind for r in catalogue)
        expected = None
        if kind in ('electron_ionization', 'electron_elastic', 'electron_momentum_transfer', 'photoionization', 'radiative_recombination'):
            expected = len(levels)
        elif kind in ('electron_excitation', 'electron_deexcitation'):
            expected = sum((b['energy_cm1'] >= a['energy_cm1'] if kind == 'electron_excitation' else b['energy_cm1'] < a['energy_cm1']) for a in levels for b in levels if a['id'] != b['id'])
        status = 'available' if count else 'missing'
        if count and expected is not None and count < expected:
            status = 'partial'
        note = 'Imported data only; not a kinetic-closure certificate.'
        if kind == 'maxwellian_upsilon':
            note = 'NOT sigma(E); cannot be used with an arbitrary EEDF.'
        if kind == 'electron_elastic':
            note = 'Integral elastic XS; not automatically a momentum-transfer or differential XS.'
        if kind in ('radiative_lines',):
            note = 'Positive A from the selected ADF04; omitted weak/metastable decays need separate data.'
        result.append(dict(process=kind, status=status, channels=count, reference_count=expected if expected is not None else '', note=note))
    return result


def native_tables(cache,entry,registry,item,nmax):
    """Canonical source adapter for any element with actual supplied data.

    Explicit URLs/hashes in atom.json; no downloaded Python is executed.
    Electronic energies must be relative to this ion's ground state.
    """
    tables={}
    for name,key in entry['native_csv'].items():
        if name not in ('levels','lines','electron_cross_sections','maxwellian_upsilon'):raise ValueError('Неизвестная native-таблица '+name)
        raw=cache.get(key+'.csv',registry['sources'][key])
        tables[name]=list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig')),delimiter=';'))
    levels=tables.get('levels',[])
    if not levels:raise ValueError('Нужны реальные уровни для '+item.label)
    for row in levels:
        if set(LEVEL_FIELDS)-set(row):raise ValueError('Неполные поля levels.csv')
        if row['species']!=item.label or row['energy_reference']!='ground_of_this_ion':raise ValueError('Неверные вид/шкала энергии native CSV')
        for key in ('n','multiplicity','L','g','source_index'):row[key]=int(row[key])
        for key in ('energy_cm1','mass_u'):row[key]=float(row[key])
        if row['g']<=0 or row['energy_cm1']<0 or row['mass_u']<=0 or not np.all(np.isfinite([row['g'],row['energy_cm1'],row['mass_u']])):raise ValueError('Неверный уровень '+row['id'])
        if abs(row['mass_u']/entry['mass_u']-1)>1e-9:raise ValueError('Разные массы в native-данных')
    if min(r['energy_cm1'] for r in levels)!=0:raise ValueError('Основной native уровень должен иметь E=0')
    all_ids={r['id'] for r in levels}
    if len(all_ids)!=len(levels):raise ValueError('Повтор исходного native-уровня')
    for row in tables.get('lines',[]):
        if row.get('lower') not in all_ids or row.get('upper') not in all_ids:raise ValueError('Линия ссылается на неизвестный исходный уровень')
    for row in tables.get('electron_cross_sections',[]):
        if row.get('initial') not in all_ids or (row.get('final') not in all_ids and row.get('kind')!='electron_ionization'):raise ValueError('sigma ссылается на неизвестный исходный уровень')
    levels=[r for r in levels if r['n']<=nmax];ids={r['id'] for r in levels}
    if len(ids)!=len(levels):raise ValueError('Повтор native-уровня')
    lines=[r for r in tables.get('lines',[]) if r['lower'] in ids and r['upper'] in ids]
    for row in lines:
        if set(LINE_FIELDS)-set(row) or not np.all(np.isfinite([float(row['A_s']),float(row['wavelength_nm'])])) or float(row['A_s'])<=0 or float(row['wavelength_nm'])<=0:raise ValueError('Неверные native A/длины волн')
    rates=[r for r in tables.get('maxwellian_upsilon',[]) if r['initial'] in ids and r['final'] in ids]
    for row in rates:
        if row.get('units')!='dimensionless_Maxwellian_effective_collision_strength' or float(row['T_K'])<=0 or float(row['value'])<0 or not np.all(np.isfinite([float(row['T_K']),float(row['value'])])):raise ValueError('Неверная native upsilon(T); это не sigma(E)')
    xs=[r for r in tables.get('electron_cross_sections',[]) if r['initial'] in ids and (r['final'] in ids or r['kind']=='electron_ionization')]
    groups={}
    for row in xs:
        if set(XS_FIELDS)-set(row) or row['kind'] not in ('electron_excitation','electron_deexcitation','electron_ionization','electron_elastic','electron_momentum_transfer'):raise ValueError('Неизвестный native sigma канал')
        groups.setdefault(row['id'],[]).append(row)
    catalogue=[]
    by={r['id']:r for r in levels}
    for ident,rows in groups.items():
        tab=np.array([(float(r['energy_eV']),float(r['sigma_cm2'])) for r in rows])
        if len(tab)<2 or not np.all(np.isfinite(tab)) or np.any(tab<0) or np.any(np.diff(tab[:,0])<=0):raise ValueError('Неверная native sigma(E): '+ident)
        row=rows[0];a,b=row['initial'],row['final'];gap=(by[b]['energy_cm1']-by[a]['energy_cm1'])*HC_EV_CM if b in by else entry['ionization_eV']-by[a]['energy_cm1']*HC_EV_CM
        catalogue.append(dict(id=ident,kind=row['kind'],initial=a,final=b,energy_min_eV=tab[0,0],energy_max_eV=tab[-1,0],model_gap_eV=gap,source_first_energy_eV=tab[0,0],source=row['source'],source_sha256='recorded_in_manifest'))
    threshold=entry.get('ionization_eV')
    if not threshold or threshold<=0:raise ValueError('Нужен положительный порог ионизации native-стадии')
    return levels,rates,lines,xs,catalogue,dict(parents=[dict(parent=1,threshold_cm1=threshold/HC_EV_CM,term='ground')],unsupported_transition_codes={},native_csv=True)


def prepare_species(config, root=ROOT, progress=None):
    root = Path(root)
    progress = progress or Progress()
    options = config.get('atomic_data', {})
    allowed = {'species', 'catalogue_file', 'cache_directory', 'prepared_directory', 'offline', 'nmax', 'discover','element_directory'}
    if set(options) - allowed:
        raise ValueError('Неизвестные atomic_data: ' + str(set(options) - allowed))
    species = selected_species(config)
    nmax = options.get('nmax', 3)
    if isinstance(nmax, bool) or not isinstance(nmax, int) or not 1 <= nmax <= 9:
        raise ValueError('atomic_data.nmax: целое 1..9')
    for flag in ('offline', 'discover'):
        if flag in options and not isinstance(options[flag], bool):
            raise ValueError('atomic_data.' + flag + ' должен быть bool')
    registry = load_catalogue(root,options)
    registry_digest=hashlib.sha256(json.dumps(registry,sort_keys=True).encode()).hexdigest()
    if registry.get('schema') != 'AtomicSources-1':
        raise ValueError('Нужна схема AtomicSources-1')
    output = root / options.get('prepared_directory', 'atomic_data')
    output.mkdir(parents=True, exist_ok=True)
    cache = SourceCache(root / options.get('cache_directory', '.atomic_cache'), options.get('offline', False), progress)
    summaries, all_coverage = [], []
    progress.start('Атомные данные', 'Поиск и проверка источников', ', '.join(s.label for s in species))
    for i, item in enumerate(species):
        progress.module('Подготовка ' + item.label, force=True)
        destination = pack_directory(output,item,nmax,options)
        destination.parent.mkdir(parents=True,exist_ok=True)
        entry = registry['species'].get(item.key)
        if destination.exists():
            manifest=json.loads((destination/'manifest.json').read_text(encoding='utf-8'))
            fingerprint=manifest.get('catalogue_selection_sha256')
            changed=(fingerprint!=catalogue_selection_hash(registry,item.key)) if fingerprint else (manifest.get('catalogue_sha256')!=registry_digest)
            if changed:
                raise ValueError(f'{destination}: изменился каталог источников; выберите другую prepared_directory')
            for name,digest in manifest['files_sha256'].items():
                if not (destination/name).is_file() or hashlib.sha256((destination/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError(f'{destination/name}: ручные правки сохранены; для новой автоподготовки выберите другую prepared_directory')
            def read_rows(name):
                with (destination/name).open(encoding='utf-8-sig',newline='') as stream:return list(csv.DictReader(stream,delimiter=';'))
            cov=read_rows('process_coverage.csv');levels=read_rows('levels.csv');lines=read_rows('lines.csv');catalogue=read_rows('cross_section_catalogue.csv')
            summaries.append(dict(species=item.label,element=item.element,charge=item.charge,levels=len(levels),lines=len(lines),cross_section_channels=len(catalogue),
                missing_processes=sum(r['status']=='missing' for r in cov),status='incomplete' if manifest.get('gaps') or any(r['status']!='available' for r in cov) else 'imported',
                runtime_adapter=manifest.get('runtime_adapter') or 'not_connected',directory=str(destination),download_errors=manifest.get('download_errors',0)))
            all_coverage.extend(dict(species=item.label,**r) for r in cov)
            cache.receipts.update(manifest['sources'])
            progress.update(i+1,len(species),item.label+' (локальный набор)',force=True)
            continue
        levels, rates, lines, xs, catalogue, gaps, candidates = [], [], [], [], [], [], []
        manifest = dict(schema='AtomicPack-1', species=item.label, element=item.element, charge=item.charge,
                        nmax=nmax, runtime_adapter=entry.get('runtime_adapter') if entry else None,
                        status='data_preparation_only', energy_reference='ground_of_this_ion', sources={}, gaps=[])
        try:
            if entry and entry.get('native_csv'):
                levels,rates,lines,xs,catalogue,adf_meta=native_tables(cache,entry,registry,item,nmax)
                manifest['adf04']=adf_meta
            elif entry and entry.get('adf04'):
                key = entry['adf04']
                info = registry['sources'][key]
                raw = cache.get(key + '.dat', info)
                levels, rates, lines, adf_meta = parse_adf04_ls3(raw, item, entry['mass_u'], info['url'], nmax)
                manifest['adf04'] = adf_meta
                provider = entry.get('electron_cross_sections')
                if provider in ('ccc_li_form','ccc_alkali_form'):
                    xs, catalogue, more_gaps = alkali_tables(cache, levels,item.element)
                    gaps.extend(more_gaps)
                elif provider:
                    info = registry['sources'][provider]
                    raw = cache.get(provider + '.zip', info)
                    xs, catalogue = ccc_zip_tables(raw, levels, entry.get('ccc_labels'), info['url'])
            else:
                gaps.append('Нет проверенного импортёра/набора для ' + item.label)
            if options.get('discover', True):
                candidates = discover_adas(cache, item)
        except (OSError, ValueError, KeyError, StopIteration, zipfile.BadZipFile) as exc:
            gaps.append(type(exc).__name__ + ': ' + str(exc))
            manifest['download_error'] = True
        cov = coverage(levels, rates, lines, catalogue)
        manifest['sources'] = dict(cache.receipts)
        manifest['gaps'] = gaps
        manifest['download_errors'] = int(manifest.get('download_error',False))+sum('download error' in gap for gap in gaps)
        manifest['catalogue_sha256'] = registry_digest
        manifest['catalogue_selection_sha256'] = catalogue_selection_hash(registry,item.key)
        manifest['kinetic_ready'] = False
        manifest['kinetic_ready_note'] = 'MixtureSlab consumes these actual levels, lines and energy cross sections; kinetic closure still requires explicit coverage checks.'
        if destination.exists():
            # Never replace edited tables, nor silently combine them with a new source.
            raise ValueError(f'{destination} уже существует. Используйте другую prepared_directory, чтобы сохранить ручные правки.')
        with tempfile.TemporaryDirectory(prefix='atomic_pack_', dir=output) as temporary:
            stage = Path(temporary)
            for name, rows, fields in [('levels.csv', levels, LEVEL_FIELDS), ('lines.csv', lines, LINE_FIELDS),
                                       ('maxwellian_upsilon.csv', rates, RATE_FIELDS), ('electron_cross_sections.csv', xs, XS_FIELDS),
                                       ('cross_section_catalogue.csv', catalogue, CATALOG_FIELDS)]:
                write_csv(stage/name, rows, fields)
            write_csv(stage/'process_coverage.csv', cov, ['process', 'status', 'channels', 'reference_count', 'note'])
            write_json(stage/'candidates.json', dict(species=item.label, candidates=candidates, policy='metadata only; no automatic scientific source selection'))
            manifest['files_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in stage.iterdir() if p.is_file()}
            write_json(stage/'manifest.json', manifest)
            # Atomic rename exposes only a complete pack, with no partially written CSV.
            stage.rename(destination)
        summaries.append(dict(species=item.label, element=item.element, charge=item.charge, levels=len(levels),
                              lines=len(lines), cross_section_channels=len(catalogue), missing_processes=sum(r['status'] == 'missing' for r in cov),
                              status='incomplete' if gaps or any(r['status'] != 'available' for r in cov) else 'imported',
                              runtime_adapter=manifest['runtime_adapter'] or 'not_connected', directory=str(destination),download_errors=manifest['download_errors']))
        all_coverage.extend(dict(species=item.label, **r) for r in cov)
        progress.update(i+1, len(species), item.label, force=True)
    write_csv(output/'species_summary.csv', summaries, list(summaries[0]))
    write_csv(output/'process_coverage.csv', all_coverage, ['species', 'process', 'status', 'channels', 'reference_count', 'note'])
    write_json(output/'preparation_report.json', dict(schema='AtomicPreparation-1', species=summaries, sources=cache.receipts,
        all_requested_species_imported=all(r['levels'] > 0 for r in summaries), multi_element_solver_connected=True))
    progress.done('Таблицы и список недостающих процессов сохранены')
    return summaries
