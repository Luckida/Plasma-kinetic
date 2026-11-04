"""Spectrum import shared by the GUI, driving light and passive probe.

Text files need two numeric columns, not prescribed column names. Units remain
explicit: auto uses a recognised header, otherwise nm (reported in metadata).
The numerical core still receives wavelength in vacuum Angstrom and cgs/A.
"""
from dataclasses import dataclass
from pathlib import Path
import csv
import io
import math
import re
import zipfile
import xml.etree.ElementTree as ET

HC_EV_A=6.62607015e-27*2.99792458e10/1.602176634e-12*1e8
EV_ERG=1.602176634e-12
C_A_S=2.99792458e18


@dataclass
class Spectrum:
    wavelength_A: tuple
    values: tuple
    metadata: dict

    def canonical_rows(self):
        return [dict(wavelength_A=x,four_pi_J_lambda_erg_cm2_s_A=y) for x,y in zip(self.wavelength_A,self.values)]


def _number(value):
    value=str(value).strip().strip('"').replace('\u2212','-').replace('\u00a0','').replace(',','.').replace('D','e').replace('d','e')
    return float(value)


def _xlsx_rows(path):
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            strings=[''.join(n.itertext()) for n in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',ns)]
        book=ET.fromstring(z.read('xl/workbook.xml'));sheet=book.find('m:sheets/m:sheet',ns)
        rel=sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        relationships=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        target=next(n.attrib['Target'] for n in relationships if n.attrib.get('Id')==rel)
        name=target.lstrip('/') if target.startswith('/') else 'xl/'+target
        for row in ET.fromstring(z.read(name)).findall('.//m:sheetData/m:row',ns):
            values={}
            for cell in row.findall('m:c',ns):
                letters=re.match('[A-Z]+',cell.attrib['r'])[0];col=0
                for letter in letters:col=26*col+ord(letter)-64
                value=cell.find('m:v',ns);text=value.text if value is not None else ''
                if cell.attrib.get('t')=='s':text=strings[int(text)]
                elif cell.attrib.get('t')=='inlineStr':text=''.join(cell.find('m:is',ns).itertext())
                values[col-1]=text
            yield [values.get(k,'') for k in range(max(values,default=-1)+1)]


def _text_rows(path,delimiter='auto'):
    raw=Path(path).read_bytes()
    if raw.startswith((b'\xff\xfe',b'\xfe\xff')):text=raw.decode('utf-16')
    else:
        try:text=raw.decode('utf-8-sig')
        except UnicodeDecodeError:text=raw.decode('cp1251')
    for line in text.splitlines():
        stripped=line.strip()
        if not stripped:continue
        comment=stripped.startswith(('#','!','//'))
        if comment:stripped=re.sub(r'^(?:#|!|//)\s*','',stripped)
        sep=delimiter
        if sep=='auto':
            words=stripped.split();numeric_words=False
            if len(words)>=2:
                try:_number(words[0]);_number(words[1]);numeric_words=True
                except ValueError:pass
            sep=';' if ';'in stripped else '\t' if '\t'in stripped else ' ' if numeric_words else ',' if ','in stripped else ' '
        if sep in ('space','whitespace',' '):row=re.split(r'\s+',stripped)
        else:row=next(csv.reader([stripped],delimiter=sep))
        yield row,comment


def _header_unit(text):
    s=text.lower().replace('å','angstrom').replace('Å','angstrom')
    if re.search(r'\bev\b|energy_e?v|энерги',s):return 'eV'
    if re.search(r'\bhz\b|frequency|частот|гц',s):return 'Hz'
    if re.search(r'\bum\b|µm|μm|micron|мкм',s):return 'um'
    if re.search(r'\bnm\b|_nm|нм',s):return 'nm'
    if re.search(r'angstrom|wavelength_a\b|lambda_a\b|\[a\]|\(a\)|ангстр',s):return 'A'
    return None


def load_spectrum(path,options=None,purpose='probe'):
    """Return sorted, duplicate-averaged pairs and an explicit import receipt.

    x_column/y_column are zero based. y_kind: relative, energy, photons.
    Physical densities are per unit of the *input* x coordinate; Jacobians are
    applied for nm/A/um/eV/Hz, with cgs or W/m2 energy units. No extrapolation.
    """
    path=Path(path);opts=dict(options or {})
    allowed={'x_unit','y_kind','energy_unit','angular_kind','mu','x_column','y_column','delimiter','negative_policy'}
    if set(opts)-allowed:raise ValueError('Неизвестные настройки спектра: '+', '.join(sorted(set(opts)-allowed)))
    if not path.is_file():raise ValueError('Файл спектра не найден: '+str(path))
    if path.suffix.lower()=='.xls':raise ValueError('Старый XLS: сохраните таблицу как XLSX, TXT или CSV; названия столбцов произвольные.')
    ix=int(opts.get('x_column',0));iy=int(opts.get('y_column',1))
    if min(ix,iy)<0 or ix==iy:raise ValueError('Укажите два разных столбца спектра (в интерфейсе нумерация с 1).')
    source=((r,False) for r in _xlsx_rows(path)) if path.suffix.lower()=='.xlsx' else _text_rows(path,opts.get('delimiter','auto'))
    pairs=[];headers=[];inferred=None;canonical=False;warnings=[]
    for line_no,(row,comment) in enumerate(source,1):
        if len(row)<=max(ix,iy):
            if row:headers.append(' '.join(row))
            continue
        try:x=_number(row[ix])
        except (ValueError,TypeError):
            text=' '.join(row);headers.append(text);inferred=inferred or _header_unit(row[ix]) or _header_unit(text)
            canonical|='four_pi_J_lambda_erg_cm2_s_A' in text
            continue
        if comment:continue
        try:y=_number(row[iy])
        except (ValueError,TypeError):raise ValueError(f'Спектр, строка {line_no}: интенсивность «{row[iy]}» не число. Выберите правильные столбцы/разделитель.') from None
        if not math.isfinite(x+y) or x<=0:raise ValueError(f'Спектр, строка {line_no}: координата должна быть >0, оба числа конечны.')
        pairs.append((x,y))
    if len(pairs)<2:raise ValueError('В спектре нужно минимум две числовые строки: длина волны/энергия и интенсивность. Поддержаны TXT, CSV, TSV, DAT, ASC и XLSX.')
    unit=opts.get('x_unit','auto')
    if unit=='auto':
        unit=inferred or 'nm'
        if not inferred:warnings.append('Единицы X не указаны в заголовке: выбран nm. Их можно изменить в настройках импорта.')
    aliases={'нм':'nm','angstrom':'A','Å':'A','мкм':'um','µm':'um','μm':'um','эв':'eV','Гц':'Hz'};unit=aliases.get(unit,unit)
    if unit not in ('A','nm','um','eV','Hz'):raise ValueError('Единицы спектра: A, nm, um, eV или Hz.')
    ykind=opts.get('y_kind','auto')
    if ykind=='auto':ykind='energy' if purpose=='driving' or canonical else 'relative'
    if ykind not in ('relative','energy','photons'):raise ValueError('Вид интенсивности: relative, energy или photons.')
    if purpose=='driving' and ykind=='relative':raise ValueError('Для накачки нужна абсолютная энергетическая/фотонная интенсивность. Относительные отсчёты используйте как пробный спектр.')
    energy_unit=opts.get('energy_unit','cgs')
    if energy_unit not in ('cgs','W_m2'):raise ValueError('Единицы мощности: cgs или W_m2.')
    angular=opts.get('angular_kind','four_pi_J');mu=float(opts.get('mu',1.))
    if angular not in ('four_pi_J','normal_flux') or not math.isfinite(mu) or not 0<mu<=1:raise ValueError('Для направленного света: 0<mu<=1, normal_flux или four_pi_J.')
    converted=[];negative=0
    for x,y in pairs:
        lam=x*{'A':1.,'nm':10.,'um':1e4}.get(unit,1.) if unit in ('A','nm','um') else HC_EV_A/x if unit=='eV' else C_A_S/x
        if y<0:
            negative+=1
            if opts.get('negative_policy','error')=='clip':y=0.
            elif ykind!='relative':raise ValueError('В спектре есть отрицательная физическая интенсивность. Исправьте фон или явно выберите «Обрезать отрицательные».')
        if ykind!='relative':
            jac=1./{'A':1.,'nm':10.,'um':1e4}[unit] if unit in ('A','nm','um') else HC_EV_A/lam**2 if unit=='eV' else C_A_S/lam**2
            y*=jac
            if ykind=='photons':y*=HC_EV_A/lam*EV_ERG
            elif energy_unit=='W_m2':y*=1e3
            if purpose=='driving' and angular=='normal_flux':y/=mu
        if not math.isfinite(lam+y):raise ValueError('После перевода единиц получены бесконечные значения спектра. Проверьте масштаб и единицы.')
        converted.append((lam,y))
    if negative:warnings.append(f'Отрицательных отсчётов: {negative}; '+('заменены нулём по выбранному правилу.' if opts.get('negative_policy')=='clip' else 'сохранены как относительные отсчёты, не физическая мощность.'))
    if any(converted[k][0]>converted[k+1][0] for k in range(len(converted)-1)):warnings.append('Точки отсортированы по длине волны.')
    grouped={}
    for x,y in converted:grouped.setdefault(x,[]).append(y)
    duplicates=len(converted)-len(grouped)
    if duplicates:warnings.append(f'Повторных координат: {duplicates}; интенсивности усреднены.')
    data=sorted((x,sum(ys)/len(ys)) for x,ys in grouped.items())
    if len(data)<2:raise ValueError('После объединения повторов в спектре осталась одна длина волны.')
    import hashlib
    return Spectrum(tuple(x for x,y in data),tuple(y for x,y in data),dict(file=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),input_x_unit=unit,input_y_kind=ykind,energy_unit=energy_unit,angular_kind=angular,points=len(data),wavelength_min_nm=data[0][0]/10,wavelength_max_nm=data[-1][0]/10,warnings=warnings,headers=headers[:8],purpose=purpose))


def driving_rows(path,options=None):
    return load_spectrum(path,options,purpose='driving').canonical_rows()
