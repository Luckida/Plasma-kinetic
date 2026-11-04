"""Explicit network/element import: no code downloaded, no guessed level mapping.

Commands: prepare, network (UMIST Rate22), register (atomic descriptor),
lxcat (export with an explicit mapping), validate (build the complete model).
"""
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
import tempfile
from atomic_core import SYMBOLS,Species,selected_species
from atomic_database import SourceCache,prepare_species,load_catalogue,write_csv,write_json,ROOT


UMIST_URL='https://umistdatabase.uk/files/rate22_final.rates'


def elements_in(name):
    if name in ('e','e-','PHOTON','CRP','CRPHOT',''):return set()
    cleaned=re.sub(r'^(?:c-|l-|o-|p-|m-)+','',name).rstrip('+-0123456789')
    found=set(re.findall(r'[A-Z][a-z]?',cleaned))
    if not found or not found<=set(SYMBOLS):raise ValueError('Неизвестная формула '+name)
    return found


def parse_umist(raw,selected_elements):
    """Keep every temperature interval, method, uncertainty and original ID.

    Reader honours quoted URLs/notes containing ':'; no split(':') shortcut.
    Photochemical/cosmic-ray/thermal/dipole laws are NOT conflated.
    """
    groups={};rejected=[];total=0;selected=0
    reader=csv.reader(io.StringIO(raw.decode('utf-8-sig')),delimiter=':',quotechar='"')
    for fields in reader:
        if not fields or not fields[0].strip() or fields[0].startswith('#'):continue
        total+=1
        try:
            number=int(fields[0]);kind=fields[1];ntr=int(fields[8])
            if ntr<1 or len(fields)<9+9*ntr:raise ValueError('Повреждены температурные интервалы')
            all_elements=set().union(*(elements_in(n) for n in fields[2:8]))
            if not all_elements or not all_elements<=set(selected_elements):continue
            selected+=1;key='--'.join(sorted(all_elements));group=groups.setdefault(key,dict(reactions=[],rates=[],unsupported=[]))
            left=Counter('e' if n=='e-' else n for n in fields[2:4] if n)
            right=Counter('e' if n=='e-' else n for n in fields[4:8] if n)
            ident='UMIST22:'+str(number)
            source='UMIST Rate22; '+fields[16]+'; '+fields[17]
            reason=''
            if kind in ('CP','CR','PH','RA','RR','REA'):reason='Requires radiation/cosmic-ray operator or overlaps level-resolved atomic recombination; not a binary chemical rate.'
            elif sum(left.values())!=2 or sum(right.values()) not in (2,3):reason='Particle operator currently supports 2 -> 2/3 only.'
            elif right.get('e',0)>1:reason='Use level-resolved electron-impact operator.'
            elif kind not in ('AD','CD','CE','DR','IN','MN','NN'):reason='Rate law not reviewed for this importer: '+kind
            row=dict(id=ident,reactants=dict(left),products=dict(right),kind={'CE':'charge_exchange','DR':'dissociative_recombination','MN':'mutual_neutralization'}.get(kind,'heavy_chemistry'),
                     rate_id=ident,coefficient_convention='event_rate',source=source,original_type=kind,requires_validation=True)
            if reason:group['unsupported'].append(dict(id=ident,reason=reason,source=source))
            else:
                group['reactions'].append(row)
                for i in range(ntr):
                    part=fields[9+9*i:18+9*i]
                    values=list(map(float,part[:5]))
                    if values[0]<0 or not 0<values[3]<=values[4]:raise ValueError('Неверные коэффициенты/диапазон')
                    group['rates'].append(dict(id=ident,model='arrhenius',alpha_cm3_s=values[0],beta=values[1],gamma_K=values[2],Tmin_K=values[3],Tmax_K=values[4],
                        reference_T_K=300.,source='UMIST Rate22; '+part[7],method=part[5],accuracy=part[6],notes=part[8]))
        except (ValueError,IndexError) as exc:rejected.append(dict(line=reader.line_num,reason=str(exc)))
    return groups,dict(source_records=total,selected_records=selected,rejected=rejected)


def save_network(raw,elements,out,url='',replace=False):
    groups,report=parse_umist(raw,elements);out=Path(out)
    if out.exists() and any(out.iterdir()) and not replace:raise ValueError('Папка уже непуста; выберите новую или --replace')
    out.mkdir(parents=True,exist_ok=True)
    for key,group in groups.items():
        folder=out/key;folder.mkdir(exist_ok=True)
        write_json(folder/'reactions.json',dict(schema='ParticleReactions-1',reviewed=False,elements=key.split('--'),selectors={},reactions=group['reactions'],
            note='Imported candidates. Supply ground energies/levels, check collisions and duplicate channels before setting reviewed=true.'))
        rows=group['rates'];fields=['id','model','alpha_cm3_s','beta','gamma_K','Tmin_K','Tmax_K','reference_T_K','source','method','accuracy','notes']
        write_csv(folder/'rate_coefficients.csv',rows,fields)
        write_csv(folder/'unsupported.csv',group['unsupported'],['id','reason','source'])
        write_csv(folder/'cross_sections.csv',[],['id','energy_eV','sigma_cm2','energy_frame','units','source'])
    report.update(url=url,sha256=hashlib.sha256(raw).hexdigest(),retrieved_utc=datetime.now(timezone.utc).isoformat(),selected_elements=elements,groups=len(groups),
                  candidate_reactions=sum(len(g['reactions']) for g in groups.values()),status='candidates; not auto-connected',
                  source_kind='rate_coefficients_only; contains no energy-resolved cross sections')
    write_json(out/'import_report.json',report)
    return report


def lxcat_blocks(raw):
    """Parse common LXCat/BOLSIG text without guessing which level is meant."""
    lines=raw.decode('utf-8-sig').splitlines();blocks=[];i=0
    kinds={'ELASTIC','EFFECTIVE','EXCITATION','IONIZATION','ATTACHMENT'}
    while i<len(lines):
        kind=lines[i].strip();i+=1
        if kind not in kinds:continue
        header=[]
        while i<len(lines) and not re.fullmatch(r'-{5,}',lines[i].strip()):header.append(lines[i]);i+=1
        i+=1;points=[]
        while i<len(lines) and not re.fullmatch(r'-{5,}',lines[i].strip()):
            text=lines[i].strip();i+=1
            if not text or text.startswith('#'):continue
            columns=text.replace('D','E').split()
            if len(columns)!=2:raise ValueError('Ожидались два столбца в LXCat')
            points.append(tuple(map(float,columns)))
        if len(points)<2:raise ValueError('Недостаточно точек в LXCat-блоке')
        blocks.append(dict(block_index=len(blocks)+1,type=kind,header='\n'.join(header),points=points))
    if not blocks:raise ValueError('Не найдены LXCat/BOLSIG-блоки; формат не угадан автоматически')
    return blocks


def convert_lxcat(raw,mapping,source):
    if mapping.get('projectile')!='electron':raise ValueError('Явно укажите projectile=electron в соответствии LXCat')
    blocks=lxcat_blocks(raw);rows=[];used=set()
    allowed={'ELASTIC':{'electron_elastic','electron_momentum_transfer'},'EXCITATION':{'electron_excitation','electron_deexcitation'},'IONIZATION':{'electron_ionization'}}
    for item in mapping['channels']:
        number=item['block_index']
        if number in used or not isinstance(number,int) or not 1<=number<=len(blocks):raise ValueError('Повтор/неизвестный LXCat-блок')
        used.add(number);block=blocks[number-1]
        if item['kind'] not in allowed.get(block['type'],set()):raise ValueError('Несовместимые типы LXCat/canonical')
        energy_scale=float(mapping['energy_scale_to_eV']);area_scale=float(mapping['area_scale_to_cm2'])
        if min(energy_scale,area_scale)<=0:raise ValueError('Явно задайте положительные преобразования единиц')
        last=-1.
        for energy,sigma in block['points']:
            energy*=energy_scale;sigma*=area_scale
            if not energy>last or sigma<0 or not math.isfinite(sigma):raise ValueError('Неверные точки LXCat')
            last=energy
            rows.append(dict(id=item['id'],kind=item['kind'],initial=item['initial'],final=item['final'],energy_eV=energy,sigma_cm2=sigma,source=source+'; block '+str(number)))
    return rows,[dict(block_index=b['block_index'],type=b['type'],header=b['header']) for b in blocks if b['block_index'] not in used]


def register_element(descriptor,element_directory,replace=False):
    """Validate source tables before making a new element discoverable."""
    document=json.loads(Path(descriptor).read_text(encoding='utf-8-sig'))
    element=document['element']
    if document.get('schema')!='ElementData-1' or element not in SYMBOLS:raise ValueError('Нужен ElementData-1')
    stages=document['ground_stages']
    if [r['charge'] for r in stages]!=list(range(SYMBOLS.index(element)+2)):raise ValueError('Нужны последовательные пороги всех ионных стадий до голого ядра')
    for row in stages[:-1]:
        if row.get('ionization_eV',0)<=0:raise ValueError('Нужны положительные пороги каждой стадии')
    for key,entry in document['species'].items():
        if not key.startswith(element+'_') or entry.get('mass_u',0)<=0:raise ValueError('Несогласованный элемент/масса')
    for info in document['sources'].values():
        if not info['url'].startswith('https://') or not re.fullmatch('[0-9a-f]{64}',info.get('sha256','')):raise ValueError('Для нового источника нужны HTTPS и точный SHA256')
    destination=Path(element_directory)/element/'atom.json'
    if destination.exists() and not replace:raise ValueError('Элемент уже зарегистрирован; --replace нужен для явного обновления')
    with tempfile.TemporaryDirectory(prefix='element_check_') as temporary:
        root=Path(temporary);folder=root/'elements'/element;folder.mkdir(parents=True);write_json(folder/'atom.json',document)
        cfg=dict(atomic_data=dict(species=[Species(element,int(k.rsplit('_',1)[1])).label for k in document['species']],element_directory='elements',nmax=9,discover=False))
        result=prepare_species(cfg,root)
        if any(not r['levels'] or r['download_errors'] for r in result):raise ValueError('Не прошла проверка импортёра нового элемента: '+str(result))
    destination.parent.mkdir(parents=True,exist_ok=True);write_json(destination,document)
    return destination


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    prepare=sub.add_parser('prepare');prepare.add_argument('--config',type=Path,required=True);prepare.add_argument('--offline',action='store_true')
    network=sub.add_parser('network');network.add_argument('--elements',nargs='+',required=True);network.add_argument('--url',default=UMIST_URL);network.add_argument('--file',type=Path);network.add_argument('--sha256');network.add_argument('--out',type=Path,required=True);network.add_argument('--replace',action='store_true')
    register=sub.add_parser('register');register.add_argument('--descriptor',type=Path,required=True);register.add_argument('--element-directory',type=Path,default=ROOT/'data/elements');register.add_argument('--replace',action='store_true')
    lxcat=sub.add_parser('lxcat');lxcat.add_argument('--file',type=Path);lxcat.add_argument('--url');lxcat.add_argument('--sha256');lxcat.add_argument('--mapping',type=Path,required=True);lxcat.add_argument('--out',type=Path,required=True)
    validate=sub.add_parser('validate');validate.add_argument('--config',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.command in ('prepare','validate'):
        cfg=json.loads(args.config.read_text(encoding='utf-8-sig'));root=args.config.resolve().parent
        if args.command=='prepare':
            from mixture_data import MixtureData
            if args.offline:cfg.setdefault('atomic_data',{})['offline']=True
            bank=MixtureData(cfg,root)
            print('Подготовлены:',', '.join(s.label for s in bank.imported),'; состояний:',len(bank.states))
        else:
            from mixture_kinetics import MixtureSlab
            model=MixtureSlab(cfg,root)
            print('Построено:',model.ns,'состояний;',model.nr,'каналов;',len(model.chemistry.reactions),'межкомпонентных реакций')
            print('Это проверка структуры, не доказательство полноты данных или сходимости.')
    elif args.command=='register':print(register_element(args.descriptor,args.element_directory,args.replace))
    elif args.command=='network':
        if not set(args.elements)<=set(SYMBOLS):raise ValueError('Неизвестный элемент')
        raw=args.file.read_bytes() if args.file else SourceCache(args.out.parent/'source_cache').get('rate22.rates',dict(url=args.url,**({'sha256':args.sha256} if args.sha256 else {})))
        if args.sha256 and hashlib.sha256(raw).hexdigest()!=args.sha256:raise ValueError('Неверный SHA256 сети')
        result=save_network(raw,args.elements,args.out,args.url,replace=args.replace)
        print(json.dumps(result,ensure_ascii=False,indent=2))
    else:
        if bool(args.file)==bool(args.url):raise ValueError('Выберите ровно один источник: --file или --url')
        if args.url and not args.sha256:raise ValueError('Для сетевого LXCat-экспорта нужен --sha256')
        raw=args.file.read_bytes() if args.file else SourceCache(args.out.parent/'source_cache').get('lxcat.txt',dict(url=args.url,sha256=args.sha256))
        if args.sha256 and hashlib.sha256(raw).hexdigest()!=args.sha256:raise ValueError('Неверный SHA256 LXCat')
        rows,unmapped=convert_lxcat(raw,json.loads(args.mapping.read_text(encoding='utf-8-sig')),args.url or args.file.name)
        args.out.parent.mkdir(parents=True,exist_ok=True)
        write_csv(args.out,rows,['id','kind','initial','final','energy_eV','sigma_cm2','source'])
        write_json(args.out.with_suffix('.provenance.json'),dict(sha256=hashlib.sha256(raw).hexdigest(),unmapped=unmapped,rows=len(rows),mapping=json.loads(args.mapping.read_text(encoding='utf-8-sig'))))
        print('Записаны реальные sigma(E):',args.out,'; не сопоставлено блоков:',len(unmapped))
    return 0


if __name__=='__main__':
    try:sys.exit(main())
    except (ValueError,OSError,KeyError) as exc:sys.exit('Ошибка импорта: '+str(exc))
