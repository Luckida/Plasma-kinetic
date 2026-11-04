"""Offline particle/level choices read from the same installed data bank."""
from pathlib import Path
import csv,json


def particle_choices(root,config):
    root=Path(root);folder=root/config.get('atomic_data',{}).get('prepared_directory','../data/elements');choices=[];seen=set()
    for file in sorted(folder.glob('*/atom.json')):
        doc=json.loads(file.read_text(encoding='utf-8-sig'));element=doc.get('element',file.parent.name)
        if isinstance(element,dict):element=element.get('symbol',file.parent.name)
        choices.append(str(element))
    for file in sorted(folder.glob('*/*/levels.csv')):
        with file.open(encoding='utf-8-sig',newline='') as stream:
            for row in csv.DictReader(stream,delimiter=';'):
                ident=row['id']
                if ident not in seen:choices.append(ident+' · '+row.get('label',ident));seen.add(ident)
    components=root/config.get('chemistry',{}).get('components_directory','../data/components')
    for file in sorted(components.glob('*.json')):
        for row in json.loads(file.read_text(encoding='utf-8-sig')).get('components',[]):choices.append(row['id'])
    # Short stage aliases are convenient, and still sourced from descriptors.
    for file in sorted(folder.glob('*/atom.json')):
        doc=json.loads(file.read_text(encoding='utf-8-sig'))
        for row in doc.get('ground_stages',[]):
            charge=int(row.get('charge',0));choices.append(file.parent.name+'+'*charge)
    return list(dict.fromkeys(choices))


def input_name(value):return value.split(' · ',1)[0].strip()
