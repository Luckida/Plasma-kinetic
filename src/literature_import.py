"""Reproduce numeric Li--He data extraction from the authors' printed curves.

Optional preparation dependency: PyMuPDF. Runtime uses the resulting CSV only.
The PDF is not redistributed. Thin model curves are explicitly excluded.
Input: the 2015 paper reproduced on pages 121,123 of Augustovicova's thesis.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from atomic_database import write_csv,write_json

SOURCE='Belyaev, Rodionov, Augustovicova, Soldan & Kraemer 2015; doi:10.1093/mnras/stv391; vector extraction of Figs 5-7 (not CDS numeric tables)'


def points(drawing):
    out=[]
    for segment in drawing['items']:
        if segment[0]!='l':raise ValueError('Expected the original piecewise-linear PDF path')
        out.extend([(float(segment[1].x),float(segment[1].y)),(float(segment[2].x),float(segment[2].y))])
    xy=np.unique(np.asarray(out),axis=0)
    return xy[np.argsort(xy[:,0])]


def extract(pdf,destination):
    import fitz
    destination=Path(destination);old={};old_cross=[];old_rates=[]
    if (destination/'reactions.json').exists():old=json.loads((destination/'reactions.json').read_text(encoding='utf-8-sig'))
    for filename,target in [('cross_sections.csv',old_cross),('rate_coefficients.csv',old_rates)]:
        path=destination/filename
        if path.exists():
            with path.open(encoding='utf-8-sig',newline='') as stream:target.extend(r for r in csv.DictReader(stream,delimiter=';') if not r['id'].startswith('BR15:'))
    doc=fitz.open(pdf);page=doc[120];draw=page.get_drawings();groups=defaultdict(list)
    # Drawing ids and calibrated axes belong to the actual primary PDF,
    # not an image fit. 2..7 are the asymptotes of Table 1 of the paper.
    specs=[(3,2,4,0),(4,2,6,0),(5,2,7,0),(109,4,6,1),(110,4,7,1),
           (261,6,7,2),(262,6,7,2),(467,2,3,3),(468,2,5,3),(469,2,6,3),
           (593,3,5,4),(594,3,6,4),(766,5,6,5),(767,5,6,5)]
    axes=[(194.184,334.578,105.293,204.076,0.,2.),(334.578,474.887,105.293,204.076,0.,2.),
          (194.184,334.578,204.076,302.887,0.,2.),(194.184,334.580,398.872,497.655,-1.,2.),
          (334.580,474.892,398.872,497.655,-1.,2.),(194.184,334.580,497.655,596.466,-1.,2.)]
    audit=[]
    for idx,a,b,panel in specs:
        item=draw[idx]
        if not 1.2<item['width']<2.2:raise ValueError('PDF layout changed: source curve identification failed')
        xy=points(item);xmin,xmax,ymin,ymax,lo,hi=axes[panel]
        ee=10**(lo+(xy[:,0]-xmin)/(xmax-xmin)*(hi-lo))
        # Printed top/bottom axes are 10^3 and 10^-7 in units 10^-16 cm2.
        sigma=10**(3.-(xy[:,1]-ymin)/(ymax-ymin)*10.)*1e-16
        groups[(a,b)].append(np.column_stack((ee,sigma)))
        audit.append(dict(drawing=idx,initial_channel=a,final_channel=b,panel=panel,width=float(item['width']),dashes=item['dashes'],points=len(xy),color=list(item['color'])))
    energy={2:19.19567,3:19.81961,4:20.61577,5:20.96410,6:21.04352,7:21.21802}
    weights={2:4,3:3,4:1,5:9,6:12,7:3}
    channels={2:{'Li(2s)':1,'He+':1},3:{'Li+':1,'He(2^3S)':1},4:{'Li+':1,'He(2^1S)':1},
              5:{'Li+':1,'He(2^3P)':1},6:{'Li(2p)':1,'He+':1},7:{'Li+':1,'He(2^1P)':1}}
    selectors={'Li(2s)':dict(element='Li',charge=0,n=2,multiplicity=2,L=0),'Li(2p)':dict(element='Li',charge=0,n=2,multiplicity=2,L=1)}
    for mult,orb in [(3,0),(1,0),(3,1),(1,1)]:selectors[f'He(2^{mult}{"SP"[orb]})']=dict(element='He',charge=0,n=2,multiplicity=mult,L=orb)
    rows=list(old_cross);reactions=[r for r in old.get('reactions',[]) if not r['id'].startswith('BR15:')]
    for (a,b),tables in sorted(groups.items()):
        gap=energy[b]-energy[a]
        nodes=np.unique(np.r_[gap,np.concatenate([tab[:,0] for tab in tables])])
        nodes=nodes[(nodes>=gap)&(nodes<=100.0002)];nodes[-1]=min(nodes[-1],100.)
        nodes=np.unique(nodes)
        values=np.zeros_like(nodes)
        for table in tables:
            values+=np.interp(nodes,table[:,0],table[:,1],left=0.,right=table[-1,1])
        values[0]=0.
        # Statistical factors are already included in the primary curves.
        # Reverse cross sections follow microscopic detailed balance, including
        # the full atomic LS degeneracies; reduced masses are identical here.
        for initial,final,ee,ss in [(a,b,nodes,values),(b,a,nodes[1:]-gap,weights[a]/weights[b]*nodes[1:]/(nodes[1:]-gap)*values[1:])]:
            ident=f'BR15:{initial}->{final}'
            for e,s in zip(ee,ss):rows.append(dict(id=ident,energy_eV=float(e),sigma_cm2=float(s),energy_frame='centre_of_mass',units='cm2',source=SOURCE))
            kind='charge_exchange' if ('He+' in channels[initial])!=('He+' in channels[final]) else 'heavy_quenching'
            reactions.append(dict(id=ident,kind=kind,reactants=channels[initial],products=channels[final],cross_section_id=ident,
                coefficient_convention='unordered_cross_section',source=SOURCE,
                model_note='Printed-curve digitization. Values below the printed floor (~1e-23 cm2) are unavailable and set to zero; threshold and inverse detailed balance are explicit. Above 100 eV no data. No extra statistical averaging.'))
    destination.mkdir(parents=True,exist_ok=True)
    write_csv(destination/'cross_sections.csv',rows,['id','energy_eV','sigma_cm2','energy_frame','units','source'])
    # Independent thermal rate curves are useful both for checking the
    # digitized cross sections and for an explicitly chosen moment closure.
    rate_specs=[(6,2,4,0),(7,2,6,0),(8,2,7,0),(9,2,3,0),(10,2,5,0),(11,2,6,0),
                (152,3,2,1),(153,3,5,1),(154,3,6,1),(295,4,2,2),(296,4,6,2),(297,4,7,2),
                (491,5,2,3),(492,5,3,3),(493,5,6,3),(494,5,6,3),(597,6,2,4),(598,6,4,4),
                (599,6,7,4),(600,6,2,4),(601,6,3,4),(602,6,5,4),(603,6,5,4),(604,6,7,4),
                (703,7,2,5),(704,7,4,5),(705,7,6,5),(706,7,6,5)]
    rate_axes=[(194.104645,288.232697,104.766594,204.178223),(288.232605,382.332123,104.766533,204.178162),
               (382.332031,476.461029,104.766533,204.178162),(194.104294,288.232330,204.178223,303.590759),
               (288.232208,382.331726,204.178223,303.590759),(382.331726,476.460632,204.178223,303.590759)]
    rate_groups=defaultdict(list);rate_draw=doc[122].get_drawings()
    for idx,a,b,panel in rate_specs:
        xy=points(rate_draw[idx]);xmin,xmax,ymin,ymax=rate_axes[panel]
        temp=10**(2+(xy[:,0]-xmin)/(xmax-xmin)*4);value=10**(-7-(xy[:,1]-ymin)/(ymax-ymin)*10)
        good=(value>1.02e-17)&(temp<=1.00001e6)
        table=np.column_stack((temp[good],value[good]));table[-1,0]=min(table[-1,0],1e6)
        rate_groups[(a,b)].append(table)
    fit_rows=list(old_rates)
    for (a,b),tables in sorted(rate_groups.items()):
        lo=min(t[0,0] for t in tables);hi=min(t[-1,0] for t in tables)
        temperatures=np.unique(np.r_[lo,np.geomspace(lo,hi,81),hi])
        value=sum(np.exp(np.interp(np.log(temperatures),np.log(t[:,0]),np.log(t[:,1]),left=-np.inf)) for t in tables)
        ident=f'BR15:{a}->{b}'
        for t,k in zip(temperatures,value):fit_rows.append(dict(id=ident,model='tabulated_loglog',alpha_cm3_s=float(k),beta=0.,gamma_K=0.,
            Tmin_K=float(lo),Tmax_K=float(hi),reference_T_K=float(t),source=SOURCE+'; Fig.7',method='digitized_quantum_curves',accuracy='plot/interpolation; omitted partials below 1e-17 cm3/s'))
        for reaction in reactions:
            if reaction['id']==ident:reaction['rate_id']=ident
    write_csv(destination/'rate_coefficients.csv',fit_rows,['id','model','alpha_cm3_s','beta','gamma_K','Tmin_K','Tmax_K','reference_T_K','source','method','accuracy'])
    write_json(destination/'reactions.json',dict(schema='ParticleReactions-1',reviewed=True,elements=['He','Li'],selectors=dict(old.get('selectors',{}),**selectors),reactions=reactions))
    write_json(destination/'provenance.json',dict(source=SOURCE,source_pdf_sha256=hashlib.sha256(Path(pdf).read_bytes()).hexdigest(),
        primary_copy='https://www.mff.cuni.cz/cs/vedecka-rada/habilitacni-prace-dr-augustovicova.pdf',pages_1based=[121,123],
        method='vector paths, not raster digitization; axes calibrated against PDF tick positions',sigma_axes=axes,curves=audit,
        missing='No CDS numbers were received. All vectors below the plotted floor are missing. Accuracy is limited by original plotting and interpolation, not the underlying quantum solver.'))
    print(f'{len(reactions)} reactions, {len(rows)} sigma(E) points saved',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pdf',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();extract(args.pdf,args.out)
