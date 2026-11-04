"""SYNTHETIC UI/network fixture. Never used as physical plasma data."""
import json
from pathlib import Path
import numpy as np
from result_store import FILE_STORE


class FixtureData:
    def __init__(self,all_zero=False,empty=False):
        self.states=[]
        for ident,element,species,q,n,mult,L,energy,caption in [
            ('H0','H','H I',0,1,2,0,0.,'H I 1²S'),('H1','H','H II',1,0,1,0,13.,'H⁺'),
            ('He0','He','He I',0,1,1,0,0.,'He I 1¹S'),('HeS','He','He I',0,2,3,0,19.,'He I 2³S'),
            ('HeP','He','He I',0,2,3,1,20.,'He I 2³P'),('He1','He','He II',1,1,2,0,24.,'He II 1²S'),
            ('Li0','Li','Li I',0,2,2,0,0.,'Li I 2²S'),('Li1','Li','Li II',1,1,1,0,5.,'Li II 1¹S')]:
            self.states.append(dict(id=ident,element=element,species=species,charge=q,n=n,multiplicity=mult,L=L,energy_cm1=energy*1000,label=caption,display_label=caption,nuclei={element:1},molecular=False,mass_u=1.))
        recipes=[('EX','excitation',{'He0':1,'e':1},{'HeS':1,'e':1},9.),('DX','deexcitation',{'HeS':1,'e':1},{'He0':1,'e':1},8.),
            ('PI_H','penning',{'HeS':1,'H0':1},{'He0':1,'H1':1,'e':1},7.),('PI_Li','penning',{'HeS':1,'Li0':1},{'He0':1,'Li1':1,'e':1},0.),
            ('SELF','associative_ionization',{'HeS':2},{'He0':1,'He1':1,'e':1},6.),('A_S','spontaneous',{'HeS':1},{'He0':1,'photon':1},5.),
            ('PUMP','photoexcitation',{'He0':1,'photon':1},{'HeP':1},.3),('CX','charge_exchange',{'H0':1,'H1':1},{'H1':1,'H0':1},4.),
            ('A_P','spontaneous',{'HeP':1},{'HeS':1,'photon':1},.7),('PI_He','penning',{'HeS':1,'He0':1},{'He1':1,'He0':1,'e':1},1.)]
        if empty:recipes=[]
        labels={s['id']:s['label'] for s in self.states};labels.update(e='e⁻',photon='hν');self.reactions=[]
        def side(items):return ' + '.join((str(count)+' ' if count>1 else '')+labels[name] for name,count in items.items())
        for ident,kind,left,right,rate in recipes:self.reactions.append(dict(id=ident,kind=kind,reactants=left,products=right,left=side(left),right=side(right),reaction=side(left)+' → '+side(right),source='SYNTHETIC UI TEST ONLY; not a physical coefficient'))
        matrix=np.zeros((len(self.states),len(recipes)));electrons=np.zeros(len(recipes))
        for k,r in enumerate(self.reactions):
            for j,s in enumerate(self.states):matrix[j,k]=r['products'].get(s['id'],0)-r['reactants'].get(s['id'],0)
            electrons[k]=r['products'].get('e',0)-r['reactants'].get('e',0)
        rates=np.broadcast_to(np.array([r[-1] for r in recipes]),(3,2,len(recipes))).copy()
        if recipes:rates[0,:,0]=30.;rates[-1,:,1]=40.;rates[:,1,6]=50.
        if all_zero:rates[:]=0.
        self.meta=dict(version='2.6.0',states=self.states,reactions=self.reactions,components={name:[j for j,s in enumerate(self.states) if s['species']==name] for name in sorted({s['species'] for s in self.states})},sources=[],config={'processes':{}},quality={'data_issues':[]})
        self.arrays=dict(time_s=np.array([0.,.5,1.]),position_cm=np.array([.25,.75]),reaction_rates=rates,reaction_integrals=rates.copy(),stoichiometry=matrix,electron_stoichiometry=electrons,n_cm3=np.ones((3,2,len(self.states))),electrons_cm3=np.ones((3,2,16,4)))

    def __getitem__(self,key):return self.arrays[key]

    def save(self,folder):np.savez_compressed(Path(folder)/FILE_STORE,metadata_json=np.asarray(json.dumps(self.meta,ensure_ascii=False)),**self.arrays)
