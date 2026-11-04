"""Kinetic wall fluxes at x=0,L, with quenching and ion neutralization.

Only particles reaching a surface react. Reflection is specular; quenching
and neutralization transfer internal energy to the surface. Probabilities are
material inputs, not purported universal experimental constants. A surface
charge/current reservoir is explicitly included in the global balance.
"""
import numpy as np
from atomic_database import write_csv
from particle_kinetics import EV


def advect_absorbing(values,vz,dt,dz):
    """Positive upwind finite-volume step, CFL<=1; outgoing phase columns."""
    if np.max(abs(vz))*dt/dz>1+1e-12:raise ValueError('wall transport: CFL>1')
    flux=np.zeros((len(values)+1,*vz.shape))
    flux[1:-1]=np.maximum(vz,0)[None,:,:]*values[:-1]+np.minimum(vz,0)[None,:,:]*values[1:]
    flux[0]=np.minimum(vz,0)*values[0];flux[-1]=np.maximum(vz,0)*values[-1]
    out=values-dt/dz*(flux[1:]-flux[:-1])
    # Roundoff at exactly CFL=1 must not create artificial particles.
    if np.any(out < -1e-13*np.maximum(values.max(),1.)):raise ValueError('Negative wall transport')
    return np.maximum(out,0.),(-dt*flux[0],dt*flux[-1])


class WallBoundaries:
    def __init__(self,model):
        self.model=model;self.enabled=model.bound.get('kind')=='walls';self.surfaces=[]
        self.fields=('incident_cm2','quenched_cm2','neutralized_cm2','stuck_cm2','returned_cm2')
        if not self.enabled:
            if 'walls' in model.bound:raise ValueError('particle_boundary.walls требует kind=walls')
            self.reset();return
        surfaces=model.bound.get('walls',{})
        if set(surfaces)-{'left','right'}:raise ValueError('walls: left и right')
        for side in ('left','right'):
            row=surfaces.get(side,{})
            if set(row)-{'quenching_probability','neutralization_probability','sticking_probability','electron_absorption_probability','overrides'}:raise ValueError('Неизвестные параметры стенки '+side)
            base=dict(quenching_probability=1.,neutralization_probability=1.,sticking_probability=0.,electron_absorption_probability=1.)
            base.update({k:v for k,v in row.items() if k!='overrides'})
            self.check(base);overrides={}
            for name,values in row.get('overrides',{}).items():
                if set(values)-set(base):raise ValueError('Неизвестные параметры wall override')
                opts=dict(base,**values);self.check(opts);overrides[model.data.resolve(name)]=opts
            self.surfaces.append((base,overrides))
        self.ground={}
        for j,state in enumerate(model.states):
            candidates=[k for k,s in enumerate(model.states) if s['element']==state['element'] and s['charge']==0 and s.get('is_ground',False)]
            self.ground[j]=min(candidates,key=lambda k:model.internal[k]) if candidates else j
        self.reset()

    @staticmethod
    def check(row):
        for key,value in row.items():
            if isinstance(value,bool) or not isinstance(value,(float,int)) or not np.isfinite(value) or not 0<=value<=1:raise ValueError('Стенка: '+key+' должно лежать в [0,1]')

    def reset(self):
        self.counts=np.zeros((2,self.model.ns+1,len(self.fields)));self.heat=np.zeros(2);self.charge=np.zeros(2)

    def checkpoint(self):return self.counts.copy(),self.heat.copy(),self.charge.copy()
    def restore(self,state):self.counts,self.heat,self.charge=state

    def apply(self,dt):
        m=self.model;speed=max(float(abs(m.ge.vz).max()),max(float(abs(g.vz).max()) for g in m.gh.values()))
        count=max(1,int(np.ceil(dt*speed/(.8*m.dz))));step=dt/count
        before=m.P.sum(axis=(2,3))
        for _ in range(count):
            nextF,outF=advect_absorbing(m.F,m.ge.vz,step,m.dz)
            nextP=np.empty_like(m.P);outP=[]
            for j,s in enumerate(m.states):
                nextP[:,j],out=advect_absorbing(m.P[:,j],m.gh[s['element']].vz,step,m.dz);outP.append(out)
            for side,z in [(0,0),(1,m.cells-1)]:
                base,overrides=self.surfaces[side]
                outgoing=outF[side];absorbed=base['electron_absorption_probability']*outgoing;returned=outgoing-absorbed
                nextF[z]+=returned[:,::-1]/m.dz
                amount=float(absorbed.sum());energy=float((absorbed*m.ge.e[:,None]).sum())*EV
                m.electron_out+=amount;m.energy_out+=energy;self.heat[side]+=energy;self.charge[side]-=amount
                self.counts[side,m.ns]+=[outgoing.sum(),0.,0.,amount,returned.sum()]
                for j,state in enumerate(m.states):
                    opts=overrides.get(j,base);outgoing=outP[j][side];incident=float(outgoing.sum())
                    if incident==0:continue
                    sticking=opts['sticking_probability'];reflected=(1-sticking)*outgoing
                    probability=opts['neutralization_probability'] if state['charge']>0 else opts['quenching_probability'] if not state.get('is_ground',False) and not state.get('molecular',False) and state['charge']==0 else 0.
                    target=self.ground[j]
                    if target==j:probability=0.
                    changed=probability*reflected;unchanged=reflected-changed
                    nextP[z,j]+=unchanged[:,::-1]/m.dz;nextP[z,target]+=changed[:,::-1]/m.dz
                    # Signed species crossing columns include state conversions;
                    # this closes nuclei/charge without an algebraic correction.
                    m.species_out[j]+=incident-float(unchanged.sum());m.species_out[target]-=float(changed.sum())
                    grid=m.gh[state['element']]
                    incoming_energy=float((outgoing*(grid.e[:,None]+m.internal[j])).sum())
                    return_energy=float((unchanged*(grid.e[:,None]+m.internal[j])).sum()+(changed*(grid.e[:,None]+m.internal[target])).sum())
                    energy=(incoming_energy-return_energy)*EV;m.energy_out+=energy;self.heat[side]+=energy
                    self.charge[side]+=state['charge']*incident-state['charge']*float(unchanged.sum())-m.states[target]['charge']*float(changed.sum())
                    converted=float(changed.sum());self.counts[side,j]+=[incident,converted if state['charge']==0 else 0.,converted if state['charge']>0 else 0.,sticking*incident,float(reflected.sum())]
            m.P=nextP;m.F=nextF
        m.transport+=m.P.sum(axis=(2,3))-before

    def export(self,folder):
        rows=[]
        if self.enabled:
            for snap in self.model.snapshots:
                for side,name in enumerate(('left','right')):
                    for j,ident in enumerate([s['id'] for s in self.model.states]+['e']):
                        rows.append(dict(time_s=snap['time_s'],side=name,position_cm=0. if side==0 else self.model.length,id=ident,
                            **{key:float(v) for key,v in zip(self.fields,snap['wall_counts'][side,j])},
                            surface_total_heat_erg_cm2=float(snap['wall_heat'][side]),surface_total_charge_e_cm2=float(snap['wall_charge'][side])))
        write_csv(folder/'wall_history.csv',rows,['time_s','side','position_cm','id',*self.fields,'surface_total_heat_erg_cm2','surface_total_charge_e_cm2'])
