"""Energetic level selection by expected reaction visits, with leakage audit.

The full available atomic bank is distinct from the active basis. Each
initially occupied level is a separate seed, protecting trace metastables.
Large A does not hide a short-lived but heavily used intermediate state.
"""
import numpy as np
from scipy.linalg import expm


def positive_metzler_exponential(matrix,dt):
    """Uniformization + scaling/squaring, without clipping evolved values.

    Used when a generic Padé exponential introduces negative roundoff in an
    otherwise positive reaction semigroup. The shifted Taylor matrix is
    nonnegative; 32 terms at shift <= 1/2 are below double precision.
    """
    rate=float(max(0.,-np.diag(matrix).min()))
    if rate==0:return np.eye(len(matrix))+matrix*dt
    scaling=max(0,int(np.ceil(np.log2(max(rate*dt/.5,1.)))))
    small=dt/2.**scaling;shift=rate*small
    x=matrix*small+np.eye(len(matrix))*shift
    term=np.eye(len(matrix));result=term.copy()
    for k in range(1,33):
        term=term@x/k;result+=term
    result*=np.exp(-shift)
    for _ in range(scaling):result=result@result
    return result


def frozen_evolution(initial,channels,rates,dt):
    """Positive frozen one-heavy-particle network, exact residence integrals.

    channels contains i,j. Initial can be [states, phase_nodes]. One block
    exponential supplies the same residence times to all reaction diagnostics.
    """
    initial=np.asarray(initial,dtype=float);states=initial.shape[0]
    if np.any(np.asarray(rates)<0) or not np.all(np.isfinite(rates)):raise ValueError('Неверные частоты реакций')
    Q=np.zeros((states,states))
    for r,k in zip(channels,rates):
        Q[r['j'],r['i']]+=k;Q[r['i'],r['i']]-=k
    matrix=np.zeros((2*states,2*states));matrix[:states,:states]=Q;matrix[states:,:states]=np.eye(states)
    e=expm(matrix*dt)
    if np.any(e<0):e=positive_metzler_exponential(matrix,dt)
    final=e[:states,:states]@initial;residence=e[states:,:states]@initial
    if not np.all(np.isfinite(final)) or not np.all(np.isfinite(residence)) or np.any(final<0) or np.any(residence<0):raise ValueError(f'Отрицательная/неконечная матричная эволюция: min(n)={final.min():.3e}, min(residence)={residence.min():.3e}, max(rate)*dt={max(rates,default=0)*dt:.3e}')
    return final,residence


class LevelSelection:
    def __init__(self,states,channels,initial,rates,end,options=None):
        self.states=states;self.channels=channels;self.options=options or {}
        if set(self.options)-{'mode','relative_tolerance','forced_levels','max_restarts','missing_data'}:
            raise ValueError('Неизвестные level_selection')
        self.tolerance=float(self.options.get('relative_tolerance',1e-6))
        if not 0<self.tolerance<.1:raise ValueError('level_selection.relative_tolerance: (0,0.1)')
        self.mode=self.options.get('mode','auto')
        if self.mode not in ('auto','all'):raise ValueError('level_selection.mode: auto или all')
        seeds=np.flatnonzero(np.max(initial,axis=0)>0)
        self.reference={}
        for element in sorted({s['element'] for s in states}):
            values=initial[:,[i for i,s in enumerate(states) if s['element']==element]]
            positive=values[values>0];self.reference[element]=float(positive.min()) if len(positive) else 1.
        self.active=np.zeros(len(states),dtype=bool);self.active[seeds]=True
        for name in self.options.get('forced_levels',[]):
            ids=[i for i,s in enumerate(states) if s['id']==name]
            if not ids:raise ValueError('Неизвестный forced_level '+name)
            self.active[ids]=True
        self.visit=np.zeros(len(states));self.omitted_integrals=np.zeros(len(states));self.restarts=0
        if self.mode=='all':self.active[:]=True
        elif len(seeds):
            # Scaled seed columns avoid abundance-dependent pruning.
            seed=np.eye(len(states))[:,seeds]
            for cell_rates in rates:
                final,residence=frozen_evolution(seed,channels,cell_rates,end)
                visits=seed.copy()
                for r,k in zip(channels,cell_rates):visits[r['j']]+=k*residence[r['i']]
                self.visit=np.maximum(self.visit,np.maximum(final.max(axis=1),visits.max(axis=1)))
            self.active|=self.visit>=self.tolerance*.1
        self.close_radiative_cascades()

    def close_radiative_cascades(self):
        # Keep every known spontaneous descendant of an active level, even
        # when it has a tiny instantaneous population.
        changed=True
        while changed:
            changed=False
            for r in self.channels:
                if r['kind'] in ('spontaneous','two_photon') and self.active[r['i']] and not self.active[r['j']]:
                    self.active[r['j']]=True;changed=True

    def channel_mask(self):
        return np.array([self.active[r['i']] and self.active[r['j']] for r in self.channels])

    def leakage(self,populations,rates,dt):
        additions=np.zeros(len(self.states))
        for r,k in zip(self.channels,rates):
            if self.active[r['i']] and not self.active[r['j']]:additions[r['j']]+=k*populations[r['i']]*dt
        return additions

    def candidates(self,integrals):
        return np.array([not self.active[i] and integrals[i]>self.tolerance*self.reference[s['element']]
                         for i,s in enumerate(self.states)])

    def expand(self,mask):
        self.active|=mask;self.close_radiative_cascades();self.omitted_integrals[:]=0.;self.restarts+=1
        if self.restarts>self.options.get('max_restarts',12):raise RuntimeError('Слишком много расширений атомного базиса; используйте level_selection.mode=all')

    def report(self):
        return [dict(id=s['id'],label=s['label'],element=s['element'],charge=s['charge'],n=s.get('n',''),
            excitation_eV=s['excitation_eV'],active=bool(self.active[i]),predicted_visits_per_seed=float(self.visit[i]),
            omitted_inflow_integral_cm3=float(self.omitted_integrals[i]),
            omitted_fraction_of_smallest_seed=float(self.omitted_integrals[i]/self.reference[s['element']]))
                for i,s in enumerate(self.states)]
