"""Dark-ionization regressions; synthetic coefficients test operators only."""
from unittest.mock import patch
import unittest
import numpy as np
from atomic_database import HC_EV_CM
from mixture_kinetics import MixtureSlab,ROOT,HPLANCK,EV,isolated_line_profile,conservative_profile
from tests_mixture import config,synthetic_bank


def dark_model(line_eV=1.,ion_edge_eV=4.,electrons=0.,impact=False,normalised_zero=False):
    cfg=config();cfg['initial_cm3']={'Hup':1e8};cfg['electron_initial']['density_cm3']=electrons
    cfg['processes'].update(spontaneous=True,photo_bound_free=True,electron_ionization=impact)
    cfg['radiation'].update(enabled=False,scale=0.)
    if normalised_zero:cfg['radiation'].update(enabled=True,normalization='integral',multiplier_W_cm2=0.,spectrum_file='missing.txt')
    cfg['solver'].update(t_end_s=1e-8,max_step_s=1e-9,min_step_s=1e-18)
    bank,photo=synthetic_bank(upper=True);bank.states[1]['energy_cm1']=ion_edge_eV/HC_EV_CM
    bank.states[4].update(energy_cm1=line_eV/HC_EV_CM,excitation_eV=line_eV)
    bank.lines=[dict(id='SYNTHETIC:line',lower='H0',upper='Hup',A_s=1e7,pump=False,source='SYNTHETIC operator test')]
    for name,edge in [('H0',ion_edge_eV),('Hup',ion_edge_eV-line_eV)]:
        photo.tables[name]=np.array([[edge,1e-18],[100.,1e-18]]);photo.sources[name]='SYNTHETIC operator test'
    edge=ion_edge_eV-line_eV
    bank.cross[('Hup','H1','ionization')]=np.array([[edge,0.],[edge+.001,1e-17],[100.,1e-17]])
    bank.cross_source[('Hup','H1','ionization')]='SYNTHETIC operator test'
    with patch('mixture_kinetics.MixtureData',return_value=bank),patch('mixture_kinetics.PhotoBank',return_value=photo):
        return MixtureSlab(cfg,ROOT)


class DarkIonization(unittest.TestCase):
    def test_nonionizing_decay_has_exactly_zero_ionization_without_electrons(self):
        m=dark_model();field=m.field(0.)['mean_4piJ'][0]
        self.assertGreater(field.max(),0.)  # Internally emitted light remains.
        self.assertEqual(field[m.energy>=3.].max(),0.)
        m.solve()
        self.assertEqual(m.n[0,m.index['H1']],0.);self.assertEqual(m.F.sum(),0.)
        self.assertLess(m.n[0,m.index['Hup']],1e8)
        self.assertLess(abs(m.total_energy()-m.photon_work-m.initial_energy)/m.initial_energy,1e-8)

    def test_zero_normalised_power_does_not_load_file_or_ionize(self):
        m=dark_model(normalised_zero=True);self.assertEqual(m.source.max(),0.)
        m.solve();self.assertEqual(m.n[0,m.index['H1']],0.)

    def test_ionizing_spontaneous_photons_still_ionize_excited_atoms(self):
        m=dark_model(line_eV=4.,ion_edge_eV=6.)
        rate=m.rates(0,m.field(0.)['mean_4piJ'][0])[0]
        k=next(k for k,r in enumerate(m.channels) if r['id']=='PI:Hup')
        self.assertGreater(rate[k],0.)
        m.solve();self.assertGreater(m.n[0,m.index['H1']],0.)

    def test_electron_impact_still_ionizes_with_zero_external_light(self):
        m=dark_model(electrons=1e7,impact=True)
        rate=m.rates(0,np.zeros_like(m.energy))[0]
        k=next(k for k,r in enumerate(m.channels) if r['id']=='ION:Hup')
        self.assertGreater(rate[k],0.)
        m.solve();self.assertGreater(m.n[0,m.index['H1']],0.)

    def test_profile_keeps_count_energy_and_finite_support_on_sparse_continuum(self):
        from radiation import weights
        center=1.;sigma=2e8;gamma=1e6;width=300*sigma*HPLANCK/EV
        e=np.unique(np.r_[np.geomspace(.001,100.,100),center+np.linspace(-width,width,401),center+np.array([-1.,1.])*width*(1+1e-6)])
        w=weights(e*EV/HPLANCK);profile=isolated_line_profile(e,center,sigma,gamma)
        fixed=conservative_profile(profile,w,e,center)
        self.assertEqual(fixed[abs(e-center)>width].max(),0.)
        self.assertAlmostEqual(float(fixed@w),1.,places=13)
        self.assertAlmostEqual(float((fixed*w)@e),center,places=13)


if __name__=='__main__':unittest.main()
