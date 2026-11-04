"""Independent analytic check of the stored real H 2p decay result."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import numpy as np

ROOT=Path(__file__).resolve().parent.parent


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',type=Path,default=ROOT/'outputs/study_01_decay');args=parser.parse_args()
    with np.load(args.folder/'данные_для_графиков.npz',allow_pickle=False) as data:
        meta=json.loads(str(data['metadata_json']));states=[s['id'] for s in meta['states']];t=data['time_s'];n=data['n_cm3'][:,0,states.index('H_0:L003')]
        error=float(np.max(abs(n-1e8*np.exp(-627000000.*t)))/1e8)
        ids=[r['id'] for r in meta['reactions']]
        if ids!=['A:A:H_0:L003:H_0:L001']:raise ValueError('Не изолирован первый тест: '+str(ids))
        if error>1e-8:raise ValueError('Кривая распада не совпала с экспонентой: '+str(error))
        quality=meta['quality'];required=['mass_status','nuclei_status','charge_status','energy_status','kinetic_status','reaction_photon_status']
        if any(quality[k]!='passed' for k in required):raise ValueError('Провален баланс: '+str({k:quality[k] for k in required}))
        version=(ROOT/'VERSION.txt').read_text().strip()
        if meta['version']!=version or quality['version']!=version:raise ValueError('Номер версии метаданных расходится с VERSION.txt')
        report=dict(version=version,python=sys.version.split()[0],platform=platform.platform(),configuration='configs/study_01_decay.json',configuration_sha256=hashlib.sha256((ROOT/'configs/study_01_decay.json').read_bytes()).hexdigest(),maximum_population_error_over_initial=error,upper_population_final_cm3=float(n[-1]),quality=quality,scope='Real stored H coefficients; analytic evolution test, not an experimental validation of A.')
    (ROOT/'validation').mkdir(exist_ok=True);out=ROOT/'validation'/('first_test_'+version.replace('.','_')+'.json');out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(f'T01 passed: {quality["accepted_steps"]} steps, max |n-ref|/n0={error:.3g}; all six balance statuses passed.')


if __name__=='__main__':main()
