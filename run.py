"""Single project entry point: model, GUI, data import, plotting or tests."""
from pathlib import Path
import os
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))


def main():
    args=sys.argv[1:]
    if args and args[0]=='data':
        from data_import import main as execute
        return execute(args[1:])
    if args and args[0]=='plot':
        if len(args)>1 and (Path(args[1])/'данные_для_графиков.npz').is_file():
            from result_plot_cli import main as execute
            return execute(args[1:])
        return subprocess.call([sys.executable,str(ROOT/'src/plot_mixture.py'),*args[1:]],cwd=ROOT)
    if args and args[0]=='test':
        environment=dict(os.environ)
        environment['PYTHONPATH']=os.pathsep.join([str(ROOT/'src'),str(ROOT/'tests'),environment.get('PYTHONPATH','')])
        return subprocess.call([sys.executable,'-m','unittest','discover','-s',str(ROOT/'tests'),'-p','tests_*.py','-v'],cwd=ROOT,env=environment)
    if args==['--gui']:
        from launcher import main as execute
        execute();return 0
    from mixture_kinetics import main as execute
    return execute(args)


if __name__=='__main__':
    try:sys.exit(main())
    except (ValueError,OSError,KeyError,RuntimeError) as exc:sys.exit('Ошибка: '+str(exc))
