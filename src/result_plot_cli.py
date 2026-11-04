"""Optional console counterpart of the selected result checkboxes."""
import argparse
from pathlib import Path
from result_store import ResultData
from result_plots import PlotOptions,build_plots,save_pages,LABELS


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);p.add_argument('--metrics',default='components');p.add_argument('--log-x',action='store_true');p.add_argument('--log-y',action='store_true');p.add_argument('--font-size',type=int,default=16);p.add_argument('--space',action='store_true');p.add_argument('--csv',action='store_true');p.add_argument('--levels',default='');p.add_argument('--components',default='');p.add_argument('--target',default='all');p.add_argument('--sink',action='store_true');p.add_argument('--destination',type=Path);a=p.parse_args(argv)
    metrics=a.metrics.split(',')
    if any(key not in LABELS for key in metrics):p.error('metrics: '+', '.join(LABELS))
    data=ResultData(a.output);options=PlotOptions(metrics=metrics,levels=[s for s in a.levels.split(',') if s] or None,components=[s for s in a.components.split(',') if s] or None,target=a.target,direction='sink' if a.sink else 'source',axis='space' if a.space else 'time',x_scale='log' if a.log_x else 'linear',y_scale='log' if a.log_y else 'linear',fontsize=a.font_size)
    files=save_pages(build_plots(data,options),a.destination or a.output/'рисунки',a.csv)
    print('Сохранено выбранных графиков:',len(files));return 0


if __name__=='__main__':raise SystemExit(main())
