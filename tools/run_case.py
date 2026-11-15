#!/usr/bin/env python3
"""Launch an explicit, checked-in case; no hidden physical defaults are added."""
import argparse,json,os,shlex,subprocess
from pathlib import Path

def main():
    """Expand a checked-in effective case, apply explicit overrides and launch MPI."""
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('case',choices=['smoke','dipole','hk-amr2','hk-amr3','hk-night-demo'])
    p.add_argument('--binary',type=Path,default=Path('build/esf3d'))
    p.add_argument('--ranks',type=int,default=1)
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--mpiexec',default='mpiexec')
    p.add_argument('--end',type=float,help='Absolute model stop time in seconds, not continuation duration')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--restart',type=Path)
    p.add_argument('--night-restart',type=Path)
    p.add_argument('--dry-run',action='store_true')
    a=p.parse_args();root=Path(__file__).resolve().parents[1]
    if a.ranks<1 or a.threads<1:p.error('ranks and threads must be positive')
    c=json.loads((root/'examples'/(a.case+'.json')).read_text());args=c['arguments'].copy()
    def put(key,value):
        """Replace one existing CLI value or append one explicitly requested option."""
        if key in args:args[args.index(key)+1]=str(value)
        else:args.extend([key,str(value)])
    for key in ['--coeff','--neutral-profile','--collision-profile','--electric-driver']:
        if key in args:
            i=args.index(key)+1;args[i]=str((root/args[i]).resolve())
    put('--out',a.out.resolve())
    if a.end is not None:put('--end',a.end)
    if a.restart:put('--restart',a.restart.resolve())
    if a.night_restart:put('--night-restart',a.night_restart.resolve())
    if bool(a.night_restart)!=bool(a.restart) and '--night-chemistry' in args and args[args.index('--night-chemistry')+1]=='1':p.error('Night continuation requires BOTH native and molecular checkpoints')
    maxr=int(args[args.index('--ny')+1]);
    if a.ranks>maxr:p.error('Row decomposition supports at most the physical root ny ranks')
    cmd=[a.mpiexec,'-n',str(a.ranks),str(a.binary.resolve()),*args]
    env=os.environ.copy();env.update(OMP_NUM_THREADS=str(a.threads),OMP_DYNAMIC='false',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    print(shlex.join(cmd),flush=True)
    print('MPI ranks =',a.ranks,'OMP threads/rank =',a.threads,'case =',c['description'],flush=True)
    if a.dry_run:return
    if not a.binary.is_file():p.error('Build the executable first')
    if (a.out/'metadata.json').exists():p.error('Choose a NEW output directory; existing runs are preserved')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    subprocess.run(cmd,env=env,check=True)

if __name__=='__main__':main()
