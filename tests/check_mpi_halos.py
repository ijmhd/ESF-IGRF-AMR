"""Exercise ownership cuts through rows and halos spanning several MPI peers."""
import argparse,hashlib,json,subprocess
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('binary',type=Path);p.add_argument('--work',type=Path,required=True);p.add_argument('--report',type=Path,required=True);p.add_argument('--ranks',type=int,nargs='+',default=[1,3,7]);a=p.parse_args()
a.work.mkdir(exist_ok=False,parents=True);report={'status':'running','binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest(),'runs':[]}
a.report.parent.mkdir(exist_ok=True,parents=True)
def save():a.report.write_text(json.dumps(report,indent=2)+'\n')
save()
try:
 for periodic in (False,True):
  reference=None
  for petsc in (0,1):
   for ranks in a.ranks:
    out=a.work/f'periodic{int(periodic)}_petsc{petsc}_p{ranks}'
    options=['--rtol','1e-11','--nx','10','--ny','13','--ns','32','--end','20','--output-every','20',
     '--parallel-flux','relaxation','--parallel-order','2','--petsc-potential',str(petsc)]
    if periodic:options+=['--field','aligned-dipole','--lon-deg','0','--apex-min-km','200','--apex-max-km','945.968717793','--background-profile','telf','--E-peak','0','--periodic-x','1','--open-apex','1','--open-ends','1']
    else:options+=['--background-east-mVm','.6','--seed-east-km','45']
    cmd=(['mpiexec','-n',str(ranks)] if ranks>1 else [])+[str(a.binary.resolve()),*options,'--out',str(out)]
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=90);(a.work/(out.name+'.log')).write_text(r.stdout+'\n'+r.stderr);assert r.returncode==0,(out.name,r.stderr[-1000:])
    c=json.loads((out/'completion.json').read_text());assert c['maximum_potential_relative_residual']<1.06e-11
    state=np.fromfile(out/'state_0001.bin',dtype='<f8');phi=np.fromfile(out/'phi_0001.bin',dtype='<f8')
    if reference is None:reference=(state.copy(),phi.copy())
    errors=[float(np.max(abs(x-y))/max(1.,np.max(abs(y)))) for x,y in zip((state,phi),reference)]
    assert max(errors)<2e-8,(out.name,errors)
    if not petsc:assert np.array_equal(state,reference[0]),(out.name,'native state changed')
    d=np.atleast_1d(np.genfromtxt(out/'diagnostics.csv',delimiter=',',names=True));budget=float(max(abs(d['relative_mass_error'])));assert budget<5e-12
    report['runs'].append({'periodic':periodic,'petsc':bool(petsc),'ranks':ranks,'state_and_potential_relative_max_difference':errors,'native_state_bitwise_equal':bool(np.array_equal(state,reference[0])),'maximum_mass_budget_relative_error':budget});save()
 report['status']='passed'
except BaseException as e:report.update(status='failed',exception=type(e).__name__,message=str(e));raise
finally:save()
print(json.dumps(report,indent=2))
