"""Independent checks of the spatial-validation extensions, including real MPI."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('binary',type=Path);p.add_argument('--mpi',action='store_true');p.add_argument('--petsc',action='store_true')
p.add_argument('--skip-time-refinement',action='store_true',help='omit order fitting during performance-only regression checks')
p.add_argument('--work',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
a=p.parse_args();root=Path(__file__).resolve().parents[1];exe=a.binary.resolve()
a.work.mkdir(parents=True,exist_ok=True);a.report.parent.mkdir(parents=True,exist_ok=True)
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
report={'status':'running','binary_sha256':sha(exe),'source_sha256':{str(f.relative_to(root)):sha(f) for d in ('src','include') for f in (root/d).glob('*') if f.is_file()}}
def save():a.report.write_text(json.dumps(report,indent=2)+'\n')
save();oldhook=sys.excepthook
def failed(kind,error,tb):
    report.update(status='failed',exception=kind.__name__,message=str(error));save();oldhook(kind,error,tb)
sys.excepthook=failed
def run(name,options,ranks=1):
    out=a.work/name
    cmd=(['mpiexec','-n',str(ranks)] if ranks>1 else [])+[str(exe),'--rtol','1e-10',*(['--petsc-potential','1'] if a.petsc else []),*map(str,options),'--out',str(out)]
    result=subprocess.run(cmd,text=True,capture_output=True,timeout=300)
    (a.work/(name+'.log')).write_text(result.stdout+'\n'+result.stderr)
    assert result.returncode==0,(name,result.returncode,result.stdout[-1500:],result.stderr[-1500:])
    c=json.loads((out/'completion.json').read_text());assert c['maximum_potential_relative_residual']<1.06e-10
    d=np.atleast_1d(np.genfromtxt(out/'diagnostics.csv',delimiter=',',names=True));assert np.min(d['min_rho'])>=0
    assert np.max(abs(d['relative_mass_error']))<5e-12
    return out
def metadata(out):return json.loads((out/'metadata.json').read_text())
def last(out,kind='state'):
    f=np.atleast_1d(np.genfromtxt(out/'frames.csv',delimiter=',',names=True,dtype=None,encoding='utf8'))[-1]
    return np.fromfile(out/str(f[kind]),dtype='<f8')
options=['--nx',8,'--ny',12,'--ns',48,'--end',120,'--output-every',60,
         '--background-east-mVm',.6,'--seed-east-km',45,'--parallel-order',2,'--parallel-flux','relaxation']
one=run('relaxation1',options)
report['relaxation_coupled']=json.loads((one/'completion.json').read_text())
comparisons=[]
if a.mpi:
    for ranks in (2,4):
        other=run('relaxation'+str(ranks),options,ranks)
        x,y=last(one),last(other);err=float(np.max(abs(x-y))/max(1.,np.max(abs(x))))
        assert err<2e-9;comparisons.append({'ranks':ranks,'relative_state_max_difference':err})
report['relaxation_MPI']=comparisons
periodic=['--field','aligned-dipole','--lon-deg',0,'--span-deg',2.5,'--nx',12,'--ny',20,'--ns',48,
          '--apex-min-km',200,'--apex-max-km',945.968717793,'--background-profile','telf','--E-peak',0,
          '--periodic-x',1,'--open-apex',1,'--open-ends',1,'--parallel-order',2,'--parallel-flux','relaxation',
          '--amplitude',.15,'--seed-width-km',48.083261120685,'--end',200,'--output-every',100]
po=run('periodic1',periodic);m=metadata(po);phi=last(po,'phi').reshape(m['ny']+1,m['nx']+1)
seam=float(np.max(abs(phi[:,0]-phi[:,-1])));assert seam<1e-7
report['periodic_seam_potential_error_V']=seam
if a.mpi:
    for ranks in (2,4):
        other=run('periodic'+str(ranks),periodic,ranks)
        err=float(np.max(abs(last(po)-last(other)))/max(1.,np.max(abs(last(po)))))
        assert err<2e-9;report.setdefault('periodic_MPI',[]).append({'ranks':ranks,'relative_state_max_difference':err})
control=run('periodic_control',periodic+['--amplitude',0]);m=metadata(control)
n=last(control).reshape(2,m['ny'],m['nx'],m['ns'])[0]
err=float(np.max(np.ptp(n,axis=1))/np.max(n));assert err<1e-9
report['unseeded_axisymmetric_longitude_relative_error']=err
coarse=run('calibration_coarse',['--nx',8,'--ny',12,'--ns',32,'--end',0,'--background-east-mVm',.6,'--seed-east-km',45])
fine=run('calibration_fine',['--nx',16,'--ny',24,'--ns',64,'--end',0,'--background-east-mVm',.6,'--seed-east-km',45])
gradients=[metadata(x)['background_potential_gradient_V'] for x in (coarse,fine)]
assert gradients[0]==gradients[1]
report['grid_independent_field_calibration_V']=gradients
aligned=run('calibration_analytic',['--field','aligned-dipole','--lon-deg',0,'--nx',8,'--ny',12,'--ns',32,'--end',0,'--background-east-mVm',.6])
exact=-.0006*np.deg2rad(2.5)*(6371200+350000)
err=abs(metadata(aligned)['background_potential_gradient_V']/exact-1);assert err<1e-8
report['aligned_apex_field_calibration_relative_error']=float(err)
if not a.skip_time_refinement:
    time_runs=[run('ap_dt'+str(i),options+['--end',20,'--output-every',20,'--dt-max',dt]) for i,dt in enumerate((1,.5,.25))]
    m=metadata(time_runs[0]);count=m['nx']*m['ny']*m['ns'];v=np.fromfile(time_runs[0]/'geometry_cells.bin',dtype='<f8').reshape(count,6)[:,5]
    states=[last(x)[:count] for x in time_runs]
    errors=[float(np.sum(abs(states[i]-states[i+1])*v)/np.sum(v)) for i in range(2)]
    order=float(np.log2(errors[0]/errors[1]));assert order>1.7,(order,errors)
    report['relaxation_density_time_refinement']={'successive_L1':errors,'order':order}
else:report['relaxation_density_time_refinement']='omitted for performance-only regression'
report['scope']='smooth/reduced isothermal velocity dynamics; shock-conservative Euler and full nonlinear spatial convergence are separate checks'
report['status']='passed';save();print(json.dumps(report,indent=2))
