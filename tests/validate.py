"""Coupled IGRF invariants, time/geometry refinement and real MPI checks."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('--mpi',action='store_true')
p.add_argument('--work',required=True,type=Path);p.add_argument('--report',required=True,type=Path);a=p.parse_args()
root=Path(__file__).resolve().parents[1];exe=a.binary.resolve();a.work.mkdir(parents=True,exist_ok=True);a.report.parent.mkdir(parents=True,exist_ok=True)
report={'status':'running','binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),
 'source_sha256':{str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for d in ('src','include') for f in (root/d).glob('*') if f.is_file()}}
def save():a.report.write_text(json.dumps(report,indent=2)+'\n')
save();old_hook=sys.excepthook
def failed(kind,error,traceback):
 report.update(status='failed',exception=kind.__name__,message=str(error));save();old_hook(kind,error,traceback)
sys.excepthook=failed
def run(name,*options,ranks=1):
 out=a.work/name
 cmd=(['mpiexec','-n',str(ranks)] if ranks>1 else [])+[str(exe),'--rtol','1e-10',*map(str,options),'--out',str(out)]
 result=subprocess.run(cmd,text=True,capture_output=True,timeout=180)
 (a.work/(name+'.log')).write_text(result.stdout+'\n'+result.stderr)
 assert result.returncode==0,(name,result.returncode,result.stdout[-2000:],result.stderr[-2000:])
 completion=json.loads((out/'completion.json').read_text());assert completion['maximum_potential_relative_residual']<1.06e-10
 d=np.atleast_1d(np.genfromtxt(out/'diagnostics.csv',delimiter=',',names=True));assert np.min(d['min_rho'])>=0
 assert np.max(abs(d['relative_mass_error']))<5e-12
 return out
def final(out):
 frames=np.atleast_1d(np.genfromtxt(out/'frames.csv',delimiter=',',names=True,dtype=None,encoding='utf8'))
 return np.fromfile(out/str(frames[-1]['state']),dtype='<f8')
def columns(out):
 frames=np.atleast_1d(np.genfromtxt(out/'frames.csv',delimiter=',',names=True,dtype=None,encoding='utf8'))
 return np.atleast_1d(np.genfromtxt(out/str(frames[-1]['columns']),delimiter=',',names=True))
options=['--nx',8,'--ny',12,'--ns',48,'--end',120,'--output-every',60]
one=run('coupled1',*options)
report['coupled']=json.loads((one/'completion.json').read_text())
report['geometry']=json.loads((one/'metadata.json').read_text())['geometry_quality']
if a.mpi:
 comparisons=[]
 for ranks in (2,4):
  other=run('coupled'+str(ranks),*options,ranks=ranks);x,y=final(one),final(other);assert x.shape==y.shape
  err=float(np.max(abs(x-y))/max(1.,np.max(abs(x))));assert err<2e-9,err
  # Geometry is broadcast from rank zero; changing ownership cannot change it.
  assert (one/'geometry_cells.bin').read_bytes()==(other/'geometry_cells.bin').read_bytes()
  comparisons.append({'ranks':ranks,'relative_state_max_difference':err})
 report['MPI_comparisons']=comparisons
deep=run('deep',*options,'--amplitude',.99,'--seed-width-km',3000)
baseE=columns(one)['pedersen_E_S'];deepE=columns(deep)['pedersen_E_S']
assert np.array_equal(baseE,deepE) and np.min(deepE)>0
report['deep_depletion_fixed_E']={'amplitude':.99,'minimum_E_line_Pedersen_S':float(np.min(deepE)),'E_conductance_unchanged':True}
coarse=run('geometry_coarse','--nx',6,'--ny',8,'--ns',32,'--end',0)
fine=run('geometry_fine','--nx',12,'--ny',16,'--ns',64,'--end',0)
cg=json.loads((coarse/'metadata.json').read_text())['geometry_quality'];fg=json.loads((fine/'metadata.json').read_text())['geometry_quality']
for key in ('field_direction_sine_error_max','adjacent_parallel_face_magnetic_flux_relative_jump_max'):
 assert fg[key]<.7*cg[key],(key,cg[key],fg[key])
report['geometry_refinement']={'coarse':cg,'fine':fg}
time_runs=[run('dt'+str(i),'--nx',8,'--ny',12,'--ns',48,'--end',20,'--output-every',20,'--dt-max',dt) for i,dt in enumerate((1,.5,.25))]
m=json.loads((time_runs[0]/'metadata.json').read_text());count=m['nx']*m['ny']*m['ns'];volume=np.fromfile(time_runs[0]/'geometry_cells.bin',dtype='<f8').reshape(count,6)[:,5]
states=[final(d)[:count] for d in time_runs]
errors=[float(np.sum(abs(states[i]-states[i+1])*volume)/np.sum(volume)) for i in range(2)]
order=float(np.log2(errors[0]/errors[1]));assert order>1.7,(order,errors)
report['ETD2_temporal_density_refinement']={'dt_s':[1,.5,.25],'successive_volume_weighted_L1':errors,'observed_order':order}
driven=run('driven',*options,'--background-east-mVm',.6,'--seed-east-km',45)
driven_diag=np.atleast_1d(np.genfromtxt(driven/'diagnostics.csv',delimiter=',',names=True))
assert np.max(abs(driven_diag['integrated_boundary_particle_flux']))>0
report['imposed_background_E_open_budget']={'east_mVm':.6,'maximum_relative_particle_budget_error':float(np.max(abs(driven_diag['relative_mass_error']))),'initial_max_up_m_s':float(driven_diag[0]['max_up_m_s'])}
if a.mpi:
 driven4=run('driven4',*options,'--background-east-mVm',.6,'--seed-east-km',45,ranks=4)
 x,y=final(driven),final(driven4);err=float(np.max(abs(x-y))/max(1.,np.max(abs(x))));assert err<2e-9
 report['imposed_background_E_open_budget']['MPI4_relative_state_max_difference']=err
report['status']='passed';save();print(json.dumps(report,indent=2))
