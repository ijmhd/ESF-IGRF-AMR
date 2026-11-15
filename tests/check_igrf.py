"""Independent official Fortran comparison; Fortran is a test dependency only.

The unmodified reference contains rounded/single-precision numeric literals.
A separately labeled precision-normalized reference isolates that difference.
"""
import argparse,hashlib,json,math,random,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True)
p.add_argument('--official-fortran',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
p.add_argument('--report',type=Path,required=True);a=p.parse_args()
raw=a.official_fortran.read_bytes()
assert hashlib.sha256(raw).hexdigest()=='0c4fe0664d837b86c1d1fb4d2a354d6a14a29bbf028b550ab1342e1b828dd2c0'
s=raw.decode();s=s[s.index('      subroutine igrf14syn'):]
b=a.build.resolve();b.mkdir(parents=True,exist_ok=True)
(b/'driver.f90').write_text('''program check
 implicit none
 integer :: ios
 real(8) :: year,lat,lon,alt,x,y,z,f
 do
 read(*,*,iostat=ios) year,lat,lon,alt
 if(ios/=0)exit
 call igrf14syn(0,year,1,alt,90d0-lat,lon,x,y,z,f)
 write(*,'(4ES25.16)') x,y,z,f
 end do
end program
''')
precise=s.replace('0.017453292','0.017453292519943295d0').replace('40680631.6',format(6378.137**2,'.16g')+'d0').replace('40408296.0',format((6378.137*(1-1/298.257223563))**2,'.16g')+'d0')
cases=[(2025.,22.3,114.2,350.),(2020.,0.,0.,0.),(2030.,-45.,240.,1200.),(1900.,60.,90.,100.)]
rng=random.Random(20261110)
for i in range(60):cases.append((rng.uniform(1900,2030),rng.uniform(-89.99,89.99),rng.uniform(0,360),rng.uniform(0,2000)))
coeff=Path(__file__).resolve().parents[1]/'data/igrf14coeffs.txt'
actual=[list(map(float,subprocess.check_output([str(a.binary.resolve()),str(coeff),*map(str,c),'13'],text=True).split())) for c in cases]
report={'points':len(cases),'reference_fortran_sha256':hashlib.sha256(raw).hexdigest(),'coefficient_sha256':hashlib.sha256(coeff.read_bytes()).hexdigest()}
for name,code,flags,limit in [('official_unmodified',s,[],.05),('precision_normalized',precise,['-fdefault-real-8','-fdefault-double-8'],1e-5)]:
 source=b/(name+'.f');source.write_text(code);exe=b/name
 subprocess.run(['gfortran','-O2','-std=legacy',*flags,str(source),str(b/'driver.f90'),'-o',str(exe)],check=True)
 result=subprocess.run([str(exe)],input=''.join(' '.join(map(str,c))+'\n' for c in cases),text=True,capture_output=True,check=True)
 reference=[list(map(float,line.split())) for line in result.stdout.splitlines() if len(line.split())==4]
 assert len(reference)==len(cases)
 error=max(abs(x-y) for row,ref in zip(actual,reference) for x,y in zip(row,ref))
 report[name]={'maximum_difference_nT':error,'threshold_nT':limit,'passed':error<limit}
report['precision_normalization']='Only radians-per-degree and WGS84 squared axes replaced by full precision; real literals compiled in double precision. This second reference is explicitly modified.'
report['hong_kong_2025_350km_ned_nT']=actual[0]
report['status']='passed' if all(report[k]['passed'] for k in ('official_unmodified','precision_normalized')) else 'failed'
a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));assert report['status']=='passed'
