"""Exercise slice/movie frame selection and unchanged complete-frame rendering.

Run on allocated compute. Synthetic fixtures and all artifacts/logs are retained
in a new work directory. No model execution is needed for these reader checks.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--baseline-tools',type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    sys.path.insert(0,str(root/'tools'))
    from frame_reader import SavedRun
    args.work.mkdir(parents=True,exist_ok=False)
    commands=[]

    def execute(name,script,*extra,success=True):
        command=[sys.executable,str(script),*map(str,extra)]
        with (args.work/f'{name}.log').open('w') as log:
            result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        commands.append({'name':name,'command':command,'exit_code':result.returncode})
        (args.work/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
        assert (result.returncode==0)==success, commands[-1]

    def fixture(name,night,levels):
        run=args.work/name
        run.mkdir()
        meta={'fine_nx':4,'fine_ny':4,'ns':2,'root_nx':2,'root_ny':2,
              'night_chemistry':night,'amr_physical_levels':levels,'E_peak_m3':1e10,
              'longitude_center_deg':112.25,'longitude_span_deg':8,
              'apex_min_m':250000,'apex_max_m':1600000}
        (run/'metadata.json').write_text(json.dumps(meta)+'\n')
        geometry=np.zeros((4,4,2,6),dtype='<f8')
        geometry[...,3]=110000
        geometry.tofile(run/'geometry_cells.bin')
        for frame in (0,1):
            np.full((4,4,2),.4+frame*.1,dtype='<f8').tofile(run/f'display_density_{frame:04d}.bin')
            if levels==3:
                np.ones((4,4,2),dtype='u1').tofile(run/f'amr_level_{frame:04d}.bin')
        if night:
            np.full((2,4,4,2),.01,dtype='<f8').tofile(run/'molecular_0000.bin')
        else:
            np.ones((2,2,2),dtype='u1').tofile(run/'density_mask_0000.bin')
            (run/'structure_frames.csv').write_text('frame,time_s,mesh_epoch,density_mask\n0,0,0,density_mask_0000.bin\n1,1,1,density_mask_0001.bin\n')
        (run/'frames.csv').write_text('frame,time_s\n0,0\n1,1\n,\n2,nan\n')
        return run

    two=fixture('two-level',False,2)
    night=fixture('night-three-level',True,3)
    inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
            for directory in (two,night) for p in directory.iterdir() if p.is_file()}
    assert [int(x['frame']) for x in SavedRun(two).complete_frames('Oplus')]==[0]
    assert [int(x['frame']) for x in SavedRun(night).complete_frames('Oplus')]==[0,1]
    assert [int(x['frame']) for x in SavedRun(night).complete_frames('electron')]==[0]

    for name,run,quantity,wanted in [('two-fallback',two,'electron',0),
        ('night-Oplus',night,'Oplus',1),('night-electron',night,'electron',0)]:
        prefix=args.work/name
        execute(name,root/'tools/plot_slice.py',run,'--quantity',quantity,'--out',prefix)
        assert json.loads(prefix.with_suffix('.json').read_text())['time_s']==wanted
        movie=args.work/f'{name}-movie'
        execute(f'{name}-movie',root/'tools/make_movie.py',run,'--view','slice',
                '--quantity',quantity,'--frames-only','--out',movie)
        m=json.loads((movie/'movie_manifest.json').read_text())
        expected=[0,1] if name=='night-Oplus' else [0]
        assert [x['frame'] for x in m['selected_frames']]==expected
        assert m['rendered']==len(expected) and m['status']=='frames_complete'
        # Compare complete-frame PNGs against the original plotting code.
        # The original code cannot parse a partial CSV row; use an immutable
        # fixture copy with only complete rows for this paired rendering check.
        complete=args.work/f'{name}-complete-input'
        complete.mkdir()
        for p in run.iterdir():
            if p.name!='frames.csv':
                (complete/p.name).write_bytes(p.read_bytes())
        (complete/'frames.csv').write_text('frame,time_s\n0,0\n')
        old_prefix=args.work/f'{name}-old'
        new_prefix=args.work/f'{name}-new'
        execute(f'{name}-old',args.baseline_tools/'plot_slice.py',complete,'--time','0',
                '--quantity',quantity,'--out',old_prefix)
        execute(f'{name}-new',root/'tools/plot_slice.py',complete,'--time','0',
                '--quantity',quantity,'--out',new_prefix)
        assert old_prefix.with_suffix('.png').read_bytes()==new_prefix.with_suffix('.png').read_bytes()

    execute('refuse-incomplete-two',root/'tools/plot_slice.py',two,'--time','1',
            '--out',args.work/'refused-two',success=False)
    execute('refuse-incomplete-Ne',root/'tools/plot_slice.py',night,'--time','1',
            '--out',args.work/'refused-Ne',success=False)
    protected=args.work/'protected'
    protected.with_suffix('.json').write_text('original provenance\n')
    execute('refuse-provenance-overwrite',root/'tools/plot_slice.py',two,
            '--out',protected,success=False)
    assert protected.with_suffix('.json').read_text()=='original provenance\n'
    assert not protected.with_suffix('.png').exists()
    for name,h in inputs.items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==h
    report={'status':'passed','complete_frame_PNG_parity':'byte-identical for all three paired fixtures',
            'checks':['two-level missing mask fallback','unfinished CSV row ignored',
                      'night Oplus without molecular file','night Ne requires molecular file',
                      'slice/movie selection agreement','exact incomplete time refused',
                      'provenance overwrite refused','input hashes unchanged'],
            'commands':len(commands)}
    (args.work/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    main()
