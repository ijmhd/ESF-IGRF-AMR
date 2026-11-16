"""Targeted AMR source identity, unchanged fields and retry-history regression.

Run on allocated compute. Every command/log/output is retained under a NEW work
directory; fixture edits affect test copies only, never the saved originals.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import time

import numpy as np


def checksum(data):
    value = 14695981039346656037
    for byte in data:
        value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return struct.pack('<Q', value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--mutant', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    args.work.mkdir(parents=True, exist_ok=False)
    records, checks = [], []

    def run(name, case, binary=None, extra=(), success=True, message=None, ranks=1):
        command_args = json.loads((root/'examples'/f'{case}.json').read_text())['arguments']
        for key in ('--coeff', '--neutral-profile'):
            if key in command_args:
                i = command_args.index(key)+1
                command_args[i] = str(root/command_args[i])
        # Explicit qualification inputs, never changes to executable/examples defaults.
        overrides = ['--nx','4','--ny','4','--ns','16','--end','.2',
                     '--output-every','.1','--progress-steps','0',*extra]
        output = args.work/name
        command = ['mpiexec','-n',str(ranks),str(binary or args.binary),
                   *command_args,*overrides,'--out',str(output)]
        record = {'name':name,'command':command,'status':'running'}
        records.append(record)
        (args.work/'commands.json').write_text(json.dumps(records,indent=2)+'\n')
        start = time.monotonic()
        with (args.work/f'{name}.log').open('w') as log:
            process = subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        record.update(exit_code=process.returncode,full_wait_s=time.monotonic()-start,status='complete')
        (args.work/'commands.json').write_text(json.dumps(records,indent=2)+'\n')
        assert (process.returncode == 0) == success, record
        if message:
            assert message in (args.work/f'{name}.log').read_text(), record
        return output

    def edited(name, data):
        path = args.work/f'{name}.bin'
        path.write_bytes(data[:-8]+checksum(data[:-8]))
        return str(path)

    def parity(left, right, tag='0002'):
        errors = {}
        for pattern in ('state','display_density','phi','molecular','molecular_momentum'):
            p, q = left/f'{pattern}_{tag}.bin', right/f'{pattern}_{tag}.bin'
            if not p.exists() and not q.exists():
                continue
            x, y = np.fromfile(p,'<f8'), np.fromfile(q,'<f8')
            assert x.shape == y.shape
            error = float(np.max(np.abs(x-y))/max(float(np.max(np.abs(x))),1e-300))
            assert np.isfinite(error) and error < 1e-10, (pattern,error)
            errors[pattern] = error
        for name in ('amr_level','density_mask'):
            p, q = left/f'{name}_{tag}.bin', right/f'{name}_{tag}.bin'
            if p.exists() or q.exists():
                assert p.read_bytes() == q.read_bytes(), name
        return errors

    source = json.loads((root/'validation/audit_fixes/source_identity.json').read_text())['source_ID']
    for label,case in [('two','hk-amr2'),('three','hk-amr3'),('night','hk-night-demo')]:
        baseline = run(f'{label}-baseline',case,args.baseline)
        current = run(f'{label}-fresh',case)
        errors = parity(baseline,current)
        for path in (baseline,current):
            completion=json.loads((path/'completion.json').read_text())
            assert completion['relative_mass_error'] < 1e-10
            assert completion['max_potential_relative_residual'] <= 1.05e-10
        native = (current/'checkpoint_latest.bin').read_bytes()
        assert native[:8] == b'ESFSRC01' and native[8:72].decode() == source
        assert native[-8:] == checksum(native[:-8])
        extra = ['--restart',str(current/'checkpoint_latest.bin'),'--end','.4']
        if label == 'night':
            extra += ['--night-restart',str(current/'night_state_latest.bin')]
        continued = run(f'{label}-continued',case,extra=extra,ranks=2)
        direct = run(f'{label}-direct',case,extra=['--end','.4'],ranks=2)
        # The final continuation frame is 2, the uninterrupted final frame is 4.
        continuation_errors = {}
        for pattern in ('state','display_density','phi','molecular','molecular_momentum'):
            p,q = continued/f'{pattern}_0002.bin',direct/f'{pattern}_0004.bin'
            if not p.exists() and not q.exists():
                continue
            x,y=np.fromfile(p,'<f8'),np.fromfile(q,'<f8')
            error=float(np.max(np.abs(x-y))/max(float(np.max(np.abs(y))),1e-300))
            assert np.isfinite(error) and error < 1e-9, (label,pattern,error)
            continuation_errors[pattern]=error
        changed = bytearray(native)
        changed[8] = ord('0') if changed[8] != ord('0') else ord('1')
        bad = edited(f'{label}-wrong-ID',changed)
        paired = ['--night-restart',str(current/'night_state_latest.bin')] if label=='night' else []
        run(f'{label}-reject-ID',case,extra=['--restart',bad,*paired],success=False,
            message='numerical source identity mismatch')
        run(f'{label}-reject-source',case,args.mutant,extra=['--restart',str(current/'checkpoint_latest.bin'),*paired],
            success=False,message='numerical source identity mismatch')
        # Actual unbound checkpoints emitted by the preserved baseline must fail.
        run(f'{label}-reject-unbound',case,extra=['--restart',str(baseline/'checkpoint_latest.bin'),*paired],
            success=False,message='lacks numerical source identity')
        damaged=bytearray(native)
        damaged[-9] ^= 1
        damaged_path=args.work/f'{label}-bad-checksum.bin'
        damaged_path.write_bytes(damaged)
        run(f'{label}-reject-checksum',case,extra=['--restart',str(damaged_path),*paired],
            success=False,message='checkpoint checksum')
        if label != 'night':
            changed=bytearray(native)
            struct.pack_into('<Q',changed,72+8+6*8,7)
            copied=edited(f'{label}-retry-history',changed)
            history=run(f'{label}-history',case,extra=['--restart',copied])
            counters=json.loads((history/'performance_counters.json').read_text())
            completion=json.loads((history/'completion.json').read_text())
            assert counters['start_retries']==counters['retries_total']==completion['retries']==7
            assert counters['reported_failure_retries']==counters['segment_failure_retries']==0
            assert counters['timestep_attempts']==0
        checks.append({'case':label,'baseline_field_errors':errors,
                       'continuation_errors':continuation_errors})
        (args.work/'checks_progress.json').write_text(json.dumps(checks,indent=2)+'\n')

    # Cover the legacy non-structure native AMR path as well as structure 2/3.
    legacy=run('legacy-fresh','hk-amr2',extra=['--amr-criteria','legacy','--amr-slab-cells','0','--seed-lon-deg','112.25'])
    run('legacy-continued','hk-amr2',extra=['--amr-criteria','legacy','--amr-slab-cells','0',
        '--seed-lon-deg','112.25','--restart',str(legacy/'checkpoint_latest.bin'),'--end','.3'])
    report={'status':'passed','checks':checks,'source_ID':source,
            'binary_SHA256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (args.binary,args.baseline,args.mutant)}}
    (args.work/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__ == '__main__':
    main()
