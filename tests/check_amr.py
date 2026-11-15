"""Small entry-point checks for the fixed two-level field-aligned solver."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile

p=argparse.ArgumentParser()
p.add_argument("binary",type=Path)
a=p.parse_args()
binary=a.binary.resolve()
env=os.environ.copy()
env.update(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")

with tempfile.TemporaryDirectory(prefix="esf-amr3-") as directory:
    root=Path(directory)
    def run(case,extra):
        out=root/case
        command=[str(binary),"--field","aligned-dipole","--periodic-x","1",
                 "--amr-mode","fixed","--amr-case",case,"--nx","6","--ny","8",
                 "--ns","16","--rtol","1e-10","--pdm-beta","1","--cfl",".25",
                 "--out",str(out)]+extra
        result=subprocess.run(command,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(result.stdout)
        return out
    out=run("potential",["--end","0"])
    m=json.loads((out/"completion.json").read_text())
    assert m["manufactured_potential_max_error"]<1e-8
    assert m["max_potential_relative_residual"]<1.1e-10
    out=run("remap",["--end","0"])
    m=json.loads((out/"remap_checks.json").read_text())
    assert m["max_mass_and_momentum_relative_error"]<1e-11
    assert m["max_constant_or_equilibrium_relative_error"]<1e-10
    assert m["max_constant_velocity_error"]<1e-10
    out=run("translation",["--end","1","--output-every",".5"])
    m=json.loads((out/"completion.json").read_text())
    assert m["constant_tracer_rhs_max"]<1e-10
    assert m["relative_mass_error"]<1e-11 and m["translation_relative_L1"]<.03
    with (out/"interfaces_rank0000.csv").open() as stream:
        interfaces=list(csv.DictReader(stream))
    assert interfaces and all(int(row["left_owned"])==int(row["right_owned"])==1 for row in interfaces)
    with (out/"diagnostics.csv").open() as stream:
        for row in csv.DictReader(stream):
            assert float(row["min_density"])>0
    print("PASS composite tensor potential, physical-volume density/momentum remaps, and 3-D interface transport")
