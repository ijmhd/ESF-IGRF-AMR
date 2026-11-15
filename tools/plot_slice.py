#!/usr/bin/env python3
"""Plot a complete saved transverse slice and actual density levels, without model execution.

This is a fixed field-label-index slice, not a geographic horizontal cut.
Input arrays are decoded only after metadata/frame/size checks. Displayed fine
variation in a coarse O+ region is conservative reconstruction, not an evolved DOF.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.collections import LineCollection

def sha(path):
    """Return the exact SHA256 of a plotted input for figure provenance."""
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def main():
    """Select a complete frame, decode normalized fields and render two linked panels."""
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path);p.add_argument('--time',type=float)
    p.add_argument('--field-index',type=int)
    p.add_argument('--quantity',choices=['electron','Oplus'],default='electron')
    p.add_argument('--out',type=Path,required=True,help='New output prefix without suffix')
    a=p.parse_args();root=a.run;meta=json.loads((root/'metadata.json').read_text())
    fx,fy,ns=[int(meta[k]) for k in ['fine_nx','fine_ny','ns']];size=fx*fy*ns
    frames=list(csv.DictReader((root/'frames.csv').open()));records=[]
    for row in frames:
        num=int(row['frame']);tag=f'{num:04d}';display=root/f'display_density_{tag}.bin'
        mol=root/f'molecular_{tag}.bin';lev=root/f'amr_level_{tag}.bin'
        if not display.exists() or display.stat().st_size!=8*size:continue
        if meta['night_chemistry'] and (not mol.exists() or mol.stat().st_size!=16*size):continue
        if int(meta['amr_physical_levels'])==3 and (not lev.exists() or lev.stat().st_size!=size):continue
        records.append(row)
    if a.time is not None:records=[r for r in records if abs(float(r['time_s'])-a.time)<1e-8]
    if not records:raise ValueError('No complete saved frame at requested time')
    row=max(records,key=lambda r:float(r['time_s']));num=int(row['frame']);tag=f'{num:04d}'
    k=ns//2 if a.field_index is None else a.field_index
    if not 0<=k<ns:raise ValueError('field-index out of range')
    files=[root/'metadata.json',root/'frames.csv',root/f'display_density_{tag}.bin']
    density=np.fromfile(files[-1],'<f8').reshape(fy,fx,ns)
    if a.quantity=='electron':
        if meta['night_chemistry']:
            f=root/f'molecular_{tag}.bin';density+=np.fromfile(f,'<f8').reshape(2,fy,fx,ns).sum(axis=0);files.append(f)
        else:
            f=root/'geometry_cells.bin';g=np.memmap(f,'<f8',mode='r',shape=(fy,fx,ns,6));z=(g[...,3]-110000.)/12000.
            density+=meta['E_peak_m3']/1e12*np.exp(1-z-np.exp(-z));files.append(f)
    lev=root/f'amr_level_{tag}.bin'
    if lev.exists():level=np.fromfile(lev,'u1').reshape(fy,fx,ns);files.append(lev)
    else:
        structure=list(csv.DictReader((root/'structure_frames.csv').open()));s=next(r for r in structure if int(r['frame'])==num)
        f=root/s['density_mask'];mask=np.fromfile(f,'u1').reshape(meta['root_ny'],meta['root_nx'],ns)
        level=mask.repeat(2,axis=0).repeat(2,axis=1);files.extend([root/'structure_frames.csv',f])
    if not np.isfinite(density).all() or np.any(density<0):raise ValueError('Invalid saved density')
    lon=meta['longitude_center_deg'];span=meta['longitude_span_deg'];h0=meta['apex_min_m']/1000;h1=meta['apex_max_m']/1000
    x=lon+span*((np.arange(fx)+.5)/fx-.5);y=h0+(h1-h0)*(np.arange(fy)+.5)/fy
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),constrained_layout=True)
    extent=[lon-span/2,lon+span/2,h0,h1]
    im=axes[0].imshow(density[:,:,k],origin='lower',extent=extent,aspect='auto',vmin=0,vmax=1,cmap='viridis',interpolation='nearest')
    # Actual digital cell-edge segments, not interpolated contour curves and
    # not every fine-grid line. Adjacent equal-level samples have no outline.
    xe=np.linspace(extent[0],extent[1],fx+1);ye=np.linspace(h0,h1,fy+1);lev2=level[:,:,k];segments=[]
    for j,i in zip(*np.where(lev2[:,1:]!=lev2[:,:-1])):segments.append([(xe[i+1],ye[j]),(xe[i+1],ye[j+1])])
    for j,i in zip(*np.where(lev2[1:,:]!=lev2[:-1,:])):segments.append([(xe[i],ye[j+1]),(xe[i+1],ye[j+1])])
    axes[0].add_collection(LineCollection(segments,colors='white',linewidths=.65))
    axes[0].set_title('Electron density' if a.quantity=='electron' else 'O+ density')
    fig.colorbar(im,ax=axes[0],label=r'$n/(10^{12}\,\mathrm{m}^{-3})$')
    cm=ListedColormap(['#e7eaf0','#5c89b5','#e18e42']);norm=BoundaryNorm([-.5,.5,1.5,2.5],cm.N)
    im2=axes[1].imshow(level[:,:,k],origin='lower',extent=extent,aspect='auto',cmap=cm,norm=norm,interpolation='nearest')
    axes[1].set_title('Actual physical density refinement level');cb=fig.colorbar(im2,ax=axes[1],ticks=[0,1,2]);cb.ax.set_yticklabels(['L0','L1','L2'])
    for ax in axes:ax.set_xlabel('Apex geographic longitude (degrees E)');ax.set_ylabel('Apex spherical altitude (km)')
    time=float(row['time_s']);fig.suptitle(f't = {time:g} s; fixed field-label index k={k} of {ns}\nNot a constant geographic altitude or equatorial slice',fontsize=11)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    for suffix in ['.pdf','.png']:
        out=a.out.with_suffix(suffix)
        if out.exists():raise FileExistsError(out)
        fig.savefig(out,dpi=220)
    plt.close(fig)
    a.out.with_suffix('.json').write_text(json.dumps({'time_s':time,'field_index':k,'quantity':a.quantity,'slice_coordinate':'fixed field label; apex-longitude/apex-altitude axes','AMR_outline':'boundaries between actual physical levels; no individual cell grid','display_limits_normalized':[0,1],'density_floor_in_model_added':False,'input_SHA256':{p.name:sha(p) for p in files},'plotter_SHA256':sha(Path(__file__)),'numpy_version':np.__version__,'matplotlib_version':matplotlib.__version__},indent=2)+'\n')

if __name__=='__main__':main()
