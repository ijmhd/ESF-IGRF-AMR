#!/usr/bin/env python3
"""Plot a complete saved transverse slice and actual density levels, without model execution.

This is a fixed field-label-index slice, not a geographic horizontal cut.
Input arrays are decoded only after metadata/frame/size checks. Displayed fine
variation in a coarse O+ region is conservative reconstruction, not an evolved DOF.
"""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.collections import LineCollection
from frame_reader import SavedRun, sha256 as sha

def main():
    """Select a complete frame, decode normalized fields and render two linked panels."""
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path);p.add_argument('--time',type=float)
    p.add_argument('--field-index',type=int)
    p.add_argument('--quantity',choices=['electron','Oplus'],default='electron')
    p.add_argument('--out',type=Path,required=True,help='New output prefix without suffix')
    a=p.parse_args()
    for suffix in ('.pdf','.png','.json'):
        if a.out.with_suffix(suffix).exists():raise FileExistsError(a.out.with_suffix(suffix))
    run=SavedRun(a.run);root=run.root;meta=run.meta
    fx,fy,ns=run.fx,run.fy,run.ns
    row=run.select(a.time,a.quantity)
    k=ns//2 if a.field_index is None else a.field_index
    if not 0<=k<ns:raise ValueError('field-index out of range')
    density,level,files=run.read(row,a.quantity)
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
    a.out.with_suffix('.json').write_text(json.dumps({'time_s':time,'field_index':k,'quantity':a.quantity,'slice_coordinate':'fixed field label; apex-longitude/apex-altitude axes','AMR_outline':'boundaries between actual physical levels; no individual cell grid','display_limits_normalized':[0,1],'density_floor_in_model_added':False,'input_SHA256':{p.name:sha(p) for p in files},'plotter_SHA256':sha(Path(__file__)),'reader_SHA256':sha(Path(__file__).with_name('frame_reader.py')),'numpy_version':np.__version__,'matplotlib_version':matplotlib.__version__},indent=2)+'\n')

if __name__=='__main__':main()
