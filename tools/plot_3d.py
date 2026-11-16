#!/usr/bin/env python3
"""Draw three curved grid-label cuts, actual AMR outlines and traced field lines.

These are geometry-aware cuts through one complete saved frame, not geographic
constant-altitude/ISR cuts. All positions come from nodes.bin; this program does
not re-trace or approximate the IGRF. See docs/FROM_SOURCE_TO_MOVIE_ZH.md.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from frame_reader import SavedRun, sha256


def edge_segments(vertices, levels):
    """Keep exact cell edges where actual AMR levels differ, without cell grids.

    vertices has one extra row/column. The result uses every saved level edge
    even when surface colours are subsampled for display with --stride.
    """
    segments = []
    for j, i in zip(*np.where(levels[:, 1:] != levels[:, :-1])):
        segments.append([vertices[j, i + 1], vertices[j + 1, i + 1]])
    for j, i in zip(*np.where(levels[1:, :] != levels[:-1, :])):
        segments.append([vertices[j + 1, i], vertices[j + 1, i + 1]])
    return segments


def main():
    """Validate inputs, transform ECEF geometry and save figure plus provenance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--time', type=float)
    parser.add_argument('--quantity', choices=['electron', 'Oplus'], default='electron')
    parser.add_argument('--i', type=int, help='Fixed longitude-label cell index')
    parser.add_argument('--j', type=int, help='Fixed apex-height-label cell index')
    parser.add_argument('--k', type=int, help='Fixed along-field-label cell index')
    parser.add_argument('--stride', type=int, default=1, help='Colour display sampling only')
    parser.add_argument('--out', type=Path, required=True, help='New figure prefix')
    args = parser.parse_args()
    if args.stride < 1:
        parser.error('stride must be positive')
    for suffix in ('.png', '.pdf', '.json'):
        if args.out.with_suffix(suffix).exists():
            raise FileExistsError(args.out.with_suffix(suffix))
    run = SavedRun(args.run)
    row = run.select(args.time, args.quantity)
    density, level, inputs = run.read(row, args.quantity)
    nodes, node_path = run.nodes()
    i = run.fx // 2 if args.i is None else args.i
    j = run.fy // 2 if args.j is None else args.j
    k = run.ns // 2 if args.k is None else args.k
    if not (0 <= i < run.fx and 0 <= j < run.fy and 0 <= k < run.ns):
        parser.error('i/j/k must be valid canonical cell indices (zero based)')

    # One fixed local Cartesian basis makes curved dipole/IGRF geometry visible.
    # The vertical axis is displacement along reference radial up, NOT altitude.
    reference = np.array(nodes[run.fy // 2, run.fx // 2, run.ns // 2])
    up = reference / np.linalg.norm(reference)
    east = np.cross([0., 0., 1.], up)
    east /= np.linalg.norm(east)
    north = np.cross(up, east)
    basis = np.stack((east, north, up))

    def project(points):
        """Convert ECEF metres to east/north/radial-up km about the reference."""
        result = (points - reference) @ basis.T / 1000.
        if not np.isfinite(result).all():
            raise ValueError('Nonfinite saved geometry')
        return result

    # Midpoint vertex surfaces and cell values have exactly compatible shapes.
    # i,j,k are cell-label cuts; they are not planes of constant geographic height.
    cuts = [
        (project(.5 * (nodes[:, i, :, :] + nodes[:, i + 1, :, :])),
         density[:, i, :], level[:, i, :]),
        (project(.5 * (nodes[j, :, :, :] + nodes[j + 1, :, :, :])),
         density[j, :, :], level[j, :, :]),
        (project(.5 * (nodes[:, :, k, :] + nodes[:, :, k + 1, :])),
         density[:, :, k], level[:, :, k]),
    ]
    fig = plt.figure(figsize=(12, 9), constrained_layout=True)
    ax = fig.add_subplot(111, projection='3d')
    cmap = plt.get_cmap('viridis')
    norm = Normalize(0, 1)
    bounds = []
    for vertices, values, labels in cuts:
        # Retain both endpoints. Colours sample the leading cell of each display
        # patch; --stride never changes scientific data or actual-level outlines.
        rr = np.unique(np.r_[np.arange(0, values.shape[0], args.stride), values.shape[0]])
        cc = np.unique(np.r_[np.arange(0, values.shape[1], args.stride), values.shape[1]])
        xyz = vertices[np.ix_(rr, cc)]
        colours = cmap(norm(values[np.ix_(rr[:-1], cc[:-1])]))
        ax.plot_surface(xyz[..., 0], xyz[..., 1], xyz[..., 2], facecolors=colours,
                        rstride=1, cstride=1, linewidth=0, antialiased=False,
                        shade=False, alpha=.94)
        edges = edge_segments(vertices, labels)
        if edges:
            # A dark underlay keeps the white physical-level edges legible.
            ax.add_collection3d(Line3DCollection(edges, colors='#17202a', linewidths=1.15))
            ax.add_collection3d(Line3DCollection(edges, colors='white', linewidths=.65))
        bounds.append(vertices.reshape(-1, 3))
    # Curves are existing grid-node magnetic traces, not straight proxy lines.
    for jj in np.linspace(0, run.fy, 3, dtype=int):
        for ii in np.linspace(0, run.fx, 3, dtype=int):
            line = project(nodes[jj, ii, :, :])
            ax.plot(*line.T, color='#343434', linewidth=.8, linestyle='--', alpha=.85)
            bounds.append(line)
    limits = np.concatenate(bounds)
    low, high = limits.min(axis=0), limits.max(axis=0)
    for setter, lo, hi in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), low, high):
        pad = max(1., .035 * (hi - lo))
        setter(lo - pad, hi + pad)
    ax.set_box_aspect(np.maximum(high - low, 1))
    ax.view_init(elev=25, azim=-60)
    ax.set_xlabel('Reference east (km)', labelpad=10)
    ax.set_ylabel('Reference north (km)', labelpad=10)
    ax.set_zlabel('Reference radial-up displacement (km)', labelpad=10)
    time = float(row['time_s'])
    ax.set_title(f"{'Electron' if args.quantity == 'electron' else 'O+'} density; t = {time:g} s\n"
                 f'Curved grid-label cuts: i={i}, j={j}, k={k}; {run.levels} physical levels', pad=20)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=.55,
                 pad=.08, label=r'$n/(10^{12}\,\mathrm{m}^{-3})$')
    ax.legend(handles=[Line2D([0], [0], color='#555555', lw=1.5, label='Actual AMR level edges (white)'),
                       Line2D([0], [0], color='#343434', ls='--', label='Saved traced field lines')],
              loc='upper left', fontsize=9)
    fig.text(.03, .025, 'Cuts follow magnetic-grid labels; radial-up is not geographic altitude. '
             'Coarse O+ is reconstructed. Colour limits are display limits only.', fontsize=9)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ('.png', '.pdf'):
        fig.savefig(args.out.with_suffix(suffix), dpi=180)
    plt.close(fig)
    inputs.append(node_path)
    provenance = {
        'time_s': time, 'frame': int(row['frame']), 'quantity': args.quantity,
        'slice_indices': {'i': i, 'j': j, 'k': k}, 'surface_display_stride': args.stride,
        'AMR_outline': 'all digital edges between actual density levels; no per-cell mesh',
        'projection': 'fixed ECEF east/north/radial-up Cartesian displacement, km',
        'reference_ECEF_m': reference.tolist(), 'projection_basis_rows': basis.tolist(),
        'display_limits_normalized': [0, 1], 'density_floor_in_model_added': False,
        'input_SHA256': {str(path.relative_to(run.root)): sha256(path) for path in inputs},
        'plotter_SHA256': sha256(Path(__file__)),
        'reader_SHA256': sha256(Path(__file__).with_name('frame_reader.py')),
        'numpy_version': np.__version__, 'matplotlib_version': matplotlib.__version__,
    }
    args.out.with_suffix('.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(args.out.with_suffix('.png'), flush=True)


if __name__ == '__main__':
    main()
