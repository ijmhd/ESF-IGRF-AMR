#!/usr/bin/env python3
"""Read complete AMR visualization frames; never start or modify a model.

Run array readers on the compute system holding the simulation. Files use
little-endian float64 in [j, i, k] order; stored densities are n / 1e12 m^-3.
The canonical display grid contains reconstructed coarse O+, not independent
fine-grid degrees of freedom. Actual density refinement is read separately.
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path):
    """Stream an input checksum without loading its full payload into memory."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


class SavedRun:
    """Bind dimensions, frame times and mesh epochs to one saved run directory."""

    def __init__(self, root):
        """Read text metadata only; array access happens after frame validation."""
        self.root = Path(root).resolve()
        self.meta = json.loads((self.root / 'metadata.json').read_text())
        self.fx, self.fy, self.ns = (int(self.meta[key]) for key in
                                    ('fine_nx', 'fine_ny', 'ns'))
        self.shape = (self.fy, self.fx, self.ns)
        self.count = self.fx * self.fy * self.ns
        self.night = bool(self.meta['night_chemistry'])
        self.levels = int(self.meta['amr_physical_levels'])
        if self.levels not in (2, 3) or min(self.shape) < 1:
            raise ValueError('Expected a two- or three-level AMR run')
        self.structures = {}
        structure = self.root / 'structure_frames.csv'
        if structure.is_file():
            with structure.open(newline='') as handle:
                for row in csv.DictReader(handle):
                    try:
                        self.structures[int(row['frame'])] = row
                    except (ValueError, TypeError, KeyError):
                        continue  # An unfinished trailing CSV row is not a frame.

    def checked(self, relative, size):
        """Return an exact-sized file, refusing absent or partly written arrays."""
        path = self.root / relative
        if not path.is_file() or path.stat().st_size != size:
            raise ValueError(f'Missing/incomplete {path}; expected {size} bytes')
        return path

    def frame_files(self, row, quantity='electron'):
        """Bind required density and actual-level files to this frame's epoch."""
        tag = f"{int(row['frame']):04d}"
        files = {'density': self.checked(f'display_density_{tag}.bin', 8 * self.count)}
        if self.night and quantity == 'electron':
            files['molecular'] = self.checked(f'molecular_{tag}.bin', 16 * self.count)
        elif quantity == 'electron':
            files['geometry'] = self.checked('geometry_cells.bin', 48 * self.count)
        if self.levels == 3:
            files['level'] = self.checked(f'amr_level_{tag}.bin', self.count)
        else:
            structure = self.structures[int(row['frame'])]
            size = int(self.meta['root_ny']) * int(self.meta['root_nx']) * self.ns
            files['mask'] = self.checked(structure['density_mask'], size)
        return files

    def complete_frames(self, quantity='electron'):
        """Snapshot saved times with complete inputs; do not invent/interpolate time."""
        records = []
        with (self.root / 'frames.csv').open(newline='') as handle:
            for row in csv.DictReader(handle):
                try:
                    time = float(row['time_s'])
                    if not np.isfinite(time):
                        continue
                    self.frame_files(row, quantity)
                except (ValueError, KeyError, TypeError):
                    continue
                records.append(row)
        records.sort(key=lambda row: (float(row['time_s']), int(row['frame'])))
        return records

    def select(self, time=None, quantity='electron'):
        """Choose the latest complete frame, or the exact requested saved time."""
        records = self.complete_frames(quantity)
        if time is not None:
            records = [row for row in records if abs(float(row['time_s']) - time) < 1e-8]
        if not records:
            raise ValueError('No complete saved frame at the requested time')
        return records[-1]

    def read(self, row, quantity='electron'):
        """Return normalized density, physical levels and exact input provenance.

        Ne is O+ + NO+ + O2+ in night mode. Otherwise the prescribed physical
        E-layer Chapman population is added using the same profile as C/plot_slice.
        No density floor, smoothing or solution modification is applied.
        """
        files = self.frame_files(row, quantity)
        density = np.fromfile(files['density'], '<f8').reshape(self.shape)
        if 'molecular' in files:
            molecular = np.memmap(files['molecular'], '<f8', mode='r',
                                  shape=(2, *self.shape))
            density += molecular[0] + molecular[1]
        elif 'geometry' in files:
            geometry = np.memmap(files['geometry'], '<f8', mode='r',
                                 shape=(*self.shape, 6))
            z = (geometry[..., 3] - 110000.) / 12000.
            density += self.meta['E_peak_m3'] / 1e12 * np.exp(1 - z - np.exp(-z))
        if 'level' in files:
            level = np.fromfile(files['level'], 'u1').reshape(self.shape)
        else:
            mask = np.fromfile(files['mask'], 'u1').reshape(
                int(self.meta['root_ny']), int(self.meta['root_nx']), self.ns)
            level = mask.repeat(2, axis=0).repeat(2, axis=1)
        if (level.shape != density.shape or np.any(level >= self.levels)
                or not np.isfinite(density).all() or np.any(density < 0)):
            raise ValueError('Invalid saved density or physical AMR labels')
        inputs = [self.root / 'metadata.json', self.root / 'frames.csv', *files.values()]
        if 'mask' in files:
            inputs.append(self.root / 'structure_frames.csv')
        return density, level, inputs

    def nodes(self):
        """Memory-map actual traced ECEF vertex positions in metres, read-only."""
        shape = (self.fy + 1, self.fx + 1, self.ns + 1, 3)
        path = self.checked('nodes.bin', 24 * (self.fy + 1) * (self.fx + 1) * (self.ns + 1))
        return np.memmap(path, '<f8', mode='r', shape=shape), path
