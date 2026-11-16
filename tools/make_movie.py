#!/usr/bin/env python3
"""Render complete saved frames with public plotters, then optionally encode MP4.

The frame list is snapshotted at start. Titles use actual saved model times;
there is no interpolation or last-frame hold while the displayed clock advances.
At fixed video FPS, irregular saved time intervals are NOT uniform physical speed.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

from frame_reader import SavedRun, sha256


def main():
    """Preserve individual plots, exact frame selection and checked encoder logs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--view', choices=['slice', '3d'], default='slice')
    parser.add_argument('--quantity', choices=['electron', 'Oplus'], default='electron')
    parser.add_argument('--time-min', type=float)
    parser.add_argument('--time-max', type=float)
    parser.add_argument('--every', type=int, default=1, help='Keep every Nth complete saved frame')
    parser.add_argument('--field-index', type=int, help='slice k; also sets 3d --k')
    parser.add_argument('--stride', type=int, default=1, help='3d colour display sampling')
    parser.add_argument('--fps', type=int, default=8)
    parser.add_argument('--ffmpeg', default='ffmpeg', help='Encoder executable or absolute path')
    parser.add_argument('--frames-only', action='store_true', help='Do not require/run FFmpeg')
    parser.add_argument('--out', type=Path, required=True, help='NEW output directory')
    args = parser.parse_args()
    if min(args.every, args.stride, args.fps) < 1:
        parser.error('every/stride/fps must be positive')
    if args.out.exists():
        parser.error('Choose a new output directory; existing results are preserved')
    encoder = shutil.which(args.ffmpeg)
    if not args.frames_only and encoder is None:
        parser.error('FFmpeg not found; install it, supply --ffmpeg, or use --frames-only')
    run = SavedRun(args.run)
    frames = [row for row in run.complete_frames(args.quantity)
              if (args.time_min is None or float(row['time_s']) >= args.time_min)
              and (args.time_max is None or float(row['time_s']) <= args.time_max)][::args.every]
    if not frames:
        parser.error('No complete saved frames in the requested interval')
    args.out.mkdir(parents=True)
    plotter = Path(__file__).with_name('plot_3d.py' if args.view == '3d' else 'plot_slice.py')
    manifest = {'run': str(run.root), 'view': args.view, 'quantity': args.quantity,
                'fps': args.fps, 'time_policy': 'one frame per selected actual saved time; no interpolation',
                'plotter_SHA256': sha256(plotter), 'movie_tool_SHA256': sha256(Path(__file__)),
                'reader_SHA256': sha256(Path(__file__).with_name('frame_reader.py')),
                'selected_frames': [{'frame': int(row['frame']), 'time_s': float(row['time_s'])}
                                    for row in frames], 'rendered': 0, 'status': 'rendering'}
    manifest_path = args.out / 'movie_manifest.json'

    def record():
        """Persist progress before/after each subprocess so failures remain inspectable."""
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')

    record()
    for index, row in enumerate(frames):
        command = [sys.executable, str(plotter), str(run.root), '--time', row['time_s'],
                   '--quantity', args.quantity, '--out', str(args.out / f'frame_{index:05d}')]
        if args.field_index is not None:
            command += ['--k' if args.view == '3d' else '--field-index', str(args.field_index)]
        if args.view == '3d':
            command += ['--stride', str(args.stride)]
        with (args.out / f'frame_{index:05d}.log').open('w') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            manifest.update(status='plot_failed', failed_command=command,
                            failed_returncode=result.returncode)
            record()
            raise RuntimeError(f'Plot failed; inspect frame_{index:05d}.log; saved science is unchanged')
        manifest['rendered'] += 1
        record()
    if args.frames_only:
        manifest['status'] = 'frames_complete'
    else:
        movie = args.out / 'movie.mp4'
        command = [encoder, '-hide_banner', '-n', '-framerate', str(args.fps),
                   '-i', str(args.out / 'frame_%05d.png'), '-vf',
                   'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-c:v', 'libx264', '-crf', '18',
                   '-pix_fmt', 'yuv420p', str(movie)]
        manifest['encoder_command'] = command
        record()
        with (args.out / 'ffmpeg.log').open('w') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            manifest.update(status='encode_failed', failed_returncode=result.returncode)
            record()
            raise RuntimeError('Encoder failed; keep the PNG/PDF frames and inspect ffmpeg.log')
        manifest.update(status='complete', movie_SHA256=sha256(movie))
    record()
    print(manifest_path, flush=True)


if __name__ == '__main__':
    main()
