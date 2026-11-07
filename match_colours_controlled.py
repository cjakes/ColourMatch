#!/usr/bin/env python3
"""Match a target aerial photograph to a reference using robust, controlled statistics.

Requires: pip install numpy rasterio
Usage: python match_colours_controlled.py target.tif reference.tif matched.tif

Inputs need not overlap or have the same dimensions, CRS or resolution.
Only colours change: output geometry is copied from the target.
Default: match bands 1, 2, 3 (or band 1 for a single-band target).
Other bands, including alpha, are copied unchanged. Band order must agree.

How the match is estimated (--method):
  tails (default) : gain/offset per band from the dark and bright tails of each
                    histogram (--tail-percentile 2 uses the 2nd and 98th
                    percentiles). Shadows/deep water and bright surfaces tend
                    to be more alike between scenes than the mid-tones, where
                    land-cover proportions dominate, so this is far less
                    sensitive to different land cover than full matching.
                    Gain is bounded by --max-gain.
  cdf             : full per-band histogram (CDF) matching, as in the original
                    script. Reproduces the reference histogram shape exactly,
                    so any difference in land cover is "corrected" too.

Controls:
  --brightness-strength 0.8 : fraction of the estimated brightness change
  --brightness-limit 0.10   : maximum absolute brightness (luma) change as a
                              fraction of --value-range
  --colour-strength 0.25    : fraction of the estimated chroma change
  --colour-limit 0.03       : maximum absolute chroma-channel change as a
                              fraction of --value-range, before quantisation
  --value-range LOW HIGH    : common input encoding range (required for floats;
                              integer dtype range is the default)
  --source-stats-mask / --reference-stats-mask : optional single-band rasters on
                              the SAME pixel grid as the image; non-zero pixels
                              are used for the statistics. Use them to compare
                              like with like (e.g. exclude cloud, water or
                              crops that occur in only one image). All valid
                              pixels of the target are still corrected.
  --dry-run                 : print before/after statistics, write nothing
  --no-progress             : hide the progress messages and bar (stderr)
  --quiet                   : hide the statistics table and the progress

Progress (stage messages and a bar with percentage, elapsed time and ETA while
writing) goes to stderr, so redirecting stdout does not capture it. On a
terminal the bar redraws one line; in log files it prints a line per 10 %.

Brightness is weighted RGB luma in the stored encoding, not physical luminance.
Chroma is RGB minus luma, not perceptual Lab colour or a hue-angle constraint.
The colour limit bounds chroma change, NOT total RGB change or brightness.
A gamut guard may reduce brightness adjustment to keep the colour limit valid.
Inputs must use the same RGB encoding, numeric scale and band order.
For 12-bit values stored as uint16, specify --value-range 0 4095.
With both strengths zero, valid source values are copied unchanged.
A pixel is corrected only if every selected band is valid; with --ignore-zero,
only pixels that are zero in ALL selected bands are skipped.
This is a visual adjustment, not radiometric calibration.
"""
import argparse
import contextlib
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import ColorInterp, MaskFlags, Resampling
from rasterio.windows import Window

_LUMA = np.array([0.2126, 0.7152, 0.0722], dtype=np.float64)
_WINDOW_PIXELS = 1_000_000   # approximate pixels per processing window
_REPORT_PIXELS = 200_000     # pixels used for the before/after report


# --------------------------------------------------------------------------
# Progress reporting (dependency-free, writes to stderr)
# --------------------------------------------------------------------------
def _clock(seconds):
    seconds = int(round(seconds))
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    return f'{hours}:{minutes:02d}:{secs:02d}' if hours else f'{minutes}:{secs:02d}'


@contextlib.contextmanager
def _stage(enabled, message):
    """Print 'message ... done (1.2 s)' around a step that cannot report fractions."""
    start = time.monotonic()
    if enabled:
        print(f'{message} ...', end=' ', file=sys.stderr, flush=True)
    try:
        yield
    except BaseException:
        if enabled:
            print('failed', file=sys.stderr, flush=True)
        raise
    if enabled:
        print(f'done ({time.monotonic() - start:.1f} s)', file=sys.stderr, flush=True)


class _Progress:
    """Progress bar with percentage, elapsed time and ETA.

    On a terminal it redraws a single line. When stderr is not a terminal
    (log file, some IDE consoles) it prints one line per 10 % instead.
    """

    def __init__(self, label, total, enabled=True, width=24):
        self.label, self.total, self.enabled, self.width = label, max(1, total), enabled, width
        self.done = 0
        self.start = time.monotonic()
        self.tty = sys.stderr.isatty()
        self.last_draw = 0.0
        self.last_decile = -1
        self.drawn = 0
        self.closed = False

    def _line(self):
        fraction = min(1.0, self.done / self.total)
        elapsed = time.monotonic() - self.start
        filled = int(round(fraction * self.width))
        text = (f'{self.label} [{"#" * filled}{"-" * (self.width - filled)}] '
                f'{fraction * 100:5.1f}%  elapsed {_clock(elapsed)}')
        if 0 < fraction < 1:
            text += f'  ETA {_clock(elapsed * (1 - fraction) / fraction)}'
        return text

    def _draw(self):
        line = self._line()
        sys.stderr.write('\r' + line + ' ' * max(0, self.drawn - len(line)))
        sys.stderr.flush()
        self.drawn = len(line)

    def update(self, amount):
        self.done += amount
        if not self.enabled:
            return
        now = time.monotonic()
        if self.tty:
            if now - self.last_draw >= 0.1 or self.done >= self.total:
                self._draw()
                self.last_draw = now
        else:
            decile = int(10 * min(1.0, self.done / self.total))
            if decile > self.last_decile:
                self.last_decile = decile
                print(self._line(), file=sys.stderr, flush=True)

    def close(self):
        """Finish the line (also after an error, so messages start on a new line)."""
        if self.enabled and self.tty and not self.closed:
            self._draw()
            sys.stderr.write('\n')
            sys.stderr.flush()
        self.closed = True


# --------------------------------------------------------------------------
# Sampling
# --------------------------------------------------------------------------
def _same_grid(a, b):
    if (a.width, a.height) != (b.width, b.height):
        return False
    tol = 1e-3 * min(abs(a.transform.a), abs(a.transform.e))
    return bool(np.allclose(tuple(a.transform)[:6], tuple(b.transform)[:6],
                            rtol=0, atol=tol))


def _sample(dataset, bands, max_samples, ignore_zero, stats_mask_path):
    """Return valid sampled pixels as float64 array (len(bands), M).

    Validity matches the pixel-processing rule exactly: a pixel counts only if
    every selected band is valid, so all bands are described by the same pixels.
    """
    factor = min(1.0, math.sqrt(max_samples / (dataset.width * dataset.height)))
    height = max(1, int(dataset.height * factor))
    width = max(1, int(dataset.width * factor))
    image = dataset.read(bands, out_shape=(len(bands), height, width),
                         masked=True, resampling=Resampling.nearest)
    data = image.data.astype(np.float64)
    valid = ~np.ma.getmaskarray(image).any(axis=0) & np.isfinite(data).all(axis=0)
    if ignore_zero:
        valid &= (data != 0).any(axis=0)
    if stats_mask_path is not None:
        with rasterio.open(stats_mask_path) as mask:
            if not _same_grid(mask, dataset):
                raise ValueError(f'Statistics mask {stats_mask_path} must have exactly '
                                 f'the same pixel grid as {dataset.name}.')
            flags = mask.read(1, out_shape=(height, width), resampling=Resampling.nearest)
        valid &= np.nan_to_num(flags.astype(np.float64)) != 0
    values = data[:, valid]
    if not values.shape[1]:
        where = ' inside the statistics mask' if stats_mask_path else ''
        raise ValueError(f'No valid sampled pixels{where}: {dataset.name}')
    return values


# --------------------------------------------------------------------------
# Matching curves
# --------------------------------------------------------------------------
def _cdf_curve(source, reference):
    """Full histogram matching (empirical CDF to empirical quantiles)."""
    src_values, src_counts = np.unique(source, return_counts=True)
    ref_values, ref_counts = np.unique(reference, return_counts=True)
    src_cdf = np.cumsum(src_counts, dtype=np.float64)
    src_cdf /= src_cdf[-1]
    ref_cdf = np.cumsum(ref_counts, dtype=np.float64)
    ref_cdf /= ref_cdf[-1]
    return ('cdf', src_values.astype(np.float64), np.interp(src_cdf, ref_cdf, ref_values))


def _tail_curve(source, reference, tail, max_gain):
    """Gain/offset that aligns the dark and bright tails of the two histograms."""
    src_lo, src_hi = np.percentile(source, (tail, 100.0 - tail))
    ref_lo, ref_hi = np.percentile(reference, (tail, 100.0 - tail))
    if src_hi > src_lo:
        gain = float(np.clip((ref_hi - ref_lo) / (src_hi - src_lo), 1.0 / max_gain, max_gain))
    else:
        gain = 1.0
    offset = (ref_lo + ref_hi) / 2.0 - gain * (src_lo + src_hi) / 2.0
    return ('affine', gain, float(offset))


def _make_curves(source, reference, method, tail, max_gain):
    def build(s, r):
        return _cdf_curve(s, r) if method == 'cdf' else _tail_curve(s, r, tail, max_gain)
    band_curves = [build(s, r) for s, r in zip(source, reference)]
    luma_curve = build(_LUMA @ source, _LUMA @ reference) if len(band_curves) == 3 else None
    return band_curves, luma_curve


def _apply_curve(curve, values):
    if curve[0] == 'affine':
        return values * curve[1] + curve[2]
    return np.interp(values, curve[1], curve[2])


# --------------------------------------------------------------------------
# Pixel transform
# --------------------------------------------------------------------------
def _controlled_rgb(original, matched, luma_target, brightness_strength, colour_strength,
                    colour_limit, brightness_limit, low, high):
    """Arrays have shape (3, N). Bound chroma and brightness change vs original.

    Shared-vector scaling avoids independently clipping RGB corrections.
    The final luma is restricted to the interval compatible with the new
    chroma and encoding range. This preserves the chroma-change bound.
    """
    span = high - low
    x = (original - low) / span
    target = (matched - low) / span
    y0 = _LUMA @ x
    yt = (luma_target - low) / span
    c0 = x - y0
    ct = target - _LUMA @ target
    dc = colour_strength * (ct - c0)
    peak = np.max(np.abs(dc), axis=0)
    dc *= np.minimum(1.0, colour_limit / np.maximum(peak, 1e-15))
    # Ensure a legal luma exists: chroma channel spread must not exceed one.
    # Find the largest common correction fraction satisfying all channel pairs.
    gamut = np.ones(x.shape[1], dtype=np.float64)
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            growth = dc[i] - dc[j]
            allowed = (1.0 - (c0[i] - c0[j])) / np.maximum(growth, 1e-15)
            gamut = np.minimum(gamut, np.where(growth > 0, allowed, 1.0))
    c = c0 + dc * np.clip(gamut, 0.0, 1.0)
    y = y0 + np.clip(brightness_strength * (yt - y0), -brightness_limit, brightness_limit)
    y = np.clip(y, -np.min(c, axis=0), 1.0 - np.max(c, axis=0))
    return np.clip((y + c) * span + low, low, high)


def _adjust(original, curves, brightness_strength, colour_strength, colour_limit,
            brightness_limit, low, high):
    band_curves, luma_curve = curves
    matched = np.stack([_apply_curve(c, original[i]) for i, c in enumerate(band_curves)])
    if original.shape[0] == 3:
        luma_target = _apply_curve(luma_curve, _LUMA @ original)
        return _controlled_rgb(original, matched, luma_target, brightness_strength,
                               colour_strength, colour_limit, brightness_limit, low, high)
    limit = brightness_limit * (high - low)
    shift = np.clip(brightness_strength * (matched - original), -limit, limit)
    return np.clip(original + shift, low, high)


def _avoid_nodata(values, nodata, dtype, high):
    """Nudge matched values that would collide with the declared nodata value."""
    if nodata is None or not np.isfinite(nodata):
        return values
    hit = values == nodata
    if hit.any():
        if dtype.kind == 'f':
            values[hit] = np.nextafter(dtype.type(nodata), dtype.type(np.inf))
        else:
            values[hit] = nodata + 1 if nodata + 1 <= high else nodata - 1
    return values


# --------------------------------------------------------------------------
# Reporting and windows
# --------------------------------------------------------------------------
def _describe(values, low, high):
    x = (values - low) / (high - low)
    y = _LUMA @ x if x.shape[0] == 3 else x[0]
    p2, p50, p98 = np.percentile(y, (2, 50, 98)) * 100
    chroma = (x - y).mean(axis=1) * 100 if x.shape[0] == 3 else None
    text = f'{p2:6.1f} {p50:6.1f} {p98:6.1f}'
    if chroma is not None:
        text += '   ' + ' '.join(f'{v:+6.2f}' for v in chroma)
    return text


def _print_report(source, reference, adjusted, low, high):
    header = '                  luma p2   p50   p98'
    if source.shape[0] == 3:
        header += '    chroma R-Y   G-Y   B-Y'
    print('Statistics from sampled pixels, % of value range:')
    print(header)
    for name, values in (('target before', source), ('reference', reference),
                         ('target after', adjusted)):
        print(f'  {name:14s}  {_describe(values, low, high)}')


def _windows(dataset):
    """Windows of roughly _WINDOW_PIXELS, aligned to the file's block layout."""
    block_h, block_w = dataset.block_shapes[0]
    if block_w < dataset.width:
        step_x = block_w * max(1, 1024 // block_w)
        step_y = block_h * max(1, 1024 // block_h)
    else:
        step_x = dataset.width
        step_y = block_h * max(1, _WINDOW_PIXELS // (dataset.width * block_h))
    for row in range(0, dataset.height, step_y):
        for col in range(0, dataset.width, step_x):
            yield Window(col, row, min(step_x, dataset.width - col),
                         min(step_y, dataset.height - row))


# --------------------------------------------------------------------------
# Main routine
# --------------------------------------------------------------------------
def match_non_overlapping_rasters(source_path, reference_path, output_path,
                                  bands=None, max_samples=2_000_000,
                                  ignore_zero=False, overwrite=False,
                                  brightness_strength=0.8, colour_strength=0.25,
                                  colour_limit=0.03, brightness_limit=0.10,
                                  value_range=None, method='tails',
                                  tail_percentile=2.0, max_gain=2.0,
                                  source_stats_mask=None, reference_stats_mask=None,
                                  dry_run=False, verbose=True, progress=None):
    """Colour-match source to reference, preserving source spatial metadata.

    Raster masks/nodata determine validity. Zero is valid unless ignore_zero
    is explicitly enabled. Statistics are approximate above max_samples pixels.
    Each selected source band maps to the same numbered reference band.
    The output is written to a temporary file and renamed only on success.
    verbose controls the statistics table; progress (default: same as verbose)
    controls the stage messages and progress bar on stderr.
    """
    show_progress = verbose if progress is None else bool(progress)
    stages = 3 if dry_run else 4
    for name, value in [('brightness_strength', brightness_strength),
                        ('colour_strength', colour_strength), ('colour_limit', colour_limit),
                        ('brightness_limit', brightness_limit)]:
        if not np.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f'{name} must be finite and between 0 and 1.')
    if method not in ('tails', 'cdf'):
        raise ValueError("method must be 'tails' or 'cdf'.")
    if not np.isfinite(tail_percentile) or not 0 < tail_percentile < 50:
        raise ValueError('tail_percentile must be greater than 0 and less than 50.')
    if not np.isfinite(max_gain) or max_gain < 1:
        raise ValueError('max_gain must be at least 1.')
    if max_samples < 1:
        raise ValueError('max_samples must be positive.')
    source_path, reference_path, output_path = map(
        Path, (source_path, reference_path, output_path))
    if not dry_run:
        if output_path.resolve() in (source_path.resolve(), reference_path.resolve()):
            raise ValueError('Output must be different from both input files.')
        if output_path.exists() and not overwrite:
            raise FileExistsError(f'Output exists: {output_path}. Use --overwrite.')
        if not output_path.parent.is_dir():
            raise ValueError(f'Output directory does not exist: {output_path.parent}')

    with rasterio.open(source_path) as src, rasterio.open(reference_path) as ref:
        if bands is None:
            if src.count == 1:
                bands = [1]
            elif src.count >= 3:
                bands = [1, 2, 3]
            else:
                raise ValueError('Specify --bands for a two-band target.')
        bands = list(bands)
        if not bands or len(set(bands)) != len(bands):
            raise ValueError('Provide a non-empty list of unique band numbers.')
        if len(bands) not in (1, 3):
            raise ValueError('Select one greyscale band or three bands in R G B order.')
        for band in bands:
            if band < 1 or band > min(src.count, ref.count):
                raise ValueError(f'Band {band} must exist in both inputs.')
            for ds in (src, ref):
                if ds.colorinterp[band - 1] == ColorInterp.alpha:
                    raise ValueError(f'Band {band} is an alpha band; do not colour-match it.')
                if ds.colorinterp[band - 1] == ColorInterp.palette:
                    raise ValueError('Convert palette imagery to RGB before colour matching.')
                if np.dtype(ds.dtypes[band - 1]).kind not in 'uif':
                    raise ValueError('Only real-valued integer or floating-point bands are supported.')
            if (src.dtypes[band - 1] != ref.dtypes[band - 1] or
                    src.dtypes[band - 1] != src.dtypes[bands[0] - 1]):
                raise ValueError('Selected band dtypes differ. Convert both images to the same encoding first.')
            if (src.scales[band - 1] != ref.scales[band - 1] or
                    src.offsets[band - 1] != ref.offsets[band - 1]):
                raise ValueError('Input band scales/offsets differ. Harmonise inputs first.')
        dtype = np.dtype(src.dtypes[bands[0] - 1])
        if value_range is None:
            if dtype.kind not in 'ui':
                raise ValueError('Floating-point images require --value-range LOW HIGH.')
            info = np.iinfo(dtype)
            low, high = float(info.min), float(info.max)
        else:
            if len(value_range) != 2:
                raise ValueError('value_range requires LOW HIGH.')
            low, high = map(float, value_range)
        if not np.isfinite([low, high]).all() or high <= low:
            raise ValueError('Value range must be finite with HIGH greater than LOW.')
        if dtype.kind in 'ui':
            info = np.iinfo(dtype)
            if low < info.min or high > info.max:
                raise ValueError('Value range exceeds target dtype limits.')

        with _stage(show_progress, f'[1/{stages}] Sampling target ({src.width} x {src.height})'):
            source_sample = _sample(src, bands, max_samples, ignore_zero, source_stats_mask)
        with _stage(show_progress, f'[2/{stages}] Sampling reference ({ref.width} x {ref.height})'):
            reference_sample = _sample(ref, bands, max_samples, ignore_zero, reference_stats_mask)
        for values in (source_sample, reference_sample):
            if values.min() < low or values.max() > high:
                raise ValueError('Sampled input values exceed --value-range. Check both input scales.')
        with _stage(show_progress, f'[3/{stages}] Estimating colour match ({method})'):
            curves = _make_curves(source_sample, reference_sample, method,
                                  tail_percentile, max_gain)
        settings = (brightness_strength, colour_strength, colour_limit,
                    brightness_limit, low, high)

        if verbose or dry_run:
            step = max(1, source_sample.shape[1] // _REPORT_PIXELS)
            preview = source_sample[:, ::step]
            _print_report(preview, reference_sample[:, ::max(1, reference_sample.shape[1] // _REPORT_PIXELS)],
                          _adjust(preview, curves, *settings) if
                          (brightness_strength > 0 or colour_strength > 0) else preview,
                          low, high)
        if dry_run:
            return None

        flags = src.mask_flag_enums
        has_alpha = any(MaskFlags.alpha in f for f in flags)
        copy_mask = not has_alpha and any(MaskFlags.per_dataset in f for f in flags)
        nodata_values = [src.nodatavals[b - 1] for b in bands]

        profile = src.profile.copy()
        profile.update(driver='GTiff', compress='deflate', tiled=True,
                       blockxsize=512, blockysize=512, BIGTIFF='IF_SAFER',
                       NUM_THREADS='ALL_CPUS',
                       predictor=3 if dtype.kind == 'f' else 2)
        if str(profile.get('photometric', '')).upper() == 'YCBCR':
            profile['photometric'] = 'RGB'   # YCbCr is only valid with JPEG compression
        temp_path = output_path.with_name(output_path.name + '.part')
        bar = _Progress(f'[4/{stages}] Writing', src.width * src.height, show_progress)
        try:
            # An internal mask keeps a per-dataset mask self-contained in the GeoTIFF.
            with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
                with rasterio.open(temp_path, 'w', **profile) as dst:
                    dst.colorinterp = src.colorinterp
                    dst.scales = src.scales
                    dst.offsets = src.offsets
                    dst.update_tags(**src.tags())
                    for band in range(1, src.count + 1):
                        dst.update_tags(band, **src.tags(band))
                        if src.descriptions[band - 1]:
                            dst.set_band_description(band, src.descriptions[band - 1])
                        if src.units[band - 1]:
                            dst.set_band_unit(band, src.units[band - 1])
                    indices = np.array(bands) - 1
                    for window in _windows(src):
                        data = src.read(window=window)
                        masks = src.read_masks(bands, window=window)
                        rgb = data[indices].astype(np.float64)
                        valid = np.all(masks != 0, axis=0) & np.all(np.isfinite(rgb), axis=0)
                        if ignore_zero:
                            valid &= np.any(rgb != 0, axis=0)
                        if np.any(valid) and (brightness_strength > 0 or colour_strength > 0):
                            original = rgb[:, valid]
                            if np.any((original < low) | (original > high)):
                                raise ValueError('Source pixels exceed --value-range.')
                            adjusted = _adjust(original, curves, *settings)
                            if dtype.kind in 'ui':
                                adjusted = np.rint(adjusted)
                            for i, band in enumerate(bands):
                                if not copy_mask:
                                    adjusted[i] = _avoid_nodata(adjusted[i], nodata_values[i],
                                                                dtype, high)
                                data[band - 1][valid] = adjusted[i].astype(dtype)
                        dst.write(data, window=window)
                        if copy_mask:
                            dst.write_mask(src.dataset_mask(window=window), window=window)
                        bar.update(int(window.width) * int(window.height))
                    dst.update_tags(COLOUR_MATCH_REFERENCE=reference_path.name,
                                    COLOUR_MATCH_METHOD=f'{method} with luma/chroma controls',
                                    BRIGHTNESS_STRENGTH=str(brightness_strength),
                                    BRIGHTNESS_LIMIT=str(brightness_limit),
                                    COLOUR_STRENGTH=str(colour_strength),
                                    COLOUR_LIMIT=str(colour_limit),
                                    TAIL_PERCENTILE=str(tail_percentile),
                                    MAX_GAIN=str(max_gain),
                                    VALUE_RANGE=f'{low},{high}',
                                    COLOUR_MATCH_BANDS=','.join(map(str, bands)))
            os.replace(temp_path, output_path)
        except BaseException:
            bar.close()
            if temp_path.exists():
                temp_path.unlink()
            raise
        bar.close()
    return str(output_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('source', help='Target image whose colours will change')
    parser.add_argument('reference', help='Reference image supplying the desired colours')
    parser.add_argument('output', help='New output GeoTIFF')
    parser.add_argument('--bands', type=int, nargs='+', help='Three band numbers in R G B order, or one greyscale band')
    parser.add_argument('--method', choices=('tails', 'cdf'), default='tails',
                        help='tails: robust tail-anchored gain/offset (default); cdf: full histogram matching')
    parser.add_argument('--tail-percentile', type=float, default=2.0,
                        help='Percentile defining the dark/bright tails for --method tails (default: 2)')
    parser.add_argument('--max-gain', type=float, default=2.0,
                        help='Largest contrast gain (and its reciprocal) for --method tails (default: 2)')
    parser.add_argument('--max-samples', type=int, default=2_000_000,
                        help='Maximum approximate sample pixels per image (default: 2000000)')
    parser.add_argument('--ignore-zero', action='store_true',
                        help='Skip pixels that are zero in all selected bands (unlabelled black borders)')
    parser.add_argument('--brightness-strength', type=float, default=0.8,
                        help='Brightness correction strength, 0 to 1 (default: 0.8)')
    parser.add_argument('--brightness-limit', type=float, default=0.10,
                        help='Maximum brightness change / value range (default: 0.10)')
    parser.add_argument('--colour-strength', type=float, default=0.25,
                        help='Chroma correction strength, 0 to 1 (default: 0.25)')
    parser.add_argument('--colour-limit', type=float, default=0.03,
                        help='Maximum chroma-channel change / value range (default: 0.03)')
    parser.add_argument('--value-range', type=float, nargs=2, metavar=('LOW', 'HIGH'),
                        help='Common encoding range; required for float inputs')
    parser.add_argument('--source-stats-mask', help='Raster (same grid as target); non-zero pixels used for statistics')
    parser.add_argument('--reference-stats-mask', help='Raster (same grid as reference); non-zero pixels used for statistics')
    parser.add_argument('--dry-run', action='store_true', help='Print before/after statistics only; write nothing')
    parser.add_argument('--no-progress', action='store_true', help='Do not show progress messages or the progress bar')
    parser.add_argument('--quiet', action='store_true', help='Do not print the statistics summary or progress')
    parser.add_argument('--overwrite', action='store_true', help='Allow replacing an existing output')
    args = parser.parse_args()
    try:
        result = match_non_overlapping_rasters(
            args.source, args.reference, args.output, bands=args.bands,
            max_samples=args.max_samples, ignore_zero=args.ignore_zero,
            overwrite=args.overwrite,
            brightness_strength=args.brightness_strength, colour_strength=args.colour_strength,
            colour_limit=args.colour_limit, brightness_limit=args.brightness_limit,
            value_range=args.value_range, method=args.method,
            tail_percentile=args.tail_percentile, max_gain=args.max_gain,
            source_stats_mask=args.source_stats_mask,
            reference_stats_mask=args.reference_stats_mask,
            dry_run=args.dry_run, verbose=not args.quiet,
            progress=False if args.no_progress else None)
    except (ValueError, OSError, rasterio.errors.RasterioError) as exc:
        parser.exit(1, f'Error: {exc}\n')
    if result:
        print(f'Saved colour-matched image: {result}')


if __name__ == '__main__':
    main()
