# match_colours_controlled.py

Adjust the colours of a **target** aerial photograph so that it resembles a **reference** aerial photograph. The two images do **not** need to overlap, or to share dimensions, CRS or resolution.

Only colour changes. The output GeoTIFF has the target's original size, resolution, coordinate system and pixel alignment. The script changes appearance, not spatial detail or accuracy.

It is a visual adjustment, **not** radiometric calibration.

## Contents

- [Quick start](#quick-start)
- [How it works](#how-it-works)
- [Command-line options](#command-line-options)
- [Recipes](#recipes)
- [Reading the statistics table](#reading-the-statistics-table)
- [Inputs and validity rules](#inputs-and-validity-rules)
- [Output](#output)
- [Python API](#python-api)
- [Limitations](#limitations)
- [Testing notes](#testing-notes)
- [Troubleshooting](#troubleshooting)

## Quick start

Requirements: Python 3.8+, `numpy` and `rasterio`.

```bash
pip install numpy rasterio
# or, with conda (often the smoothest route on Windows):
conda install -c conda-forge numpy rasterio
```

Preview the effect first. Nothing is written, but an output path is still required by the command line:

```bash
python match_colours_controlled.py target.tif reference.tif matched.tif --dry-run
```

Then run it for real:

```bash
python match_colours_controlled.py target.tif reference.tif matched.tif
```

The output file must not already exist unless you add `--overwrite`.

## How it works

1. **Sample.** Each image is read as a bounded, nearest-neighbour sample (about 2 million pixels by default), so large photographs are handled quickly. Only valid pixels count, and a pixel counts only if every selected band is valid, so all bands are described by the same pixels.
2. **Estimate the match.** For each band, and for luma, the script works out how the target's values should change to resemble the reference. Two methods are available:
   - `tails` (default): a gain and offset per band, anchored on the dark and bright tails of each histogram (by default the 2nd and 98th percentiles). Gain is limited to between `1/--max-gain` and `--max-gain`.
   - `cdf`: full histogram (CDF) matching per band. This reproduces the reference histogram shape exactly.
3. **Split into brightness and colour.** The correction is separated into brightness (weighted RGB luma) and colour (chroma, which is RGB minus luma). Each is scaled by its own strength and bounded by its own limit:
   - brightness change = estimated change × `--brightness-strength`, capped at ±`--brightness-limit`
   - chroma change = estimated change × `--colour-strength`, scaled so no channel moves by more than `--colour-limit`
4. **Gamut guard.** If the requested brightness change would push any channel outside the value range, brightness is reduced for that pixel. The chroma limit is never exceeded as a result.
5. **Write.** The target is processed in windows of roughly one million pixels, so memory stays bounded, and written to a new GeoTIFF.

### Why `tails` is the default

Different land cover is the main way a colour match goes wrong. If the target is mostly water and the reference is mostly bare ground, full histogram matching "corrects" the difference in what the photos show as if it were a difference in colour. Land cover shows up mostly in the mid-tones, while shadows, deep water, cloud and bright surfaces tend to be more alike between scenes. Anchoring on the tails is therefore much less sensitive to land-cover proportions. It is not immune: see [Limitations](#limitations).

`--method cdf` is available if you want full histogram matching.

### What "brightness" and "colour" mean here

Brightness is weighted RGB luma in the stored encoding, not physical luminance. Chroma is RGB minus luma, not perceptual Lab colour and not a hue-angle constraint. The colour limit bounds the *chroma* change, not the total RGB change and not brightness. The limits apply before quantisation, so integer outputs can differ from them by about one level.

## Command-line options

```
python match_colours_controlled.py TARGET REFERENCE OUTPUT [options]
```

| Option | Default | Meaning |
|---|---|---|
| `--bands B [B ...]` | `1 2 3` (or `1` for a one-band target) | Three band numbers in R G B order, or one greyscale band. Required for a two-band target. |
| `--method {tails,cdf}` | `tails` | How the match is estimated (see above). |
| `--tail-percentile P` | `2` | Percentile defining the dark and bright tails for `tails`; uses P and 100 − P. Must be between 0 and 50. |
| `--max-gain G` | `2` | Largest contrast gain for `tails` (and its reciprocal). Must be at least 1; `1` allows a pure offset only. |
| `--brightness-strength S` | `0.8` | Fraction of the estimated brightness change to apply, 0 to 1. |
| `--brightness-limit L` | `0.10` | Maximum brightness change as a fraction of the value range, 0 to 1. |
| `--colour-strength S` | `0.25` | Fraction of the estimated chroma change to apply, 0 to 1. |
| `--colour-limit L` | `0.03` | Maximum change of any chroma channel as a fraction of the value range, 0 to 1. |
| `--value-range LOW HIGH` | integer dtype range | Common encoding range of both inputs. **Required for floating-point images.** |
| `--source-stats-mask FILE` | none | Single-band raster on the same pixel grid as the target. Non-zero pixels are used for the statistics. |
| `--reference-stats-mask FILE` | none | Same, for the reference. |
| `--ignore-zero` | off | Skip pixels that are zero in **all** selected bands (unlabelled black borders). |
| `--max-samples N` | `2000000` | Approximate maximum sample pixels per image. Above this, statistics are approximate. |
| `--dry-run` | off | Print the statistics table and write nothing. |
| `--no-progress` | off | Hide stage messages and the progress bar. |
| `--quiet` | off | Hide the statistics table and the progress. |
| `--overwrite` | off | Allow replacing an existing output file. |

With both strengths at `0`, valid pixels are copied unchanged.

Progress and errors go to stderr. The statistics table and the final `Saved colour-matched image: ...` line go to stdout. The exit status is `0` on success and `1` on error.

## Recipes

Preview, then run with defaults:

```bash
python match_colours_controlled.py target.tif reference.tif matched.tif --dry-run
python match_colours_controlled.py target.tif reference.tif matched.tif
```

Compare like with like by excluding cloud, water or crops that appear in only one image. Each mask must be a single-band raster on exactly the same grid as its image, with non-zero meaning "use this pixel":

```bash
python match_colours_controlled.py target.tif reference.tif matched.tif ^
    --source-stats-mask target_land_only.tif --reference-stats-mask reference_land_only.tif
```

(Use `\` instead of `^` for line continuation on Linux and macOS.)

Imagery that is not 8-bit:

```bash
# 12-bit values stored as uint16
python match_colours_controlled.py target.tif reference.tif matched.tif --value-range 0 4095

# floating point, scaled 0 to 1
python match_colours_controlled.py target.tif reference.tif matched.tif --value-range 0 1
```

A black border that is not declared as nodata:

```bash
python match_colours_controlled.py target.tif reference.tif matched.tif --ignore-zero
```

Bands stored as B G R, or a single greyscale band:

```bash
python match_colours_controlled.py target.tif reference.tif matched.tif --bands 3 2 1
python match_colours_controlled.py target.tif reference.tif matched.tif --bands 1
```

A gentler or stronger correction:

```bash
# gentler
python match_colours_controlled.py target.tif reference.tif matched.tif --brightness-strength 0.5 --brightness-limit 0.05

# stronger colour correction (use with care; see Limitations)
python match_colours_controlled.py target.tif reference.tif matched.tif --colour-strength 0.6 --colour-limit 0.06
```

Closest to full histogram matching, with no brightness cap:

```bash
python match_colours_controlled.py target.tif reference.tif matched.tif --method cdf --brightness-limit 1
```

## Reading the statistics table

Every run (and `--dry-run`) prints a summary computed from the sampled pixels, as a percentage of the value range:

```
Statistics from sampled pixels, % of value range:
                  luma p2   p50   p98    chroma R-Y   G-Y   B-Y
  target before     20.8   39.0   57.2    +0.07  -0.02  +0.02
  reference         32.6   50.8   68.9    +0.01  -0.01  +0.04
  target after      30.3   48.4   66.5    +0.03  -0.01  -0.02
```

- **luma p2 / p50 / p98** are the dark end, median and bright end of brightness.
- **chroma R-Y, G-Y, B-Y** are the average colour offsets of each channel from luma. They show colour cast: a positive R-Y means a red-leaning image.
- "target after" should move *toward* "reference", but usually not all the way, because of the strengths and limits. A large remaining gap is normal.

A sensible check is to look at how far "target before" moves to "target after". If the median brightness shifts by much more than you would expect from the lighting alone, the two photos probably differ in land cover, and a statistics mask is worth trying.

## Inputs and validity rules

- Any raster GDAL can read is accepted. The output is always a GeoTIFF.
- Target and reference must use the **same RGB encoding, numeric scale and band order**. The selected bands must have the same data type, and the same scale and offset values, in both files.
- Alpha and palette bands cannot be selected for matching. Any bands that are not selected, including alpha, are copied unchanged.
- A pixel is corrected only if **every selected band is valid**: not masked, not nodata, not an alpha-transparent pixel, and finite. Pixels that are not corrected are copied unchanged.
- With `--ignore-zero`, a pixel is also skipped if it is zero in **all** selected bands. A pixel with a zero in only one band is treated as real data.
- Zero is otherwise a valid value.
- Sampled values and source pixels must lie inside `--value-range`. For integer images the default is the full range of the data type.
- The statistics masks affect only which pixels are used to estimate the match. All valid pixels of the target are still corrected.

## Output

- Written as a tiled (512 × 512), deflate-compressed GeoTIFF with the predictor enabled, and `BIGTIFF=IF_SAFER`.
- Size, resolution, CRS, transform, band count, data type, colour interpretation, scales, offsets, band descriptions, units and tags are copied from the target.
- A JPEG/YCbCr source is written as an ordinary RGB GeoTIFF.
- An alpha band is copied unchanged and stays a true alpha band. A non-alpha internal mask is copied into the output.
- If a matched value would equal the declared nodata value, it is moved by one level, so valid pixels never turn into nodata. The footprint of valid pixels is preserved.
- The file is first written as `<output>.part` and renamed only on success, so a failure never leaves a half-written file under the real name.
- These metadata tags record how the file was made: `COLOUR_MATCH_REFERENCE`, `COLOUR_MATCH_METHOD`, `BRIGHTNESS_STRENGTH`, `BRIGHTNESS_LIMIT`, `COLOUR_STRENGTH`, `COLOUR_LIMIT`, `TAIL_PERCENTILE`, `MAX_GAIN`, `VALUE_RANGE`, `COLOUR_MATCH_BANDS`.

## Python API

```python
from match_colours_controlled import match_non_overlapping_rasters

path = match_non_overlapping_rasters(
    "target.tif", "reference.tif", "matched.tif",
    overwrite=True,
    brightness_strength=0.6,
    verbose=False,        # no statistics table (and, by default, no progress)
)
```

The keyword arguments mirror the command-line options (`bands`, `method`, `tail_percentile`, `max_gain`, `brightness_strength`, `brightness_limit`, `colour_strength`, `colour_limit`, `value_range`, `source_stats_mask`, `reference_stats_mask`, `ignore_zero`, `max_samples`, `overwrite`, `dry_run`). `progress` defaults to the value of `verbose`; pass `progress=False` to hide only the progress output.

The function returns the output path as a string, or `None` for a dry run. Invalid arguments and unusable inputs raise `ValueError`.

## Limitations

- **Statistics cannot tell lighting from land cover.** Matching the tails reduces the problem but does not remove it. If one image contains something the other lacks at the dark or bright end (cloud, deep shadow, water, bright roofs), the estimate is skewed. Use statistics masks, and check the result by eye.
- **Colour cast is the least reliable part.** In testing, chroma correction improved far less than brightness correction, which is why the colour defaults are conservative. Raise them gradually.
- **One global correction per image.** It does not fix vignetting, hot spots, seams or changes of lighting *within* a photo, and it does not blend neighbouring photos.
- **Works in the stored encoding.** Values are treated as they are stored (typically gamma-encoded), not converted to linear light or to Lab.
- **Sampled statistics.** Above `--max-samples` pixels the statistics are approximate. They are taken from a regular, nearest-neighbour decimation.
- **Different resolutions.** Fine detail such as shadows and specular highlights shows up more in sharper imagery, which can slightly shift the tails when the two photos have very different ground sampling distance.

## Testing notes

The script was tested on a real Landsat RGB scene split into two non-overlapping halves, with known colour, haze and exposure distortions added to the target, so results could be scored against the true colours. No aerial photographs were available for testing, so treat these as indicative only.

Mean absolute error in 8-bit levels (lower is better):

| Scenario | No correction | Full histogram matching (`--method cdf` style, original script) | Default (`tails`) |
|---|---|---|---|
| Real distortion, tails line up (cloud in both scenes) | 10.1 | 16.8 | 8.1 |
| Real distortion, cloud in only one scene | 10.1 | 13.0 | 9.8 |
| Target already correct, only land cover differs | 0 | 13.3 | 6.5 |

The default is clearly safer than full histogram matching, and it helps when the tails correspond. It is roughly neutral when they do not, and it still moves an already-correct target slightly. Before using it on a batch, run `--dry-run` on one pair, and compare with `--method cdf` on your own imagery.

The documented limits were checked numerically: chroma and brightness changes stayed within `--colour-limit` and `--brightness-limit` (plus about one level of rounding for 8-bit output), and the nodata footprint of the output was identical to the source.

## Troubleshooting

| Message | Cause and fix |
|---|---|
| `Floating-point images require --value-range LOW HIGH.` | Float images have no natural range. Pass the range, for example `--value-range 0 1`. |
| `Sampled input values exceed --value-range.` | One of the images has values outside the range you gave. Widen the range, or check that both images use the same scale. |
| `Source pixels exceed --value-range.` | A part of the target that the sample did not cover has values outside the range. Nothing is left behind. Widen the range. |
| `Selected band dtypes differ...` | Convert both images to the same data type first. |
| `Input band scales/offsets differ...` | Harmonise the two images' scale and offset metadata first. |
| `Output exists: ... Use --overwrite.` | Choose a new output name or add `--overwrite`. |
| `Output must be different from both input files.` | The output path points at one of the inputs. |
| `Statistics mask ... must have exactly the same pixel grid as ...` | A stats mask must match its image in width, height and georeferencing. Re-rasterise it against that image. |
| `No valid sampled pixels...` | Everything was masked, nodata or excluded. For a stats mask, check it has non-zero pixels over valid imagery. |
| `Band N is an alpha band` | You selected an alpha band with `--bands`. Select only colour bands. |
| `Specify --bands for a two-band target.` | Two-band images are ambiguous. Pass one band, or three. |
| Import error for `rasterio` | Install it (`pip install rasterio` or `conda install -c conda-forge rasterio`). |
