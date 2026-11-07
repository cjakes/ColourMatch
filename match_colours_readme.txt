The script adjusts a target aerial photo’s colours to resemble a reference aerial photo, using histogram matching for each colour band.
It:
compares the distribution of pixel values in each image and uses that comparison to adjust the target’s red, green and blue bands
excludes masked or nodata pixels, while preserving alpha and other unmatched bands
saves a new GeoTIFF with the target’s original dimensions, resolution, coordinate system and pixel alignment
processes the image in blocks to reduce memory use. The images do not need to overlap.
1.py
It changes appearance, not spatial detail or accuracy. Different land-cover proportions can influence the colour match.

python match_colours.py input_image.tif reference_image.tif corrected_image.tif

reference_image.tif = selected colour balanced image for that local area/vegetation type
input_image.tif = image with poor colour balance
corrected_image.tif = fixed image

python match_colours_controlled.py bad.tif ok.tif fixed.tif --overwrite

python match_colours_controlled.py bad.tif ok.tif fixed2.tif --brightness-strength 0.8 --colour-strength 0.25 --colour-limit 0.03 --overwrite

python match_colours_controlled.py bad.tif ok.tif fixed3.tif  --brightness-strength 0.8 --colour-strength 0.1 --colour-limit 0.015 --overwrite



Claude Notes
I’ve added progress reporting to match_colours_controlled.py, and the output pixels are identical to the previous version. I tested it on a terminal, in a log file, with each flag, and with an error mid-write, and the earlier regression checks still pass.

A run now looks like this:

[1/4] Sampling target (5000 x 5000) ... done (3.5 s)
[2/4] Sampling reference (5000 x 5000) ... done (1.3 s)
[3/4] Estimating colour match (tails) ... done (0.1 s)
[4/4] Writing [##############----------]  60.0%  elapsed 0:03  ETA 0:02
Terminal: The writing bar redraws a single line with percentage, elapsed time and ETA. The ETA is based on pixels written.
Log file or IDE console: It prints one line per 10% instead, so a redirected log doesn’t fill with redraw characters.
Where it goes: Progress is written to stderr. The statistics table and the final “Saved…” line stay on stdout, so > log.txt captures only those.
Windows: The bar uses plain # and -, so it should display in Anaconda Prompt and PowerShell, though I could only test it on Linux.
Sampling steps: Sampling a large photo is a single read, so those stages show a message and a timing rather than a percentage.
Errors: If something fails mid-write, the bar line is finished first, so the error message starts on a clean line and the half-written file is still removed.

There are two new flags:

--no-progress: Hides the progress output but keeps the statistics table.
--quiet: Now hides both the table and the progress. Before this change it only hid the table.

--dry-run shows the three sampling and estimating stages and then the table, with no writing bar. If you call match_non_overlapping_rasters from your own code, progress= defaults to following verbose, so scripted batch runs that pass verbose=False stay silent.