# Fluorescence lifetime imaging (FLIM)

`Lifetime/` contains the MATLAB script and the Fiji plugin of the fluorescence lifetime imaging
(FLIM) workflow for *Deinococcus radiodurans*. The workflow quantifies the spatial distribution of
the fluorescence lifetime and removes the background signal. The example used throughout this
document is a DR WT sample labelled for 2 h with **50 µM TMRDA**.

`demo/Lifetime/` is a **self-contained package** of that module: the module script, the Fiji
plugin, the entry point `runDemoLifetime.m` with its two step wrappers and the section driver, the
example data set in `demo_data/` and this Readme all live in this single folder, so the folder can
be downloaded and run on its own.

The workflow is documented against the three items of the published example data set, referred to
below as **example data 1** … **example data 3**. This package ships all three items in
`demo_data/`, under the English name of the original item; the image and metadata files themselves
are unchanged.

The module consists of four parts:

| Part | Content | Can be run unattended |
|---|---|---|
| 1 | [Raw data preparation](#1-raw-data-preparation) | No — LAS X FLIM/FCS export in the acquisition software |
| 2 | [Image stacking](#2-image-stacking) | **Yes** — covered by `runDemoLifetime.m` |
| 3 | [Background masking](#3-background-masking) | No — Fiji plugin |
| 4 | [Lifetime and intensity quantification](#4-lifetime-and-intensity-quantification) | **Yes** — covered by `runDemoLifetime.m` |

Parts 2 and 4 are the two scriptable steps of the workflow, and they are the ones the automated
demo replays (see [Automated demo](#automated-demo)).

---

## 1. Raw data preparation

The intensity and the lifetime images are exported from **LAS X FLIM/FCS** (Leica).

* Export ranges: **intensity** 0–65535, **lifetime** 0–10 (here: nanoseconds).
* Both images are saved for every field of view (FOV): the intensity image is channel 0
  (`_ch0.ome.tif`) and the lifetime image is channel 1 (`_ch1.ome.tif`).

→ produced: **example data 1** — 12 OME-TIFF images of 1024 × 1024 pixel (3 FOVs × 2 cells), plus
the LAS X `MetaData/` folder of every FOV.

## 2. Image stacking

All exported images are stacked into one intensity stack and one lifetime stack.

* MATLAB script: `lifeTcalculationbernsen.m`, **Section 1**
* **Input**: example data 1
* **Output**: example data 2
* This step runs unattended and is covered by the demo.

The section loops over the field-of-view folders and over the cells within them, reads
`<image>_ch0.ome.tif` (intensity) and `<image>_ch1.ome.tif` (lifetime) of every cell and appends
them to the intensity stack and to the lifetime stack — one page per image, in the order
field of view, then cell.

The shipped section is written for one acquisition series: its root folder, the number of images
per field of view and the field-of-view folders it loops over are hard-coded, and the file names
carry the prefix of that series. The demo wrapper replaces exactly those three settings and nothing
else (see [How the demo drives the module script](#how-the-demo-drives-the-module-script)).

## 3. Background masking

The cytoplasmic and the peripheral background signal is removed with Bernsen adaptive local
thresholding in Fiji.

* Fiji plugin: `bersenThtest` (`bersenThtest.ijm`)
* **Input**: the intensity stack (example data 2), converted to 8-bit
* **Output**: the binary mask (example data 3)
* This step needs Fiji/ImageJ and cannot be replayed by the demo.

The macro applies the *Auto Local Threshold* → *Bernsen* method to the active image; the resulting
binary mask marks the cell-wall signal (255) against the background (0) and is the mask that
step 4 applies to both stacks. The parameters recorded for the published mask differ between the
two files that describe this step, see [Note on the Bernsen parameters](#note-on-the-bernsen-parameters).

## 4. Lifetime and intensity quantification

The lifetime and the intensity of every background-free pixel are computed, together with the
lifetime distribution.

* MATLAB script: `lifeTcalculationbernsen.m`, **Section 3**
* **Input**: the stacks of example data 2 and the binary mask of example data 3
* **Output**: the background-free stacks, the lifetime histogram and the per-page samples
* This step runs unattended and is covered by the demo.

For every page of the stacks the section sets every background pixel (mask = 0) of the intensity
and of the lifetime image to zero, rescales the 16-bit lifetime grey values to nanoseconds
(`lifetime = grey / 65535 × 10 ns`), collects the pixels with a non-zero lifetime and stores them
per page. It then builds the lifetime histogram on the fixed grid `0.30 : 0.03 : 3.60 ns`
(probability normalisation), computes the mean and the standard deviation of the lifetime and of
the intensity, saves everything to a `.mat` file and draws the histogram.

The shipped section asks for the output `.mat` file with `uiputfile`, i.e. it waits for a click on
*Save*; the demo wrapper writes it to a fixed name instead (see
[How the demo drives the module script](#how-the-demo-drives-the-module-script)).

---

## Example data

The published example data set of the FLIM workflow consists of three items; this package ships
all three (35 files, ≈ 52 MB).

| Supplement item | Folder | Content | Used by the demo |
|---|---|---|---|
| 1. LAS X export data | `01_lasx_export/` | For every field of view (`0`, `1`, `2`) and every cell (1, 2) one intensity image (`_ch0.ome.tif`) and one lifetime image (`_ch1.ome.tif`), 1024 × 1024 uint16, plus the LAS X `MetaData/` folder. | step 1 (input) |
| 2. Stacked image data | `02_stacked/` | `WT50uMintensityStack.tif` and `WT50uMlifetimeStack.tif`, 6 pages of 1024 × 1024 uint16 — the output of step 2. | step 1 (reference), step 2 (input) |
| 3. Background-removed data | `03_background_masked/` | `WT50uMintensityStack-binary.tif` (the Bernsen mask, 8-bit, 0/255), `WT50uMintensityStack-filter.tif` and `WT50uMlifetimeStack-filter.tif` (the background-free stacks) — the output of step 3. | step 2 (input and reference) |

Two notes on this data set:

* The stack of example data 2 has one page per exported image, in the order field of view, then
  cell — the same order in which step 1 appends the pages, which is why the demo can compare its
  own stacks with the published ones page by page.
* `02_stacked/` and `03_background_masked/` serve a second purpose in the demo: they are the
  **reference** the wrappers compare their own output with, so the demo shows that it reproduces
  the published stacks and the published background-free images.

---

## Requirements

* MATLAB R2020a or later — **base MATLAB only** (`imread`, `imwrite`, `imquantize`, `histcounts`,
  `bar`); no toolbox is needed for either demo step.
* Fiji (ImageJ) with the *Auto Local Threshold* plugin, for the Bernsen masking (part 3).
* LAS X FLIM/FCS, for the data acquisition and export (part 1).

## How to use

1. Export the intensity and the lifetime images from LAS X FLIM/FCS (intensity 0–65535, lifetime
   0–10 ns) for every field of view.
2. Run Section 1 of `lifeTcalculationbernsen.m` to stack the exported images.
3. Convert the intensity stack to 8-bit and run `bersenThtest` in Fiji to mask the background.
4. Run Section 3 of `lifeTcalculationbernsen.m` to quantify the lifetime and the intensity of every
   background-free pixel.

## Automated demo

Parts 2 and 4 — the two steps that need no human interaction and no Fiji — are covered by the demo
that ships with this folder:

    >> cd demo/Lifetime
    >> runDemoLifetime

or, from a terminal: `matlab -batch "cd demo/Lifetime; runDemoLifetime"`. It runs
`runDemoLifetimeStacking.m` (part 2) and `runDemoLifetimeQuantify.m` (part 4) on the example data
and writes everything into `demo/Lifetime/demo_output/`. The expected console output, the expected
values and the run time are listed in the [Demo section of the repository
README](../../README.md#demo).

Neither the package nor the example data are modified by the demo: both wrappers assemble an
isolated working copy below `demo_output/`, which is not tracked by git.

## How the demo drives the module script

`lifeTcalculationbernsen.m` is shipped unchanged. It was written for one acquisition series: the
root folder, the number of images per field of view, the field-of-view folders it loops over and
the output file name (asked for with a save dialog) are hard-coded in it, and its file names carry
the `CFX` prefix of that series. The published example data set is the `WT50uM` series, whose
images are named `<FOV>-50uM-<cell>_ch<c>.ome.tif`.

Both wrappers therefore replay the section they need in an isolated working copy below
`demo_output/`:

* the published images are staged under the file names the script expects (`CFX-<cell>_ch<c>.ome.tif`
  for step 1; `CFXintensityStack.tif`, `CFXlifetimeStack.tif` and `CFXintensityStack-binary.tif`
  for step 2), and
* `lifetimeSectionDriver.m` copies the section out of the shipped script **verbatim** and replaces
  the acquisition settings of the demo. Every replacement is echoed to the console and the
  generated driver is kept in the working copy, so that what the demo runs can be inspected
  directly.

| Section | Replaced line | Replacement |
|---|---|---|
| 1 | `dirRoot = ['D:\ImageData\…'];` | the working copy |
| 1 | `fnum = 3;` | `fnum = 2;` — the number of images per field of view that the example data contains |
| 1 | `for jj = 1 : 3` | `for jj = [0 1 2]` — the field-of-view folders of the example data |
| 3 | `dirRoot = ['D:\ImageData\…'];` | the working copy |
| 3 | `filenameS = uiputfile('.mat','Save the results');` | `filenameS = 'LifetimeResults.mat';` — a fixed name, so that the demo does not wait for a click on *Save* |

The results are renamed back to the names of the published example data on collection, which also
corrects the typo `CFXXMlifetimeStack-filter.tif` (note the doubled `XM`) of the shipped script.

## Note on the Bernsen parameters

The two files that describe step 3 record different parameters:

* the section-2 comment of `lifeTcalculationbernsen.m`: window size 7 pixel, contrast threshold
  (`parameter_1`) 7;
* the shipped macro `bersenThtest.ijm`: a radius of 10–30 pixel (step 5) at a contrast threshold
  of 15.

The published mask itself is part of the example data (example data 3), so the demo does not depend
on which of the two records applies. The two records should be reconciled before the macro is used
on new data.

## Contents of this folder

| File | Role |
|---|---|
| `Readme.md` | this document |
| `lifeTcalculationbernsen.m` | module script — Section 1 image stacking, Section 2 Bernsen masking (in Fiji), Section 3 quantification |
| `bersenThtest.ijm` | Fiji macro of part 3 — Bernsen background mask |
| `runDemoLifetime.m` | demo entry point — runs the two unattended steps and prints a summary |
| `runDemoLifetimeStacking.m` | demo step 1 — stages the LAS X export and drives Section 1 on example data 1 |
| `runDemoLifetimeQuantify.m` | demo step 2 — stages the stacks and the mask and drives Section 3 |
| `lifetimeSectionDriver.m` | demo helper — copies a section of the module script verbatim into a driver and replaces the acquisition settings |
| `demo_data/` | the published example data set, items 1–3 (35 files, ≈ 52 MB) |
