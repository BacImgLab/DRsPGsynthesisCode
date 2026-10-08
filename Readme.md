# Multicolor Colocalization

`Colocalization/` contains the MATLAB scripts and the Fiji/ImageJ dependencies of the multicolor
fluorescence workflow for *Deinococcus radiodurans*. The workflow reconstructs population
demographs of the septum-associated signal along the cell cycle and quantifies the colocalization
between the fluorescence channels.

`demo/Colocalization/` is a **self-contained package** of that module: the 18 analysis scripts, the
two demo wrappers plus the entry point `runDemoColoc.m`, the published example data set in
`demo_data/` and this Readme all live in this single folder, so the folder can be downloaded and
run on its own.

The example data set used throughout this document is a colocalization experiment on the
FtsZ–FtsW pair, strain **JW425** (R1, *ftsW*::*ftsW*-linker0-ht str,
pJWK-P<sub>tet</sub>dr-*ftsZ*-linker0-mng). All folders referenced below are relative to the
`demo_data/` folder of this package, which ships the published example data set of the
Supplementary Information in six items. The folder names are the English translation of the six
original item names; the image and `.mat` files themselves are unchanged.

The module consists of four parts:

| Part | Content | Can be run unattended |
|---|---|---|
| 1 | [Data preprocessing](#1-data-preprocessing) | No — Fiji/ImageJ, Cellpose3 and DeCNN, all interactive |
| 2 | [Demograph construction](#2-demograph-construction) | No — the septal profiling needs septum lines drawn by hand |
| 3 | [S1 septum demograph analysis](#3-s1-septum-demograph-analysis) | **Yes** — covered by `runDemoColoc.m` |
| 4 | [Colocalization analysis (PCC)](#4-colocalization-analysis-pcc) | **Yes** — covered by `runDemoColoc.m` |

Parts 3 and 4 are the only two steps of the workflow that consume a plain, scriptable input, and
they are the two that the automated demo replays (see [Automated demo](#automated-demo)).

---

## 1. Data preprocessing

### Image acquisition

Four channels were acquired for every field of view:

| Channel | Excitation / label | Emission filter |
|---|---|---|
| BF | transmitted light, 1 % intensity | — |
| 488 | mNeonGreen (mNG) | ET525/50m |
| 561 | JF552-HT | ET605/70m |
| 647 | Potomac red, membrane marker | ET705/72m |

Each channel was recorded as a Z-series of five planes, centred on the midplane with a step of
0.2 µm (± 0.5 µm, i.e. a total Z depth of 1 µm per channel). In the processing scripts the
channels are numbered `channel1` = BF, `channel2` = 647, `channel3` = 561 and `channel4` = 488.

In addition, multicolor fluorescent beads were acquired in all three fluorescence channels. They
are the calibration reference for the chromatic aberration correction of step 1.2 below.

### 1.1 Drift correction and Z-projection

Lateral drift in the multichannel Z-stacks is corrected automatically and the corrected stacks are
Z-projected.

* MATLAB script: `MultiChannelDriftCorrection.m`
* Fiji plugin, called internally: `HyperStackReg`

The script initializes Fiji through the Miji interface, reads a configurable number of Z-planes per
channel (three by default), builds one Z-stack per channel and merges the channels into a
hyperstack. Drift correction is then run with `HyperStackReg` using channel 2 (647 nm) as the
reference channel. The drift-corrected grayscale planes are split and saved, and an
average-intensity Z-projection is computed and saved for channels 2 to 4
(`AVG_ChannelX_000-0_2.tif`).

The Java/Miji paths and the input folder are hard-coded at the top of the script and have to be
edited before use.

### 1.2 Chromatic aberration correction

Chromatic aberration displaces the fluorescence channels relative to each other. The displacement
is measured on the calibration beads and removed from the experimental data, aligning the 488 nm
and 647 nm channels onto the 561 nm reference channel.

* MATLAB script: `MultiChannelChromaticAberration.m`
* Dependency: `imreg2Dr.m`
* Requires user interaction: the square ROI on the bead images is selected manually.

The script loads the bead images of all positions for the three fluorescence channels, asks the
user to place a square ROI on the maximum projection of the 647 nm channel, crops and stacks that
ROI from every bead image, and computes the spatial transformation that aligns the 488 nm and
647 nm channels to the 561 nm channel. The transformation matrices are stored in `tforms.mat` and
applied to the experimental data sets, which are then written out as multi-page TIFFs.

### 1.3 Denoising

Every fluorescence channel is denoised separately to improve the signal-to-noise ratio before the
quantitative analysis.

* Fiji plugin: `PureDenoise` — run manually in Fiji.

### 1.4 Cell-cycle classification

The denoised, chromatic-aberration-corrected images are segmented and classified by cell-cycle
stage:

1. In **Cellpose3**, the membrane channel is segmented with the trained model, which yields a mask
   image for every single cell.
2. The mask of the membrane channel of every cell is classified with the neural-network model
   [BacImgLab/DeCNN](https://github.com/BacImgLab/DeCNN). The other three images of the same cell
   are stacked onto it, which gives one multichannel stack per cell and per cell-cycle stage
   (`stage1` … `stage5`).
3. The three independent experiments are merged.

→ produced: **example data 1** and **example data 2** (see [Example data](#example-data)).

---

## 2. Demograph construction

### 2.1 Channel merging

The Z-projected images of the different fluorescence channels are merged into one multichannel
image stack per cell, which is the input of the downstream analysis.

* MATLAB script: `combineSegfiles4C.m`

The script assumes four single-channel folders below a common root,
`basePath\C1\stageName`, `basePath\C2\stageName`, `basePath\C3\stageName` and
`basePath\C4\stageName`, and writes one four-page stack `basePath\stageName\Ch4_xxx.tif` in which
pages 1 to 4 are C1 to C4.

→ produced: **example data 3**. Channel merging and the demograph analysis that follows were
carried out on the stage-3, stage-4 and stage-5 cells; the shipped example data set also contains
the two stage-2 cells, which enter the workflow at the next step.

### 2.2 Septal signal profiling in individual cells

The fluorescence intensity distribution of the division septum (S1 and S0) is extracted from every
segmented single cell.

* MATLAB scripts: `DemoDR_stage35.m` (stages 3–5, septa S0, S1-left, S1-right) and
  `DemoDR_stage2.m` (stage 2, septum S0)
* Dependencies: `bf2rotateRec.m`, `DR4cShow.m`, `lineProfile.m`, `S0S1select.m`, `S0select.m`
* **Requires user interaction**: the S1 and S0 septum lines are drawn by hand in a MATLAB figure.

For every image the bright-field channel is used to compute the rotation angle and the bounding box
of the cell, and all fluorescence channels are rotated so that the cell is horizontally aligned.
The rotated multichannel stack is displayed for visual quality control and acceptance. If the cell
is accepted, the septum positions are selected interactively on the membrane channel — S0 (the
primary septum) and, in the stage-3–5 script, S1-left and S1-right (the secondary septa). The
intensity profiles along each septum are then extracted for all fluorescence channels with a
user-defined line width, together with the septal geometric parameters such as the cell width `D`
and the septum length `d`.

The output per cell is written into the `Processed` subfolder of the stage folder as
`Ch4_DR_*_data.mat` (septum geometry, S0/S1 features and intensity profiles) next to
`Ch4_DR_*_rot.tif` (the same cell with the septum rotated to vertical).

→ produced: **example data 4**.

This step is the one interactive step that sits in the middle of the workflow: everything upstream
of it runs inside Fiji, Cellpose3 or the DeCNN classifier, and everything downstream of it consumes
the per-cell profiles it produces. It is therefore not replayed by the automated demo.

### 2.3 Signal sorting and trend visualisation

All cells are sorted by the septum-related metric — the maturity of the septum — so that the
change of the fluorescence signal over the cell cycle becomes visible as a trend.

* MATLAB script: `BeforeDemoProcess.m` (note the capitalisation; the manuscript refers to it as
  `Beforedemoprocess`)

The script extracts and computes the intensity metrics, merges the left- and right-side data,
renames the fields for consistency, sorts the cells by their features, and finally extracts the
aligned and normalised channel intensity profiles. The results are written into the subfolders
`extracted_files`, `extracted_files/S1` and `extracted_files/S1/SortedData` of the working copy.

→ produced: **example data 5**. The demograph is a matrix in which every column is one cell,
sorted by septum maturity, and every row is one position along the septum profile; it is written
per side (S0, S1) and per channel (`dDratioS0Channel1..4.mat`, `dDratioS1Channel1..4.mat`).

### 2.4 Demograph smoothing

The demograph is smoothed to make the population-level trend clearer.

* MATLAB script: `demoSmooth.m`
* Dependency: `average_colN.m`

The script loads the demograph matrix of a channel, smooths it along the cell cycle by averaging
columns, normalises the result by its maximum, and saves the normalised matrix as a `.mat` file and
as a grayscale `.tif` image. It expects the variables `dsDChannel2new`, `dsDChannel3new` and
`dsDChannel4new` in the input `.mat` files and is run for channels 2, 3 and 4.

→ produced: **example data 6** (`demo_sm_norm2new.mat`, `demo_sm_norm3new.mat`,
`demo_sm_norm4new.mat`, 80 × 80 each; rows = position along the septum profile, columns = the 80
cell-cycle bins).

---

## 3. S1 septum demograph analysis

The distance between the signal peaks of two fluorescence channels in the S1 division septum is
quantified.

* MATLAB script: `DemoS1AnalysisW.m`
* Dependency: `demofitGauss2.m`
* **Input**: example data 6 — channel 2 (membrane) and channel 4 (protein) of the smoothed S1
  demograph
* **Output**: `WmemData.mat` and `Wmem.tif` — these two files are referred to as `AmemData` and
  `Amem.tif` in the manuscript text, see [Note on the file names](#note-on-the-file-names)
* This step is run unattended and is covered by the demo.

A two-component Gaussian model is fitted to every column of the two smoothed demograph matrices —
that is, to the S1 profile of every cell-cycle bin — and the position of the first Gaussian peak is
recorded for the protein channel (`mui_poi`) and for the membrane channel (`mu1_mem`). The
distance between the two ridges across the septum is the quantity reported for S1 in the
manuscript. `WmemData.mat` also contains `ParameterAll`, the fit parameters of every column.
`Wmem.tif` is a multi-page TIFF in which every page is the plot of one column with its raw data and
its two fitted curves.

Requires the Optimization Toolbox (`lsqcurvefit`) and the Image Processing Toolbox (`rgb2ind`).

---

## 4. Colocalization analysis (PCC)

The Pearson correlation coefficient (PCC) is computed between two fluorescence channels for every
single cell, as a quantitative measure of how strongly the two signals colocalize.

* MATLAB script: `PCCcaclu.m`
* **Input**: example data 4 — the per-cell septum profiles of stage2 … stage5
* **Output**: the PCC value between the different channels of every cell-cycle stage, plus the
  merged PCC value over all cells
* This step is run unattended and is covered by the demo.

The script has two parts. The first part walks the stage folders `stage2` … `stage5`, searches for
the rotated images `Ch4_DR_*_rot.tif` in each `Processed` subfolder, locates the matching unrotated
original `Ch4_DR_*.tif` in the stage folder and copies it into a `PCCcaclu` subfolder. The second
part reads the four-page stacks of those originals, masks the background (a pixel is kept only if
all of channels 2, 3 and 4 are non-zero) and computes the PCC of the pixel intensities between
channel 2 and 3 (`C2C3`), channel 2 and 4 (`C2C4`) and channel 3 and 4 (`C3C4`), together with the
mean intensity of each channel after removing the zero-value background pixels and subtracting 100.

Every stage is written to `stageN_PCC_Values.csv` with the columns
`FileName, C2C3, C2C4, C3C4, meanC2, meanC3, meanC4`, and all stages are concatenated into
`Merged_PCC_Values.csv`.

---

## Example data

`demo_data/` contains the six items of the published example data set, 90 files, ≈ 13 MB in total.

| Supplement item | Folder | Content | Read by the demo |
|---|---|---|---|
| 1. Preprocessed data | `01_preprocessed/` | Central 512 × 512 crop of the middle Z-plane of the four preprocessed channels (BF, 488, 561, 647), single precision. | — |
| 2. Cell-cycle classified data | `02_cell_cycle_stacks/` | The four single-channel images of the ten example cells, one folder per channel (`C1-BF`, `C2-647`, `C3-561`, `C4-488`) and one subfolder per cell-cycle stage (`stage1` … `stage5`), two cells per stage. 150 × 150 uint16. | — |
| 3. Channel-merged data | `03_channel_merged/` | The eight stage-2 … stage-5 cells merged into four-page stacks `Ch4_DR_*.tif` (pages 1..4 = C1..C4). | — |
| 4. Per-cell septum profiles | `04_septum_profiles/` | The output of the manual septal profiling: `Ch4_DR_*.tif` plus `Processed/Ch4_DR_*_rot.tif` (septum rotated to vertical) and `Processed/Ch4_DR_*_data.mat` (septum geometry, S0/S1 features, intensity profiles). Eight cells over stage2 … stage5. | step 2 (input) |
| 5. Sorted demograph | `05_sorted_demograph/` | The demographs: `S0/dDratioS0Channel1..4.mat` (80 × 1175) and `S1/dDratioS1Channel1..4.mat` (80 × 1100). One column per cell, sorted by septum maturity. | — |
| 6. Smoothed demograph | `06_smoothed_demograph/` | The smoothed demographs: `S0|S1/demo_sm_norm2..4new.mat` (80 × 80). | step 1 (input: `S1`, channels 2 and 4) |

`stage1` … `stage5` in the folder names above always denotes the **cell-cycle** stage. The six
example-data items are numbered instead, precisely to avoid the two meanings of the word *stage*.

Three notes on this data set:

* Only `04_septum_profiles/` and `06_smoothed_demograph/` are read by the demo; the other four
  folders are shipped because they are part of the published example data set, not because the
  demo needs them.
* The four files in `01_preprocessed/` are **crops**, not the complete preprocessing output. The
  full stacks are 1006 × 1006 × 16 slices in single precision, ≈ 62 MiB per channel — above the
  25 MiB single-file limit of the GitHub web uploader. The crop covers rows and columns 250 : 761
  of the middle Z-plane (slice 8 of 16) and keeps the original intensity values, dtype and scaling.
* `05_sorted_demograph/` and `06_smoothed_demograph/` describe the whole pooled population of the
  manuscript (1175 cells for S0, 1100 for S1), whereas only ten of those cells are shipped
  individually as images in `02`–`04`.

---

## Requirements

* MATLAB R2024a
* Fiji/ImageJ 1.54f with
  * `HyperStackReg` (version 5.7) for the drift correction,
  * `PureDenoise` for the denoising,
  * the Miji interface (ImageJ–MATLAB bridge) installed and configured for
    `MultiChannelDriftCorrection.m`
* Cellpose3 and [BacImgLab/DeCNN](https://github.com/BacImgLab/DeCNN) for the cell-cycle
  classification
* For `DemoS1AnalysisW.m` (part 3): the Optimization Toolbox (`lsqcurvefit`) and the Image
  Processing Toolbox (`rgb2ind`). Parts 1, 2 and 4 run on base MATLAB.
* Windows 10 / 11 (64-bit)

Input data are multichannel Z-stack fluorescence images (`.tif`) and the bright-field images used
for the segmentation, plus multicolor bead images for the chromatic aberration correction.

---

## How to use

1. Run `MultiChannelDriftCorrection.m` to apply the drift correction and compute the Z-projections
   (needs Fiji and Miji; edit the paths at the top of the script first).
2. Run `MultiChannelChromaticAberration.m` to correct the chromatic aberration; select the bead ROI
   when prompted.
3. Denoise every fluorescence channel with the Fiji `PureDenoise` plugin.
4. Segment the membrane channel in Cellpose3 and classify the cells into cell-cycle stages with the
   DeCNN model.
5. Merge the channels with `combineSegfiles4C.m`.
6. Extract the septal profiles from the segmented cells with `DemoDR_stage35.m` and
   `DemoDR_stage2.m`; draw the S0 and S1 septum lines when prompted, and accept or reject every
   cell in the quality-control window.
7. Sort the cells and build the demographs with `BeforeDemoProcess.m`, then smooth them with
   `demoSmooth.m`.
8. Quantify the ridge distance in the S1 septa with `DemoS1AnalysisW.m`.
9. Compute the PCC between the channels with `PCCcaclu.m`.

## Automated demo

Steps 8 and 9 — the only two steps with no user interaction — are covered by the demo that ships
with this folder:

    >> cd demo/Colocalization
    >> runDemoColoc

or, from a terminal: `matlab -batch "cd demo/Colocalization; runDemoColoc"`. It runs
`runDemoColocS1.m` (part 3) and `runDemoColocPcc.m` (part 4) on the example data and writes
everything into `demo/Colocalization/demo_output/`. The expected console output, the expected
values and the run time are listed in the [Demo section of the repository
README](../../README.md#demo).

Neither the package nor the example data are modified by the demo: both wrappers assemble an
isolated working copy below `demo_output/`, which is not tracked by git.

## Note on the file names

Two file names differ between the manuscript text and the scripts:

| In the manuscript text | Written by the script | Where |
|---|---|---|
| `AmemData`, `Amem.tif` | `WmemData.mat`, `Wmem.tif` | output of `DemoS1AnalysisW.m` |
| `dDratioS1Channel2`, `dDratioS1Channel4` | `demo_S1_smooth_mem.mat`, `demo_S1_smooth_W.mat` | input of `DemoS1AnalysisW.m` |

The demo follows the scripts: it feeds the published smoothed demographs `demo_sm_norm2new.mat`
(channel 2, membrane) and `demo_sm_norm4new.mat` (channel 4, protein) of `06_smoothed_demograph/S1/`
to part 3 under the two names the script expects, and it reports the output as `WmemData.mat` /
`Wmem.tif`.

## Contents of this folder

| File | Role |
|---|---|
| `Readme.md` | this document |
| `MultiChannelDriftCorrection.m` | 1.1 — drift correction and Z-projection (Fiji `HyperStackReg` via Miji) |
| `MultiChannelChromaticAberration.m` | 1.2 — chromatic aberration correction from bead images |
| `imreg2Dr.m` | intensity-based 2D registration, used by 1.2 |
| `splitchannelsafterCACorrection.m` | splits the corrected multi-channel TIFFs back into single-channel images |
| `combineSegfiles4C.m` | 2.1 — merges four single-channel images into one four-page stack |
| `DemoDR_stage35.m` | 2.2 — septal profiling for stage-3 to stage-5 cells (S0, S1-left, S1-right) |
| `DemoDR_stage2.m` | 2.2 — septal profiling for stage-2 cells (S0) |
| `bf2rotateRec.m` | orientation and rotated bounding box of a cell from the bright-field image |
| `DR4cShow.m` | displays a 4-channel cell for visual acceptance |
| `S0S1select.m` | interactive selection of the S0 and S1 septum lines |
| `S0select.m` | interactive selection of the S0 septum line |
| `lineProfile.m` | intensity profile along a line, with averaging over a line width |
| `BeforeDemoProcess.m` | 2.3 — metrics, sorting and demograph construction |
| `demoSmooth.m` | 2.4 — column-wise smoothing and normalisation of a demograph |
| `average_colN.m` | averages the columns of a matrix into `N` bins, used by 2.4 |
| `DemoS1AnalysisW.m` | 3 — two-Gaussian fit of every S1 column, ridge distance |
| `demofitGauss2.m` | two-component Gaussian fit, used by 3 |
| `PCCcaclu.m` | 4 — pairwise Pearson correlation between the channels of every cell |
| `runDemoColoc.m` | demo entry point — runs the two unattended steps (3 and 4) in sequence and prints a summary |
| `runDemoColocS1.m` | demo step 1 — assembles an isolated working copy and drives `DemoS1AnalysisW.m` on example data 6 |
| `runDemoColocPcc.m` | demo step 2 — assembles an isolated working copy and drives `PCCcaclu.m` on example data 4 |
| `demo_data/` | the published example data set, items 1–6 (90 files, ≈ 13 MB) |
