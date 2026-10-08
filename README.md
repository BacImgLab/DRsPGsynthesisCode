# DRsPGsynthesisCode

`DRsPGsynthesisCode` is a collection of MATLAB scripts and Fiji/ImageJ-based workflows for quantitative analysis of *Deinococcus radiodurans* imaging data.

The repository provides analysis workflows for:

* Multicolor fluorescence colocalization and demography analysis
* Fluorescence lifetime imaging microscopy (FLIM)
* Single-molecule tracking (SMT)
* 3D-dSTORM septum width analysis
* Electron tomography (ET) septum width analysis

The repository contains analysis scripts and workflows that use Fiji/ImageJ and other external software for image preprocessing, segmentation, registration, localization, trajectory analysis, and quantitative measurements.

Four self-contained, ready-to-run demo packages are provided in `demo/`, one per analysis workflow
that has unattended steps: `demo/Colocalization/` (multicolor colocalization), `demo/Lifetime/`
(fluorescence lifetime imaging), `demo/SMT/` (single-molecule tracking) and `demo/SMLM/`
(3D-dSTORM septum width). Each runs the analysis steps that need no user interaction on the
published example data set of that workflow, without Fiji/ImageJ, without acquisition hardware and
without the raw imaging data (see [Demo](#demo)).
The electron tomography workflow is interactive in all of its steps and ships as a documentation
package instead: `demo/Tomo/` contains the module script, the published example data set and the
step-by-step walkthrough, but no runnable entry point (see [Demo](#demo)).

---

## Repository structure

The main analysis modules are organized according to the type of imaging data:

| Module | Directory |
|---|---|
| 1. Multicolor Colocalization Imaging | `Colocalization/` |
| 2. Fluorescence Lifetime Imaging (FLIM) | `Lifetime/` |
| 3. Single-Molecule Tracking (SMT) | `SMTanalysis/` |
| 4. 3D-dSTORM Septum Width Analysis | `3D-dSTORM/` |
| 5. Electron Tomography Septum Width Analysis | `Tomo/` |

Each module directory contains its own `Readme.md` with module-specific instructions.

Four ready-to-run demo packages are provided in `demo/`, one for the multicolor colocalization
workflow (`demo/Colocalization/`), one for the fluorescence-lifetime-imaging workflow
(`demo/Lifetime/`), one for the single-molecule-tracking workflow (`demo/SMT/`) and one for the
3D-dSTORM septum-width workflow (`demo/SMLM/`); each ships the module scripts it drives, the
published example data set of that workflow and its own `Readme.md` (see [Demo](#demo)). A fifth,
documentation-only package, `demo/Tomo/`, ships the module script of the electron tomography
workflow together with its published example data set and walkthrough; it has no runnable entry
point because both steps of that workflow are interactive.

---

## 1. Multicolor Colocalization Imaging

Directory: `Colocalization/`

This module is used to preprocess and quantitatively analyze multicolor fluorescence images of *D. radiodurans*.

### Workflow

The general workflow includes:

1. Multichannel image preprocessing
2. Image drift correction
3. Chromatic aberration correction
4. Image denoising
5. Cell/septum segmentation
6. Extraction of septum-associated fluorescence signals and morphological features
7. Demography or cell-cycle-stage analysis
8. Quantification of fluorescence correlations between channels

### Main scripts

* `MultiChannelDriftCorrection.m`
  Performs multichannel drift correction using the Fiji/ImageJ HyperStackReg plugin.

* `MultiChannelChromaticAberration.m`
  Performs chromatic aberration correction using `imreg2Dr.m`.

* `combineSegfiles4C.m`
  Combines segmentation and quantitative data for multichannel analysis.

* `BeforeDemoProcess.m`
  Performs preprocessing before demography analysis.

* `DemoDR_stage35.m`
  Performs demography analysis according to *D. radiodurans* cell-cycle stages.

* `demoSmooth.m`
  Performs smoothing of quantitative measurements.

* `PCCcaclu.m`
  Calculates Pearson correlation coefficients between fluorescence channels.

### Fiji/ImageJ dependencies

Depending on the analysis workflow, the following Fiji/ImageJ tools may be required:

* HyperStackReg for image registration/drift correction
* PureDenoise for image denoising
* `imreg2Dr.m` for chromatic aberration correction

---

## 2. Fluorescence Lifetime Imaging (FLIM)

Directory: `Lifetime/`

This module is used to process fluorescence lifetime and intensity images exported from Leica LAS X FLIM/FCS software.

### Workflow

1. Acquire FLIM data using a Leica STELLARIS 8 confocal microscope.
2. Export intensity and lifetime images from Leica LAS X FLIM/FCS software.
3. Perform image stacking and preprocessing using MATLAB.
4. Perform background subtraction using the Fiji/ImageJ `bersenThtest` macro.
5. Calculate pixel-wise fluorescence lifetime and intensity values.
6. Export quantitative measurements for downstream statistical analysis.

### Main script

* `lifeTcalculationbernsen.m`
  Stacks exported FLIM images, performs image processing and background correction, and calculates pixel-wise fluorescence lifetime and intensity values.

### Fiji/ImageJ dependency

* `bersenThtest.ijm` for background subtraction

### Input

Intensity and lifetime image files exported from Leica LAS X FLIM/FCS software.

### Output

Processed fluorescence lifetime and intensity measurements for downstream quantitative analysis.

---

## 3. Single-Molecule Tracking (SMT)

Directory: `SMTanalysis/`

This module is used to process and analyze single-molecule imaging data, including localization, trajectory generation, trajectory segmentation, and molecular-state classification.

### Workflow

1. Prepare and crop regions of interest from raw images.
2. Perform single-molecule localization using ThunderSTORM in Fiji/ImageJ.
3. Correct chromatic aberration where required.
4. Segment cells using Cellpose3.
5. Link localized spots into molecular trajectories.
6. Segment trajectories into different movement states.
7. Classify molecular states.
8. Calculate velocity and trajectory distributions.
9. Remove trajectories or regions affected by drift where required.
10. Generate trajectory visualizations and quantitative measurements.

### Main scripts

* `SMTdataPrepare.m`
  Prepares and crops SMT image data and defines regions of interest.

* `CACorrectionofSMTsplitter.m`
  Performs chromatic aberration correction for split-channel SMT data.

* `spotsLinking.m`
  Links localized single-molecule positions into trajectories.

* `RefineTraceSegDr.mlapp`
  App Designer app for interactive trajectory segmentation and refinement.

* `statesClassifyDr.m`
  Classifies molecular trajectories into movement states.

* `dataprocessSMTWCF.m`
  Performs velocity distribution fitting and related quantitative analysis.

* `tracePlotafterseg.m`
  Generates trajectory visualizations after trajectory segmentation.

* `RemoveDrift.m`
  Removes regions or trajectories affected by sample drift.

* `MSDsingle2D.m`
  Calculates mean squared displacement (MSD) curves from trajectories.

* `CDF_logCalc.m`
  Calculates step-length cumulative distribution functions.

* `linfitR.m`
  Performs linear fitting of MSD/speed data.

### External software

The SMT workflow uses:

* Fiji/ImageJ
* ThunderSTORM for single-molecule localization
* Cellpose3 for cell segmentation

---

## 4. 3D-dSTORM Septum Width Analysis

Directory: `3D-dSTORM/`

This module is used to quantify septum width from 3D-dSTORM imaging data.

### Workflow

1. Generate 2D projections of 3D-dSTORM localization data.
2. Use Fiji/ImageJ to manually define septum regions of interest (ROIs).
3. Import the ROI and image information into MATLAB.
4. Extract fluorescence intensity profiles across the septum.
5. Calculate septum width from the intensity profile.
6. Determine the full width at half maximum (FWHM).

### Main script

* `dSTORMwidthcaclu.m`
  Calculates septum width from fluorescence intensity profiles and determines the FWHM.

### Input

* 2D projections of 3D-dSTORM images
* Manually defined septum ROIs

### Output

Quantitative septum-width measurements based on FWHM analysis.

Of the two steps of this workflow only the second one runs without a human — the line ROIs are
drawn by hand in Fiji — so the demo starts from the published example data set (the exported
intensity profiles) and replays the FWHM calculation: [`demo/SMLM/`](demo/SMLM/Readme.md).

---

## 5. Electron Tomography Septum Width Analysis

Directory: `Tomo/`

This module is used to quantify septum width from electron tomography images.

### Workflow

1. Generate or obtain 2D images from electron tomography data.
2. Open the image in Fiji/ImageJ.
3. Manually define the septum ROI.
4. Use `TomowidthCaclu.m` to extract the intensity profile along the selected ROI.
5. Calculate the septum width from the resulting intensity profile.

### Main script

* `TomowidthCaclu.m`
  Calculates septum width from the intensity profile along the manually selected septum ROI.

### Input

* 2D electron tomography images
* Manually defined septum ROIs

### Output

Quantitative measurements of septum width.

Both steps of this workflow are interactive — the septum ROI is drawn by hand in Fiji, and
`TomowidthCaclu.m` asks for two clicks in a MATLAB figure window to define the main axis of the
septum — so the module ships no automated demo. The published example data set (the tomographic
slices, the drawn ROIs and the resulting width measurements) and the step-by-step walkthrough are
available in [`demo/Tomo/`](demo/Tomo/Readme.md).

---

## Requirements

### Software

The analysis workflows use the following software and plugins, depending on the analysis module:

* MATLAB R2024a
* Fiji/ImageJ (ImageJ 1.54f)
* Cellpose3
* ThunderSTORM (dev-2016-09-04-b1)
* HyperStackReg (Version 5.7)
* PureDenoise
* For the S1 ridge-distance step of the colocalization demo (`Colocalization/DemoS1AnalysisW.m`): the Optimization Toolbox (`lsqcurvefit`) and the Image Processing Toolbox (`rgb2ind`). Every other script and demo step runs on base MATLAB.

### Operating system

* Windows 10 / 11 (64-bit)

### Tested environment

The analysis workflows were tested using **MATLAB R2024a on Windows 11 (64-bit)**, together with ImageJ 1.54f, Cellpose3 and ThunderSTORM (dev-2016-09-04-b1).

### Hardware

No non-standard hardware is required to run the analysis scripts.

The imaging instruments used to acquire the experimental datasets are described in the corresponding Methods section of the manuscript and are not required for running the analysis code on existing image datasets.

---

## Installation

### 1. Obtain the repository

Download or clone this repository to a local computer.

### 2. Install MATLAB

Install MATLAB R2024a.

### 3. Install Fiji/ImageJ

Install Fiji/ImageJ and the required plugins for the selected analysis workflow.

### 4. Install Cellpose3

Install Cellpose3 according to its standard installation procedure.

Cellpose3 is required for cell segmentation in the SMT workflow.

### 5. Install ThunderSTORM

Install the ThunderSTORM plugin in Fiji/ImageJ.

ThunderSTORM is required for single-molecule localization in the SMT workflow.

### 6. Install additional Fiji/ImageJ tools

Install the tools required by the selected workflow, including:

* HyperStackReg
* PureDenoise

Not all tools are required for every analysis module.

### 7. MATLAB scripts

No compilation is required. The scripts in `Colocalization/`, `Lifetime/`, `SMTanalysis/`, `3D-dSTORM/` and `Tomo/` are plain MATLAB `.m` files (with one App Designer app, `SMTanalysis/RefineTraceSegDr.mlapp`).

Open the MATLAB script corresponding to the desired analysis module and specify the input data path and analysis parameters before running the script.

---

## Typical installation time

Obtaining the repository and preparing the MATLAB scripts takes less than 5 minutes on a normal desktop computer, provided that MATLAB and the required Fiji/ImageJ plugins for the selected module are already installed.

Installing MATLAB, Fiji/ImageJ and the plugins from scratch may take longer; the required time depends on the user's computer, network connection and existing software environment.

Running the demos described in the [Demo](#demo) section requires only MATLAB and takes about three to four minutes in total (≈ 60 s for the colocalization demo, ≈ 10 s for the FLIM demo, ≈ 20 s for the SMT demo and ≈ 105 s for the 3D-dSTORM demo).

---

## How to use

Select the analysis module corresponding to the imaging dataset.

### General workflow

1. Prepare the input data according to the requirements of the selected module.
2. Perform the required Fiji/ImageJ preprocessing.
3. Open the corresponding MATLAB script.
4. Specify the input data path.
5. Set the required analysis parameters.
6. Run the script.
7. Inspect and export the resulting quantitative measurements.
8. Perform downstream statistical analysis as described in the manuscript.

Module-specific input requirements and processing steps are described in the corresponding `Readme.md` file of each module directory and in the MATLAB scripts themselves.

---

## Data preparation

The required input data depend on the analysis module.

### Multicolor colocalization

**Input:**

* Multichannel fluorescence images
* Corresponding segmentation information where required

**Preprocessing may include:**

* Drift correction
* Chromatic aberration correction
* Denoising
* Cell/septum segmentation

### FLIM

**Input:**

* Lifetime images
* Intensity images exported from Leica LAS X FLIM/FCS software

### SMT

**Input:**

* Single-molecule image sequences
* Localization results generated using ThunderSTORM
* Cell segmentation results generated using Cellpose3 where required

### 3D-dSTORM

**Input:**

* 2D projections of 3D-dSTORM localization data
* Manually defined septum ROIs

The published example data set of this workflow — the `X,Y` intensity profiles of 116 WT septum
ROIs, one folder per cell-cycle stage (S0: 51, S1: 65) — ships in
`demo/SMLM/demo_data/01_wt_s0_s1_profiles/`.

### Electron tomography

**Input:**

* 2D electron tomography images
* Manually defined septum ROIs

The published example data set of this workflow — two tomographic slices, the two Fiji polygon ROIs
and the resulting width measurements — ships in `demo/Tomo/demo_data/`.

---

## Output

The scripts generate quantitative measurements for downstream statistical analysis and figure preparation.

Depending on the analysis module, outputs may include:

* Fluorescence intensity measurements
* Fluorescence lifetime measurements
* Pearson correlation coefficients
* Cell-cycle-stage-associated fluorescence measurements
* Single-molecule trajectories
* Molecular movement-state classifications
* Velocity distributions
* Septum width measurements
* FWHM measurements
* Tomography-derived septum width measurements

The exact output format depends on the individual MATLAB script.

---

## Reproducibility

The analysis code in this repository was used to process imaging data described in the associated manuscript.

### Testing the code without experimental data

Four self-contained demo packages are included in `demo/`, and none of them needs Fiji/ImageJ or the acquisition hardware:

* `demo/Colocalization/runDemoColoc.m` runs the two steps of the multicolor colocalization pipeline that need no user interaction, on the **published example data set** shipped in `demo/Colocalization/demo_data/`: the S1 septum demograph analysis (`Colocalization/DemoS1AnalysisW.m`) and the pairwise Pearson correlation between the channels (`Colocalization/PCCcaclu.m`).
* `demo/Lifetime/runDemoLifetime.m` runs the two steps of the FLIM pipeline that need no user interaction, on the **published example data set** shipped in `demo/Lifetime/demo_data/`: the image stacking (`Lifetime/lifeTcalculationbernsen.m`, Section 1) and the per-pixel lifetime and intensity quantification (Section 3).
* `demo/SMT/runDemoSMT.m` runs the one step of the single-molecule-tracking pipeline that needs no user interaction, on the **published example data set** shipped in `demo/SMT/demo_data/`: the log-normal fitting of the speed distribution of the directed FtsW segments (`SMTanalysis/dataprocessSMTWCF.m`).
* `demo/SMLM/runDemoSMLM.m` runs the one step of the 3D-dSTORM septum-width pipeline that needs no user interaction, on the **published example data set** shipped in `demo/SMLM/demo_data/`: the FWHM calculation of every septum profile (`3D-dSTORM/dSTORMwidthcaclu.m`).

The upstream parts of all four workflows are interactive by design and are not replayed; see the [Demo](#demo) section below.

The electron tomography workflow has no unattended step at all — the ROI is drawn by hand in Fiji and `Tomo/TomowidthCaclu.m` asks for two clicks in a figure window — so it ships no runnable demo; its published example data set and walkthrough are in `demo/Tomo/` instead.

See the [Demo](#demo) section below for the exact commands and the expected output.

### Reproducing the quantitative results

To reproduce an analysis:

1. Obtain the corresponding imaging dataset (see the data-availability statement of the manuscript).
2. Follow the preprocessing workflow described for the relevant analysis module.
3. Run the corresponding MATLAB script.
4. Use the analysis parameters specified in the script and/or manuscript Methods.
5. Apply the same data-selection and segmentation criteria described in the manuscript.

| Quantitative result (figure / table) | Analysis module | Script(s) |
|---|---|---|
| Fig. 1b, c | Demography / cell-cycle staging | `Colocalization/BeforeDemoProcess.m` → `Colocalization/DemoDR_stage35.m` (cell-cycle classification: [BacImgLab/DeCNN](https://github.com/BacImgLab/DeCNN)) |
| Fig. 1c–e, Fig. 4a–e; Supplementary Table 5 | ET septum width | `Tomo/TomowidthCaclu.m` |
| Fig. 2a; Supplementary Fig. 7a, c | Colocalization (Pearson correlation) | `Colocalization/PCCcaclu.m` |
| Fig. 2b–d | Septal enrichment, demographs, displacement | `Colocalization/S0S1select.m`, `Colocalization/lineProfile.m`, `Colocalization/demoSmooth.m`, `Colocalization/DemoS1AnalysisW.m` |
| Fig. 3, Fig. 6a–d; Supplementary Fig. 9, 10f–g | FLIM lifetime and intensity | `Lifetime/lifeTcalculationbernsen.m` |
| Fig. 4f–g; Supplementary Table 5 | 3D-dSTORM / 3D-SMLM septum width | `3D-dSTORM/dSTORMwidthcaclu.m` |
| Fig. 5a–g; Supplementary Fig. 15–17 | SMT: localization, linking, state classification, MSD, speed distributions | `SMTanalysis/SMTdataPrepare.m` → ThunderSTORM → `SMTanalysis/spotsLinking.m` → `SMTanalysis/RefineTraceSegDr.mlapp` → `SMTanalysis/statesClassifyDr.m` → `SMTanalysis/dataprocessSMTWCF.m` and `SMTanalysis/MSDsingle2D.m` / `SMTanalysis/CDF_logCalc.m` (MSD, step-length CDF and velocity fits: `SMTanalysis/linfitR.m`) |
| Fig. 5d; Supplementary Fig. 16d–g | 2D-projection correction of FtsW speed | `SMTanalysis/unwrapTraj.m`, `SMTanalysis/DrunwrapX.m`, `SMTanalysis/DrunwrapY.m` |

Raw imaging datasets (single-molecule localization movies, FLIM photon data and electron-tomography tilt series) are not included in this repository because of their very large size; they are available from the corresponding author upon reasonable request. The multicolor colocalization, the FLIM, the single-molecule-tracking, the 3D-dSTORM and the electron-tomography modules are instead shipped with the published example data set of the Supplementary Information (see [Demo](#demo); for the tomography the shipped slices are the 2D images the ROIs were drawn on, not the tilt series).

---

## Demo

`demo/` contains four self-contained, ready-to-run demo packages, one per analysis workflow that has unattended steps, plus one documentation package for the electron tomography workflow. None of them needs acquisition hardware or raw imaging data, and none of them writes into the repository tree: all output goes to a `demo_output/` folder, which is not tracked by git.

Each package lives in its own subfolder together with the module scripts it drives, the published example data set of that workflow and the module `Readme.md`, so the folder can be downloaded and used on its own.

| Demo | Entry point | Data | MATLAB requirements | Run time |
|---|---|---|---|---|
| Multicolor colocalization | `demo/Colocalization/runDemoColoc.m` | published example data set (`demo/Colocalization/demo_data/`) | step 2 runs on base MATLAB; step 1 also needs the Optimization and Image Processing Toolboxes | ≈ 1 min |
| Fluorescence lifetime imaging | `demo/Lifetime/runDemoLifetime.m` | published example data set (`demo/Lifetime/demo_data/`) | base MATLAB only | ≈ 10 s |
| Single-molecule tracking | `demo/SMT/runDemoSMT.m` | published example data set (`demo/SMT/demo_data/`) | Optimization and Statistics and Machine Learning Toolboxes | ≈ 20 s |
| 3D-dSTORM septum width | `demo/SMLM/runDemoSMLM.m` | published example data set (`demo/SMLM/demo_data/`) | base MATLAB only | ≈ 2 min |
| Electron tomography | — (documentation only) | published example data set (`demo/Tomo/demo_data/`) | Mapping Toolbox (`polyxpoly`) if the script is run by hand | — |

### Single-molecule-tracking example data

`demo/SMT/demo_data/` holds the item of the published example data set that the SMT demo reads, in
the numbering of the Supplementary example data:

| Folder | Content | Read by the demo |
|---|---|---|
| `04_classified_traces/` | `FtsW-all.mat` — the output of the state-classification step: `SegVSRPTC` (1405 × 7 = velocity, StD, R, P, dwell time, location flag, state), the classification thresholds (`Pmin = 0.8`, `Rmax1 = 0.5`, `Rmax2 = 0.3`) and the speed and dwell-time vectors `Vd`/`Vs`/`Td`/`Ts` (487 directional and 917 stationary segments). | step 2 (input) |

The other three items of the published SMT example data set — the ThunderSTORM localisation table,
the linked trajectories and the segmented trajectories — are not read by the demo and are not
shipped; they are available from the corresponding author. See `demo/SMT/Readme.md` for the full
workflow, the description of all four items and the note on the variable names (`DataSMT` vs
`SegVSRPTC`).

### Multicolor colocalization example data

`demo/Colocalization/demo_data/` holds the example data set that accompanies the multicolor colocalization workflow, in the six items of the Supplementary example data. The folder names are the English translation of the original item names; the image and `.mat` files themselves are unchanged. 90 files, ≈ 13 MB in total.

| Folder | Content | Read by the demo |
|---|---|---|
| `01_preprocessed/` | Central 512 × 512 crop of the middle Z-plane of the four preprocessed channels (BF, 488, 561, 647), single precision. | — |
| `02_cell_cycle_stacks/` | The four single-channel images of the eight example cells, one folder per channel (`C1-BF`, `C2-647`, `C3-561`, `C4-488`) and one subfolder per cell-cycle stage (`stage1` … `stage5`). 150 × 150 uint16. | — |
| `03_channel_merged/` | The same eight cells merged into four-page stacks (`Ch4_DR_*.tif`, pages 1..4 = C1..C4). | — |
| `04_septum_profiles/` | Per-cell output of the manual septum profiling: `Ch4_DR_*.tif` plus `Processed/Ch4_DR_*_rot.tif` (septum rotated to vertical) and `Processed/Ch4_DR_*_data.mat`. | step 2 (input) |
| `05_sorted_demograph/` | The demographs: `S0/dDratioS0Channel1..4.mat` (80 × 1175) and `S1/dDratioS1Channel1..4.mat` (80 × 1100). One column per cell, sorted by septum maturity. | — |
| `06_smoothed_demograph/` | The smoothed demographs: `S0|S1/demo_sm_norm2..4new.mat` (80 × 80). | step 1 (input: `S1` channels 2 and 4) |

The folder names follow the numbering of the Supplementary example data, because the workflow itself uses "stage" for two different things: the six example-data items above, and the **cell-cycle** stages `stage1` … `stage5`, which appear as subfolders of `02_cell_cycle_stacks` and of `03`/`04`. Only the cell-cycle folders keep the word *stage*.

Three notes on this data set:

* Only `04_septum_profiles/` and `06_smoothed_demograph/` are read by `runDemoColoc`; the other four folders are shipped because they are part of the published example data set, not because the demo needs them.
* The four files in `01_preprocessed/` are **crops**, not the complete preprocessing output. The full stacks are 1006 × 1006 × 16 slices in single precision, ≈ 62 MiB per channel — above the 25 MiB single-file limit of the GitHub web uploader and far above a *small* dataset. The crop covers rows and columns 250 : 761 of the middle Z-plane (slice 8 of 16) and keeps the original intensity values, dtype and scaling.
* `05_sorted_demograph/` and `06_smoothed_demograph/` describe the whole pooled population of the manuscript (1175 cells for S0, 1100 for S1), whereas only eight of those cells are shipped individually as images in `02`–`04`.

### A note on the file names

Two file names differ between the manuscript text and the scripts:

| In the manuscript text | Written by the script | Where |
|---|---|---|
| `AmemData`, `Amem.tif` | `WmemData.mat`, `Wmem.tif` | output of `Colocalization/DemoS1AnalysisW.m` |
| `dDratioS1Channel2`, `dDratioS1Channel4` | `demo_S1_smooth_mem.mat`, `demo_S1_smooth_W.mat` | input of `Colocalization/DemoS1AnalysisW.m` |

The demo follows the scripts: it feeds the published smoothed demographs `demo_sm_norm2new.mat` (channel 2, membrane) and `demo_sm_norm4new.mat` (channel 4, protein) of `06_smoothed_demograph/S1/` to step 1 under the two names the script expects, and it reports the output as `WmemData.mat` / `Wmem.tif`.

### FLIM example data

`demo/Lifetime/demo_data/` holds the published example data set of the FLIM workflow in the three
items of the Supplementary example data. The folder names are the English translation of the
original item names; the image and metadata files themselves are unchanged. 35 files, ≈ 52 MB in
total.

| Folder | Content | Read by the demo |
|---|---|---|
| `01_lasx_export/` | The LAS X FLIM/FCS export: for every field of view (`0`, `1`, `2`) and every cell (1, 2) one intensity image (`_ch0.ome.tif`) and one lifetime image (`_ch1.ome.tif`), 1024 × 1024 uint16, plus the LAS X `MetaData/` folder. | step 1 (input) |
| `02_stacked/` | `WT50uMintensityStack.tif` and `WT50uMlifetimeStack.tif`, 6 pages of 1024 × 1024 uint16 — the output of the image-stacking step. | step 1 (reference), step 2 (input) |
| `03_background_masked/` | `WT50uMintensityStack-binary.tif` (the Bernsen mask, 8-bit, 0/255), `WT50uMintensityStack-filter.tif` and `WT50uMlifetimeStack-filter.tif` (the background-free stacks). | step 2 (input and reference) |

`02_stacked/` and `03_background_masked/` serve a second purpose in the demo: they are the reference the wrappers compare their own output with, so the demo shows that it reproduces the published stacks and the published background-free images page by page.

### Electron-tomography example data

`demo/Tomo/demo_data/` holds the two items of the published example data set of the electron
tomography workflow, under the English name of the original item. 8 files, ≈ 1.5 MB in total.

| Folder | Content | Role in the workflow |
|---|---|---|
| `01_manual_septum_roi/` | Two 2D tomographic slices of WT cells (963 × 1126 and 933 × 1143 pixel, 8-bit RGB JPEG with a 200 nm scale bar) and two Fiji polygon ROIs drawn around their septa (`...S0Roi.roi`, 54 vertices; `...S1Roi_1.roi`, 30 vertices). | step 1 (output) / step 2 (input) |
| `02_septum_width_calculation/` | The published script output for the two cells: the annotated figures (1609 × 4300 JPEG) and the result files (`Res` = distance along the middle line, local width, midpoint x and y, in nm; plus `CoordXY`). | step 2 (output) |

The two published result files describe a septum width of 27.9 ± 3.7 nm (range 22.7–51.7 nm, over
714 nm of axis, 715 samples, `S0Roi`) and 27.7 ± 8.6 nm (range 10.3–76.4 nm, over 126 nm of axis,
125 samples, `S1Roi_2`).

There is **no runnable entry point** in `demo/Tomo/`: the septum ROI is drawn by hand in Fiji, and
`TomowidthCaclu.m` itself asks for two clicks in a MATLAB figure window (`ginput(2)`) to define the
main axis of the septum, so an unattended replay would have to substitute the very step the module
exists for. `demo/Tomo/Readme.md` documents both steps in full, and carries two notes a reader
should be aware of: the shipped `...S1Roi_1.roi` and the published `...S1Roi_2.mat`/`.jpg` are not a
matched pair (the result files were produced from a second ROI that is not part of the example
data), and the published `...S1Roi_2.mat` carries five columns per row of `Res` where the current
script writes four.

### 3D-dSTORM septum-width example data

`demo/SMLM/demo_data/` holds the published example data set of the 3D-dSTORM septum-width workflow
under the English name of the original item. 116 files, ≈ 330 KB in total.

| Folder | Content | Read by the demo |
|---|---|---|
| `01_wt_s0_s1_profiles/WTS0/` | 51 septum profiles of the S0 cells: `WT_S0_width_roiN_5nmGuass_500mW_roi1_z450_650_Values.csv`, an `X,Y` header plus one row per sample (`X` = distance along the line ROI, `Y` = grey value). | step (input) |
| `01_wt_s0_s1_profiles/WTS1/` | 65 septum profiles of the S1 cells, same format. | step (input) |

These files are the Fiji export of step 1 of the workflow — the intensity profile along every
manually drawn line ROI — so they are exactly what the FWHM step consumes. All profiles are
sampled on the same 0.005 µm grid (5 nm); the line ROIs span 0.59–1.21 µm.

### Instructions to run on the demo data

All four demos are packaged in their own folder, so that each folder is self-contained and can also be
run after downloading that folder alone:

    >> cd demo/Colocalization
    >> runDemoColoc

    >> cd demo/Lifetime
    >> runDemoLifetime

    >> cd demo/SMT
    >> runDemoSMT

    >> cd demo/SMLM
    >> runDemoSMLM

All four can also be started from a terminal, e.g. `matlab -batch "cd demo/SMT; runDemoSMT"`.

`demo/Tomo/` has no entry point to run: follow the manual walkthrough in `demo/Tomo/Readme.md` to
reproduce the published tomography results by hand.

`demo/Colocalization/runDemoColoc.m` executes the two steps of the multicolor workflow that run unattended, driving the module scripts that the folder carries alongside the entry point:

1. `runDemoColocS1.m` — **S1 septum demograph analysis**. It feeds channels 2 (membrane) and 4 (protein) of the smoothed S1 demograph to `DemoS1AnalysisW.m`, which fits a two-component Gaussian model to each of the 80 columns and records the first peak position of both profiles, so that the distance between the protein ridge and the membrane ridge across the septum can be measured.
2. `runDemoColocPcc.m` — **PCC analysis**. It assembles an isolated working copy of the per-cell septum profiles (example data 4) and runs `PCCcaclu.m`, which computes the pairwise Pearson correlation coefficients of channels 2, 3 and 4 for every single cell and writes one table per cell-cycle stage plus the merged one.

`demo/Lifetime/runDemoLifetime.m` executes the two steps of the FLIM workflow that run unattended, driving the module script that the folder carries alongside the entry point:

1. `runDemoLifetimeStacking.m` — **image stacking**. It stages the published LAS X export (example data 1) under the file names `lifeTcalculationbernsen.m` expects and replays **Section 1** of that script, which reads the intensity image (`_ch0`) and the lifetime image (`_ch1`) of every cell of every field of view and appends them to an intensity and a lifetime stack.
2. `runDemoLifetimeQuantify.m` — **lifetime and intensity quantification**. It stages the stacks (example data 2) together with the published Bernsen mask (example data 3) and replays **Section 3** of the script, which sets every background pixel of both stacks to zero, rescales the lifetime grey values to nanoseconds, collects the remaining pixels per page, builds the lifetime histogram on the grid 0.30 : 0.03 : 3.60 ns and computes the mean and the standard deviation of both signals.

`demo/SMT/runDemoSMT.m` executes the one step of the SMT workflow that runs unattended:

1. `runDemoSMTVelocity.m` — **speed-distribution fitting**. It assembles an isolated working copy of the classified trajectory table (example data 4), adds the variable alias the module script expects (see the note on the file names in `demo/SMT/Readme.md`), and runs `dataprocessSMTWCF.m`, which computes the empirical CDF of the directed-segment speeds above 1 nm/s on 21 log-spaced bins between 1 and 100 nm/s, fits a single- and a double-log-normal population by least squares, bootstraps both fits over 200 resamples and reconstructs the two probability density functions for the comparison with the speed histogram.

The other steps of all three workflows are **not replayed**, because they need a human or an external GUI. For the multicolor colocalization workflow these are the drift correction and the chromatic-aberration correction in Fiji (`MultiChannelDriftCorrection.m` via `HyperStackReg`, `MultiChannelChromaticAberration.m`), the denoising with the Fiji `PureDenoise` plugin, the cell-cycle classification with Cellpose3 together with [BacImgLab/DeCNN](https://github.com/BacImgLab/DeCNN), and the septum profiling (`DemoDR_stage35.m`, `DemoDR_stage2.m`), which requires the S0 and S1 septum lines to be drawn by hand in a figure; `runDemoColoc` therefore starts from the published example data 4 and 6. For the FLIM workflow these are the LAS X FLIM/FCS export, which runs in the acquisition software, and the Bernsen background masking (`bersenThtest.ijm`), which runs in Fiji; `runDemoLifetime` therefore starts from the exported images (example data 1) and from the published mask (example data 3). For the SMT workflow these are the ROI cropping (`SMTdataPrepare.m`), the ThunderSTORM localisation, the chromatic-aberration correction (`CACorrectionofSMTsplitter.m`), the Cellpose3 segmentation, the trajectory linking (`spotsLinking.m`), the interactive trajectory segmentation in the `RefineTraceSegDr` app and the state classification (`statesClassifyDr.m`, which opens a save dialog and waits for a button press); `runDemoSMT` therefore starts from the published example data 4.

### Expected output

After `runDemoColoc`, `demo/Colocalization/demo_output/` contains:

    demo/Colocalization/demo_output/
    ├── S1_peak_distance/
    │   ├── WmemData.mat                                 mui_poi, mu1_mem, ParameterAll
    │   └── Wmem.tif                                     80 pages, 560 x 420, ≈ 55 MB
    └── PCC/
        ├── Merged_PCC_Values.csv                        all cells: FileName, C2C3, C2C4, C3C4, meanC2..4
        ├── stage2_PCC_Values.csv                        the same table per cell-cycle stage
        ├── stage3_PCC_Values.csv
        ├── stage4_PCC_Values.csv
        └── stage5_PCC_Values.csv

After `runDemoSMT`, `demo/SMT/demo_output/` contains:

    demo/SMT/demo_output/
    └── velocity_distribution/
        ├── FtsW-S1-Cef.mat                              the complete Result struct of the module script
        ├── FtsW_speed_CDF.png                           CDF of the speeds with both fits and residuals
        ├── FtsW_speed_PDF.png                           speed histogram with the two reconstructed PDFs
        └── SMT_velocity_summary.csv                     the fitted parameters in a flat two-column table

After `runDemoLifetime`, `demo/Lifetime/demo_output/` contains:

    demo/Lifetime/demo_output/
    ├── stacking/
    │   ├── WT50uMintensityStack.tif                     6 pages, 1024 x 1024, uint16
    │   ├── WT50uMlifetimeStack.tif                      6 pages, 1024 x 1024, uint16
    │   └── Lifetime_stacking_summary.csv                page counts and reference comparison
    └── quantification/
        ├── WT50uMintensityStack-filter.tif              background set to zero, 6 pages
        ├── WT50uMlifetimeStack-filter.tif               background set to zero, 6 pages
        ├── LifetimeResults.mat                          hLT, Result, LT_all, LT_mean, In_mean, In_all
        ├── Lifetime_histogram.png                       the lifetime histogram drawn by the script
        ├── Lifetime_summary.csv                         the numbers printed below
        └── Lifetime_per_page.csv                        pixel count and means per page

After `runDemoSMLM`, `demo/SMLM/demo_output/width/` contains:

    demo/SMLM/demo_output/width/
    ├── SMLM_width_summary.csv                all 116 profiles: FileName, HalfPeakWidth,
    │                                         Stage, HalfPeakWidth_nm
    ├── SMLM_width_stats.csv                  per stage: nProfiles, nNaN, nWidths, mean,
    │                                         median, sd, min, max (µm)
    ├── WTS0/
    │   ├── AllWidths.csv                     the summary as written by the module script
    │   ├── *_HalfPeakWidth.png               one annotated FWHM figure per S0 septum (51)
    │   └── work/                             the isolated working copy used for the run
    └── WTS1/
        ├── AllWidths.csv                     the summary as written by the module script
        ├── *_HalfPeakWidth.png               one annotated FWHM figure per S1 septum (65)
        └── work/                             the isolated working copy used for the run

All four demos also leave the scratch working copies that the wrappers assemble (`*/demo_output/*/work`), and `demo_output/` is not tracked by git.

Expected values for the multicolor colocalization demo (MATLAB R2024a):

* Step 1: `columns fitted : 80`, none of the 80 fits returns `NaN`. The protein ridge (channel 4, mNeonGreen) sits at 31.6 px on average (27.2–36.2 px) and the membrane ridge (channel 2, Potomac red) at 24.6 px (22.5–28.1 px), i.e. a mean ridge separation of 7.0 px (median 7.2 px, range 2.7–11.9 px). `Wmem.tif` has 80 pages of 560 × 420 px.
* Step 2: `cells analysed : 8` over the four cell-cycle stages; mean PCC(ch2–ch3) = 0.613, PCC(ch2–ch4) = 0.389 and PCC(ch3–ch4) = 0.909. The per-cell values are listed in `Merged_PCC_Values.csv`; PCC(ch3–ch4) stays within 0.862–0.979, while PCC(ch2–ch4) spreads over 0.172–0.561.
* The console lists the two steps and ends with `=== Colocalization demo finished in ≈ 60 s ===`.

Expected values for the single-molecule-tracking demo (MATLAB R2024a):

* 1405 segments in the table, of which 238 are directional and 877 stationary (location flag 3); 226 directional segments are faster than the 1 nm/s threshold and enter the CDF.
* Single log-normal fit: mean speed V = 14.40 nm/s (P = 1.000, µ = 2.200, σ = 0.967), bootstrap SEM 1.04 nm/s.
* Double log-normal fit: V1 = 2.69 nm/s (18.0 % of the population, σ = 0.31) and V2 = 15.74 nm/s (82.0 %, σ = 0.76), bootstrap SEM 4.38 and 1.49 nm/s; the two-population model is the one reported for FtsW.
* Maximum absolute CDF residual: 0.045 (single fit) and 0.025 (double fit).
* The console ends with `=== SMT demo finished in ≈ 20 s ===`.

Expected values for the FLIM demo (MATLAB R2024a):

* Step 1: `images to stack : 12`; both stacks have 6 pages of 1024 × 1024 pixel in uint16 and are **identical** to the published stacks of example data 2 (`6 of 6` pages for the intensity stack and for the lifetime stack).
* Step 2: `images analysed : 6` with **195 213 signal pixels** in total (19 908, 39 595, 19 534, 35 512, 27 848 and 52 816 per page); none of the pixels has a `NaN` lifetime. The lifetime is 1.269 ± 0.244 ns (range 0.690–2.143 ns) and the intensity 506.5 ± 119.7 a.u. (range 19–1037). The lifetime histogram has 110 bins between 0.30 and 3.60 ns and peaks at 0.0505 probability at 1.04 ns, with a second local maximum near 1.4 ns. Both background-free stacks are **identical** to the published images of example data 3 (`6 of 6` pages each).
* The console lists the two steps and ends with `=== FLIM demo finished in ≈ 10 s ===`.

Expected values for the 3D-dSTORM septum-width demo (MATLAB R2024a):

* **WTS0**: 51 profiles, none of them without an FWHM. Width 0.1097 ± 0.0206 µm (median 0.1068, range 0.0681–0.1775 µm) = **109.7 ± 20.6 nm**.
* **WTS1**: 65 profiles, none of them without an FWHM. Width 0.1062 ± 0.0309 µm (median 0.0989, range 0.0564–0.2030 µm) = **106.2 ± 30.9 nm**.
* The per-profile values are listed in `SMLM_width_summary.csv` and the per-stage statistics in `SMLM_width_stats.csv`; `WTS0/AllWidths.csv` and `WTS1/AllWidths.csv` are the tables exactly as the module script writes them.
* The console ends with `=== SMLM demo finished in ≈ 2 min ===`.

### Expected run time for the demo

`runDemoColoc`: ≈ 1 min on a normal desktop computer with MATLAB R2024a (measured 59.8 s and 63.1 s), of which about 55 s is step 1 — the 80 × 2 Gaussian fits plus one figure and one TIFF page per demograph column. Step 1 writes an uncompressed 80-page TIFF of ≈ 55 MB into `demo/Colocalization/demo_output/`; delete that folder afterwards if disk space matters. Step 2 takes 1–2 s.

`runDemoSMT`: ≈ 20 s on the same machine (measured 17.2 s), of which the 200 × 2 bootstrap fits take the larger share.

`runDemoLifetime`: ≈ 10 s on the same machine (measured 8.5 s and 9.1 s), of which the writing and re-reading of the six-page stacks takes the larger share. The demo leaves ≈ 100 MB of output and working copies in `demo/Lifetime/demo_output/`; delete that folder afterwards if disk space matters.

`runDemoSMLM`: ≈ 2 min on the same machine (measured 106 s and 129 s), almost all of it in the 116 figures the module script draws and saves. The demo leaves ≈ 11 MB of output and working copies in `demo/SMLM/demo_output/`; delete that folder afterwards if disk space matters.

The very first run after a cold start of MATLAB can take a few minutes because the graphics renderer is initialized then; subsequent runs take the values above.

### Scope of the demo

The four demos cover the numerical-analysis entry points of the multicolor colocalization, the FLIM, the single-molecule-tracking and the 3D-dSTORM septum-width workflow, i.e. the first points in the four pipelines that can be executed without human interaction. The upstream steps of all four workflows are interactive by design — they require manually drawn ROIs, per-cell visual confirmation, or external GUI tools (Fiji/ImageJ plugins, Cellpose, ThunderSTORM, the LAS X export, or the `RefineTraceSegDr` App Designer app) — and therefore cannot be replayed unattended.

| Module | Covered by the demo | Demo entry point | Interactive upstream steps not covered |
|---|---|---|---|
| Multicolor colocalization | Yes, the two unattended steps | `Colocalization/DemoS1AnalysisW.m` and `Colocalization/PCCcaclu.m` (via `demo/Colocalization/runDemoColoc.m`) | `MultiChannelDriftCorrection.m` (Miji + Fiji `HyperStackReg`), `PureDenoise`, `MultiChannelChromaticAberration.m`, the cell-cycle classification chain (Cellpose3 + [BacImgLab/DeCNN](https://github.com/BacImgLab/DeCNN)), and the manual S0/S1 point selection in `DemoDR_stage2.m` / `DemoDR_stage35.m` |
| Fluorescence lifetime imaging | Yes, the two scriptable sections | `Lifetime/lifeTcalculationbernsen.m`, Sections 1 and 3 (via `demo/Lifetime/runDemoLifetime.m`) | the LAS X FLIM/FCS export, which runs in the acquisition software, and the Bernsen masking plugin `bersenThtest.ijm`, which runs in Fiji |
| Single-molecule tracking | Yes, the speed-distribution fitting | `SMTanalysis/dataprocessSMTWCF.m` (via `demo/SMT/runDemoSMT.m`) | `SMTdataPrepare.m`, the ThunderSTORM localisation, `CACorrectionofSMTsplitter.m`, the Cellpose3 segmentation, `spotsLinking.m`, the `RefineTraceSegDr.mlapp` app and `statesClassifyDr.m` |
| 3D-dSTORM septum width | Yes, the FWHM calculation | `3D-dSTORM/dSTORMwidthcaclu.m` (via `demo/SMLM/runDemoSMLM.m`) | drawing the line ROIs across the septa in Fiji and exporting the intensity profile of each line |
| Electron tomography septum width | No — documentation package only | `Tomo/TomowidthCaclu.m` — needs an ImageJ polygon ROI plus two manual clicks to define the main axis; the example data and the walkthrough ship in `demo/Tomo/` | drawing the septum ROIs in Fiji, and the two axis clicks inside the script itself |

The electron-tomography module is the only one not covered: its numerical core consumes an input that no item of the published example data set provides — a manually drawn polygon ROI, plus two clicks inside the script itself — so an unattended replay would have to substitute that step rather than execute it; its published example data set (the slices, the drawn ROIs and the resulting width measurements) is small enough to ship, so it is documented and shipped in `demo/Tomo/` even though no step of it can be replayed. The four covered modules are precisely those whose interactive steps end at an item of the published example data set: the per-cell septum profiles and the smoothed demographs for the multicolor colocalization, the exported images and the published Bernsen mask for the FLIM workflow, the classified-trajectory table for the single-molecule tracking, and the exported intensity profiles of the drawn line ROIs for the 3D-dSTORM septum width.

For the multicolor colocalization module and for the FLIM module the interactive step sits in the middle of the workflow. In the colocalization case it is the septum profiling (example data 4), which requires the S0 and S1 septum lines to be drawn by hand and every step upstream of it runs inside Fiji, Cellpose or the DeCNN classifier; the two steps downstream of it consume the per-cell images and demograph matrices that this step produces, and are therefore exactly the two that `runDemoColoc` replays. In the FLIM case it is the Bernsen masking, which runs in Fiji between the stacking and the quantification; the published mask is example data 3, so the two sections on either side of it can be replayed from the published data. For the SMT module the interactive steps end with the state classification, whose output table is the published example data 4 that `runDemoSMT` consumes. For the 3D-dSTORM module the interactive step sits at the very beginning of the pipeline — the line-ROI drawing and the profile export in Fiji — and the published example data set is precisely its output, so the one step downstream of it, the FWHM calculation, is the one `runDemoSMLM` replays.

---

## Code functionality and documentation

The functionality of the scripts is described in this README, in the `Readme.md` file of each analysis module, and in the corresponding MATLAB scripts.

For each analysis module, the repository specifies the major processing steps, required external software, expected input data, and quantitative outputs.

Additional experimental details, imaging parameters, and data-analysis procedures are described in the associated manuscript Methods section.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

Copyright (c) 2026 BacImaging.

---

## Reference

Please cite the associated manuscript when using this code:

Class A PBPs reinforce the septal cell wall following initial synthesis by SEDS-bPBP pairs during bacterial cytokinesis.
*Author list, Journal, Year, DOI — to be completed once the manuscript is published.*

---

## Contact

For questions regarding the analysis code or its application to the datasets described in the associated manuscript, please contact the corresponding author.
