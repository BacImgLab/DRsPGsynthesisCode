# DRsPGsynthesisCode

`DRsPGsynthesisCode` is a collection of MATLAB scripts and Fiji/ImageJ-based workflows for quantitative analysis of *Deinococcus radiodurans* imaging data.

The repository provides analysis workflows for:

* Multicolor fluorescence colocalization and demography analysis
* Fluorescence lifetime imaging microscopy (FLIM)
* Single-molecule tracking (SMT)
* 3D-dSTORM septum width analysis
* Electron tomography (ET) septum width analysis

The repository contains analysis scripts and workflows that use Fiji/ImageJ and other external software for image preprocessing, segmentation, registration, localization, trajectory analysis, and quantitative measurements.

A small **simulated** demo dataset and a one-command demo workflow are provided in `demo/`; they allow the code to be tested without any experimental data, without Fiji/ImageJ and without any MATLAB toolbox (see [Demo](#demo)).

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

A small simulated demo dataset and a ready-to-run demo script covering three of these modules (multicolor colocalization, SMT and 3D-dSTORM septum width) are provided in `demo/` (see [Demo](#demo)).

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

Running the simulated demo described in the [Demo](#demo) section requires only MATLAB and takes less than a minute (≈ 20–30 s).

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

### Electron tomography

**Input:**

* 2D electron tomography images
* Manually defined septum ROIs

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

The repository includes a small simulated demo dataset and a single command that runs three representative workflows on it: the septum-width analysis of 3D-dSTORM intensity profiles, a trajectory analysis of simulated single-molecule tracking data, and the Pearson-correlation calculation of a simulated multicolor colocalization dataset. See the [Demo](#demo) section below.

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
| Fig. 2b–d | Septal enrichment, demographs, displacement | `Colocalization/S0S1select.m`, `Colocalization/lineProfile.m`, `Colocalization/demoSmooth.m` |
| Fig. 3, Fig. 6a–d; Supplementary Fig. 9, 10f–g | FLIM lifetime and intensity | `Lifetime/lifeTcalculationbernsen.m` |
| Fig. 4f–g; Supplementary Table 5 | 3D-dSTORM / 3D-SMLM septum width | `3D-dSTORM/dSTORMwidthcaclu.m` |
| Fig. 5a–g; Supplementary Fig. 15–17 | SMT: localization, linking, state classification, MSD, speed distributions | `SMTanalysis/SMTdataPrepare.m` → ThunderSTORM → `SMTanalysis/spotsLinking.m` → `SMTanalysis/RefineTraceSegDr.mlapp` → `SMTanalysis/statesClassifyDr.m` → `SMTanalysis/dataprocessSMTWCF.m` and `SMTanalysis/MSDsingle2D.m` / `SMTanalysis/CDF_logCalc.m` (MSD, step-length CDF and velocity fits: `SMTanalysis/linfitR.m`) |
| Fig. 5d; Supplementary Fig. 16d–g | 2D-projection correction of FtsW speed | `SMTanalysis/unwrapTraj.m`, `SMTanalysis/DrunwrapX.m`, `SMTanalysis/DrunwrapY.m` |

Raw imaging datasets (single-molecule localization movies, FLIM photon data and electron-tomography tilt series) are not included in this repository because of their very large size; they are available from the corresponding author upon reasonable request. The simulated demo dataset in `demo/` allows the code to be tested without them.

---

## Demo

A small **simulated** demo dataset is provided in `demo/demo_data/`, together with a single command that runs three representative analysis workflows on it. Running the demo requires only base MATLAB: no Fiji/ImageJ, no additional plugin, no MATLAB toolbox and no experimental data are needed.

### Demo dataset

| File | Content | Used by |
|---|---|---|
| `demo/demo_data/dSTORM_profile_01.csv` … `_05.csv` | Five simulated septum intensity profiles of a 2D projection of 3D-dSTORM data. Column 1: distance along the line ROI (nm); column 2: gray value. | `3D-dSTORM/dSTORMwidthcaclu.m` |
| `demo/demo_data/SMT_tracks_demo.csv` | Simulated 2D single-molecule tracking dataset: 12 molecules × 200 frames, 0.16 µm per pixel, 110 ms per frame. Columns: track, frame, x, y, intensity. | `demo/runDemoSMT.m` (uses `SMTanalysis/MSDsingle2D.m`, `SMTanalysis/CDF_logCalc.m` and `SMTanalysis/linfitR.m`) |
| `demo/demo_data/PCC/stage2 … stage5/Ch4_DR_00X.tif` | Eight simulated 4-channel (4-page) 64 × 64 uint16 TIFF stacks, two fields of view per stage folder. Channels 2, 3 and 4 share a common spatial structure with decreasing correlation. | `Colocalization/PCCcaclu.m` via `demo/runDemoPcc.m` |

All files are simulated and contain no experimental measurements. They are reproducible: the septum profiles and the colocalization stacks are generated deterministically, and the tracking dataset is regenerated by `demo/makeDemoData.m` with a fixed random seed (`rng(2026,'twister')`).

### Instructions to run on the demo data

In MATLAB:

    >> cd demo
    >> runDemo

or, from a terminal:

    matlab -batch "cd demo; runDemo"

`runDemo.m` executes four steps:

1. `makeDemoData.m` regenerates the simulated datasets in `demo/demo_data/` (septum profiles, tracking table and multicolor stacks).
2. `runDemoDstorm.m` copies the profile CSVs into `demo/demo_output/dSTORM_profiles/` and runs `3D-dSTORM/dSTORMwidthcaclu.m` on them.
3. `runDemoSMT.m` analyses `SMT_tracks_demo.csv` with `MSDsingle2D.m`, `CDF_logCalc.m` and `linfitR.m`.
4. `runDemoPcc.m` assembles an isolated working copy of the stage folders and runs `Colocalization/PCCcaclu.m` on the simulated 4-channel stacks.

### Expected output

After the run, `demo/demo_output/` contains:

    demo/demo_output/
    ├── dSTORM_profiles/
    │   ├── Plots/AllWidths.csv                        FWHM of each profile (1 header + 5 rows)
    │   └── Plots/dSTORM_profile_0X_HalfPeakWidth.png  annotated FWHM plot per profile
    ├── SMT_MSD_demo.csv                               lag time (s), mean MSD (µm²), SEM (µm²)
    ├── SMT_MSD_demo.png                               log-log MSD curve
    ├── SMT_stepCDF_demo.csv                           step length (µm), cumulative probability (101 bins)
    ├── SMT_stepCDF_demo.png                           step-length CDF plot
    ├── SMT_velocity_demo.csv                          per-trajectory linear-fit velocity (12 rows)
    └── PCC/
        ├── Merged_PCC_Values.csv                      pairwise PCC and mean intensities (8 rows)
        └── work/                                      isolated working copy used for the run

All results are written into `demo/demo_output/`, which is not tracked by git; the folder also holds a scratch working copy (`demo/demo_output/PCC/work/`) that `runDemoPcc.m` assembles so that `PCCcaclu.m` can run without modifying the repository tree.

Expected values for the delivered data (MATLAB R2024a):

* `AllWidths.csv`: FWHM ≈ 23.7, 28.4, 32.7, 37.3 and 42.1 nm for the five simulated septa (the simulated half-widths σ are 10, 12, 14, 16 and 18 nm, i.e. FWHM ≈ 2.355·σ).
* `SMT_MSD_demo.csv`: 50 lag times from 0.11 s to 5.50 s; the mean MSD increases approximately linearly from ≈ 0.015 µm² (0.11 s) to ≈ 0.59 µm² (5.50 s), corresponding to an effective diffusion coefficient of ≈ 0.027 µm²/s (the simulated values were D = 0.005, 0.02 and 0.05 µm²/s plus 30 nm localization uncertainty).
* `SMT_stepCDF_demo.csv`: 50 % of the simulated steps are shorter than ≈ 0.09 µm and 90 % shorter than ≈ 0.20 µm.
* `SMT_velocity_demo.csv`: 12 rows with a mean linear-fit speed of ≈ 0.07 µm/s. (For demonstration the linear fit is applied to whole trajectories here; in the analysis of the experimental data it is applied to trajectory segments after state segmentation, see `statesClassifyDr.m` / `RefineTraceSegDr.mlapp`.)
* `PCC/Merged_PCC_Values.csv`: 8 rows (2 fields of view × 4 stage folders) with the three pairwise Pearson correlation coefficients and the mean channel intensities. On the simulated stacks the mean values are PCC(ch2–ch3) ≈ 0.966, PCC(ch2–ch4) ≈ 0.771 and PCC(ch3–ch4) ≈ 0.745, reproducing the decreasing correlation that was used to generate the data.
* The console lists the four steps and ends with `Demo finished in ≈ 25 s`.

### Expected run time for the demo

≈ 20–30 s on a normal desktop computer with MATLAB R2024a, including writing all figures (measured 19.5 s, 23.1 s and 25.2 s on Windows 10/11 with MATLAB R2024a). All four steps together stay well below one minute. The very first run after a cold start of MATLAB can take a few minutes because the graphics renderer is initialized then; subsequent runs take ≈ 25 s.

### Scope of the demo

The demo covers the numerical-analysis entry point of each module, i.e. the first
point in every pipeline that can be executed without human interaction. The
upstream steps of all modules are interactive by design — they require manually
drawn ROIs, per-cell visual confirmation, or external GUI tools (Fiji/ImageJ
plugins, Cellpose, ThunderSTORM, the LAS X export, or the `RefineTraceSegDr` App
Designer app) — and therefore cannot be replayed unattended.

| Module | Covered by the demo | Demo entry point | Interactive upstream steps not covered |
|---|---|---|---|
| Multicolor colocalization | Yes | `Colocalization/PCCcaclu.m` (via `runDemoPcc.m`) | `MultiChannelDriftCorrection.m` (Miji + Fiji `HyperStackReg`), `PureDenoise`, `MultiChannelChromaticAberration.m`, and the S0/S1 manual point selection in `DemoDR_stage2.m` / `DemoDR_stage35.m` |
| Single-molecule tracking | Yes | `SMTanalysis/MSDsingle2D.m`, `CDF_logCalc.m`, `linfitR.m` (via `runDemoSMT.m`) | `SMTdataPrepare.m`, ThunderSTORM localisation, Cellpose3 segmentation, `CACorrectionofSMTsplitter.m`, `Spotslink.m`, the `RefineTraceSegDr.mlapp` app and `statesClassifyDr.m` |
| 3D-dSTORM septum width | Yes | `3D-dSTORM/dSTORMwidthcaclu.m` (via `runDemoDstorm.m`) | drawing the line ROIs across the septa in Fiji |
| Electron tomography septum width | No | `Tomo/TomowidthCaclu.m` — needs an ImageJ line ROI plus two manual clicks to define the main axis | the same Fiji line-ROI drawing step as 3D-dSTORM |
| FLIM | No | `Lifetime/lifeTcalculationbernsen.m` — needs an intensity/lifetime stack and a Bernsen background mask | the LAS X FLIM/FCS export and the Bernsen masking plugin, both run outside MATLAB |

The two modules that are not covered share the same reason: their numerical core
consumes an input that is produced by an interactive step — a manually drawn line
ROI in the case of tomography, a Fiji-generated background mask in the case of
FLIM — so an unattended replay would have to substitute that step rather than
execute it. The three covered modules are precisely those whose numerical core
takes a plain, scriptable input: an intensity-profile CSV, a linked-trajectory
table, and a 4-channel image stack.

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
