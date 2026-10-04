# DRimage_processing

`DRimage_processing` is a collection of MATLAB scripts and Fiji/ImageJ-based workflows for quantitative analysis of *Deinococcus radiodurans* imaging data.

The repository provides analysis workflows for:

* Multicolor fluorescence colocalization and demography analysis
* Fluorescence lifetime imaging microscopy (FLIM)
* Single-molecule tracking (SMT)
* 3D-dSTORM septum width analysis
* Electron tomography (ET) septum width analysis

The repository contains analysis scripts and workflows that use Fiji/ImageJ and other external software for image preprocessing, segmentation, registration, localization, trajectory analysis, and quantitative measurements.

---

## Repository structure

The main analysis modules are organized according to the type of imaging data:

1. Multicolor Colocalization Imaging
2. Fluorescence Lifetime Imaging (FLIM)
3. Single-Molecule Tracking (SMT)
4. 3D-dSTORM Septum Width Analysis
5. Electron Tomography Septum Width Analysis

---

## 1. Multicolor Colocalization Imaging

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
  Performs chromatic aberration correction using `imreg2Dr`.

* `combineSegfiles4C.m`
  Combines segmentation and quantitative data for multichannel analysis.

* `Beforedemoprocess.m`
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
* `imreg2Dr` for chromatic aberration correction

---

## 2. Fluorescence Lifetime Imaging (FLIM)

This module is used to process fluorescence lifetime and intensity images exported from Leica LAS X FLIM/FCS software.

### Workflow

1. Acquire FLIM data using a Leica STELLARIS 8 confocal microscope.
2. Export intensity and lifetime images from Leica LAS X FLIM/FCS software.
3. Perform image stacking and preprocessing using MATLAB.
4. Perform background subtraction using the Fiji/ImageJ `bersenThtest` plugin.
5. Calculate pixel-wise fluorescence lifetime and intensity values.
6. Export quantitative measurements for downstream statistical analysis.

### Main script

* `lifeTcalculationbernsen.m`
  Stacks exported FLIM images, performs image processing and background correction, and calculates pixel-wise fluorescence lifetime and intensity values.

### Input

Intensity and lifetime image files exported from Leica LAS X FLIM/FCS software.

### Output

Processed fluorescence lifetime and intensity measurements for downstream quantitative analysis.

---

## 3. Single-Molecule Tracking (SMT)

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

* `Spotslink.m`
  Links localized single-molecule positions into trajectories.

* `RefineTraceSegDr`
  Performs trajectory segmentation and refinement.

* `statesClassifyDr.m`
  Classifies molecular trajectories into movement states.

* `dataprocessSMTWCF.m`
  Performs velocity distribution fitting and related quantitative analysis.

* `tracePlotafterseg.m`
  Generates trajectory visualizations after trajectory segmentation.

* `RemoveDrift.m`
  Removes regions or trajectories affected by sample drift.

### External software

The SMT workflow uses:

* Fiji/ImageJ
* ThunderSTORM for single-molecule localization
* Cellpose3 for cell segmentation

---

## 4. 3D-dSTORM Septum Width Analysis

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
* `imreg2Dr`
* `bersenThtest`

### Operating system

* Operating system: Windows 10 / 11 (64-bit)

### Tested environment

The analysis workflows were tested using **MATLAB R2024a on Windows 11 (64-bit)**, together with ImageJ 1.54f, Cellpose 3 and ThunderSTORM (dev-2016-09-04-b1).

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
* `imreg2Dr`
* `bersenThtest`

Not all tools are required for every analysis module.

### 7. MATLAB scripts

No compilation is required.

Open the MATLAB script corresponding to the desired analysis module and specify the input data path and analysis parameters before running the script.

---

## Typical installation time

If MATLAB, Fiji/ImageJ, Cellpose3, and ThunderSTORM are already installed, obtaining the repository and preparing the MATLAB scripts typically takes only a few minutes.

Installation time for external software and plugins depends on the user's computer, network connection, and existing software environment.

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

Module-specific input requirements and processing steps are described in the corresponding MATLAB scripts and documentation.

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

To reproduce an analysis:

1. Obtain the corresponding imaging dataset.
2. Follow the preprocessing workflow described for the relevant analysis module.
3. Run the corresponding MATLAB script.
4. Use the analysis parameters specified in the script and/or manuscript Methods.
5. Apply the same data-selection and segmentation criteria described in the manuscript.

The repository contains the analysis scripts used for quantitative image processing. Raw imaging datasets are not currently included in this repository.

---

## Demo dataset

A separate demo dataset is not currently provided with this repository.

Users should therefore use their own imaging datasets prepared according to the input requirements described above.

---

## Code functionality and documentation

The functionality of the scripts is described in this README and in the corresponding MATLAB scripts.

For each analysis module, the repository specifies the major processing steps, required external software, expected input data, and quantitative outputs.

Additional experimental details, imaging parameters, and data-analysis procedures are described in the associated manuscript Methods section.

---

## Reference

Please cite the associated manuscript when using `DRimage_processing`:

**XXX**

---

## Contact

For questions regarding the analysis code or its application to the datasets described in the associated manuscript, please contact the corresponding author.

