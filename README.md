# DRsPGsynthesisCode

DRsPGsynthesisCode contains Matlab scripts and Fiji plugins for analyzing *Deinococcus radiodurans* imaging data, including multicolor colocalization, fluorescence lifetime imaging, single-molecule tracking (SMT), 3D-dSTORM septum width analysis, and tomography septum width analysis.

---

###   Multicolor Colocalization Imaging

Includes preprocessing steps such as drift correction (via `MultiChannelDriftCorrection.m` with Fiji plugin Hyperstackreg), chromatic aberration correction (`MultiChannelChromaticAberration.m`, depends on `imreg2Dr`), and denoising (Fiji plugin Puredenoise). Demography analysis scripts combine images, record septum features, and analyze signal trends with various Matlab scripts like `combineSegfiles4C.m`, `DemoDR_stage35.m`, `Beforedemoprocess.m`, and `demoSmooth.m`. Pearson correlation coefficients between channels are calculated by `PCCcaclu.m`.

---

###   Fluorescence Lifetime Imaging (FLIM)

Export intensity and lifetime images from Leica LAS X FLIM/FCS software, then use `lifeTcalculationbernsen.m` to stack images and process background subtraction (via Fiji plugin bersenThtest). The script computes pixel-wise lifetime and intensity values.

---

###   Single Molecule Tracking (SMT)

Preprocessing includes ROI cropping (`SMTdataPrepare.m`), localization with ThunderSTORM in Fiji, chromatic correction (`CACorrectionofSMTsplitter.m`), cell segmentation with Cellpose3, and trajectory linking (`Spotslink.m`). Trajectory segmentation and state classification are performed with `RefineTraceSegDr` and `statesClassifyDr.m`, with velocity distribution fitting via `dataprocessSMTWCF.m`. Visualization and drift region removal are included with scripts like `tracePlotafterseg.m` and `RemoveDrift.m`.

---

###   3D-dSTORM Septum Width Analysis

Manually draw septum ROIs on 2D projections of 3D-dSTORM images (~200 nm thick) in Fiji. Calculate septum width using `dSTORMwidthcaclu.m`, which measures full width at half maximum (FWHM) of intensity profiles.

---

###   Tomography Septum Width Analysis

Draw septum ROIs on 2D tomography images in Fiji. Use `TomowidthCaclu.m` to quantify septum width from the intensity profile along the ROI.

---

## Requirements

- Matlab 2020a or newer  
- Fiji (ImageJ) with required plugins  
- Cellpose3 for cell segmentation  
- Input data as specified per module

---

## How to Use

Please refer to individual module scripts and included documentation for detailed step-by-step usage instructions.

---

## Reference

Please cite this paper when using DRsPGsynthesisCode:  
XXX
