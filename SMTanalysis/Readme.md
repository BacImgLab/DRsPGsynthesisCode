# Single_Molecule_Tracking

Single_Molecule_Tracking contains Matlab scripts and Fiji plugins for analyzing single-molecule tracking (SMT) data of *Deinococcus radiodurans*. The pipeline includes raw data preprocessing, trajectory linking, state classification, and drift removal.

---

## 1. Data Preprocessing

1. **ROI Selection and Cropping**  
   Identify a region of interest (ROI) with uniform illumination across the splitter using fluorescent bead images.  
   - Matlab Script: `SMTdataPrepare.m`  
   This script performs batch ROI cropping on all raw imaging data.

2. **Single-Molecule Localization**  
   Use Fiji to localize individual molecules from raw movies using ThunderSTORM.  
   - Fiji Plugin: `MacroThurderSTORMDr.ijm`  
   This macro automates single-molecule localization and outputs molecular position data.

3. **Chromatic Aberration Correction**  
   Correct chromatic offset between fluorescence and tracking channels.  
   - Matlab Script: `CACorrectionofSMTsplitter.m`

4. **Cell Segmentation Using Cellpose3**  
   Segment cell contours from bright-field images using Cellpose3, and save the masks as binary images.

5. **Trajectory Linking**  
   Link localizations between adjacent frames into continuous trajectories.  
   - Function: `Spotslink.m`

6. **Trajectory Extraction and MSD Analysis**  
   Extract trajectories within cell masks and calculate MSD for each trajectory.  
   - Functions and Scripts:  
     - `roitextRead.m`  
     - `traceInROI.m`  
     - `MSDsingle2D.m`  
     - `traceInROI.m` (main script for generating per-cell trajectory plots)

---

## 2. Trajectory Segmentation

- Use a custom Matlab App to segment each trajectory into unidirectional segments.  
- Matlab App: `RefineTraceSegDr`  
- Dependencies:
  - `MSDsingle2D.m`  
  - `linfitR.m`

---

## 3. Trajectory State Classification

1. **Classify Segment Motion States**  
   Determine motion type (directional, diffusive, stationary) for each segmented trajectory.  
   - Matlab Script: `statesClassifyDr.m`  
   - Dependencies:  
     `rcdfCal.m`, `mergeSegDrJ.m`, `segUpdate.m`, `segPRplot.m`,  
     `histLog_xy.m`, `MSDcaclulate_2d.m`, `kusumi_xy.m`, `kelsey2.m`, `colorCodeTracePlot.m`

2. **Velocity Distribution Fitting**  
   Fit velocity distributions (single or double log-normal) for directionally moving segments and compute mean velocity and population proportion.  
   - Matlab Script: `dataprocessSMTWCF.m`  
   - Dependencies:  
     `logn1cdf.m`, `logn2cdf.m`, `CDF_logCalc.m`

---

## 4. Visualization of Segmented Trajectories

- Visualize the position of segmented, directionally moving trajectories.  
- Matlab Script: `tracePlotafterseg.m`  
- Dependency: `colorCodeTracePlot.m`

---

## 5. Drift Region Removal

- Manually identify cells with apparent internal drift (2–3 pixels) using ROI tools in Fiji.  
- Read and remove those regions from analysis.  
  - Matlab Scripts:  
    - `ReadImageJROI.m`  
    - `RemoveDrift.m`

---

## Requirements

- Matlab 2020a or later  
- Fiji with ThunderSTORM and ROI Manager  
- Cellpose3 for cell segmentation  
- Input: Single-molecule movies, bright-field images, calibration data

---

## How to Use

1. Crop ROIs using `SMTdataPrepare.m`  
2. Run `MacroThurderSTORMDr.ijm` in Fiji to localize molecules  
3. Correct channel offset with `CACorrectionofSMTsplitter.m`  
4. Segment cells with Cellpose3 and load masks  
5. Link trajectories with `Spotslink.m`  
6. Use `traceInROI.m` and `MSDsingle2D.m` to extract trajectories per cell  
7. Segment trajectories using `RefineTraceSegDr`  
8. Classify motion states with `statesClassifyDr.m`  
9. Fit speed distributions using `dataprocessSMTWCF.m`  
10. Visualize segmented trajectories with `tracePlotafterseg.m`  
11. Remove drifting cells using `ReadImageJROI.m` and `RemoveDrift.m`


