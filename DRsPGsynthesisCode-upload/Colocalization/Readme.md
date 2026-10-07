# Multicolor_Colocalization

Multicolor_Colocalization contains Matlab scripts and Fiji plugins for processing multicolor fluorescence imaging data of *Deinococcus radiodurans*, enabling the analysis of spatial protein distribution, septum-associated signal patterns, and colocalization across channels.

---

## 1. Data Preprocessing

1. **Drift Correction and Z-Projection**  
   Multichannel Z-stack images are automatically corrected for lateral drift and Z-projected using the provided Matlab script.  
   The script `MultiChannelDriftCorrection.m` uses the Miji interface to launch Fiji and run the `HyperStackReg` plugin directly within Matlab.  
   The output is one Z-projected 2D image per channel (averaged along Z).  
   - Matlab Script: `MultiChannelDriftCorrection.m`  
   - Fiji Plugin (called internally): `HyperStackReg`

2. **Chromatic Aberration Correction**  
   Channel misalignment caused by chromatic aberration is corrected using calibration images (e.g., multicolor beads).  
   - Matlab Script: `MultiChannelChromaticAberration.m`  
   - Dependency: `imreg2Dr.m`

3. **Denoising**  
   Each channel image is denoised to improve signal-to-noise ratio prior to downstream analysis.  
   - Fiji Plugin: `PureDenoise` (run manually in Fiji)

---

## 2. Demograph Construction

1. **Channel Stacking**  
   Z-projected images from different channels are merged into a multichannel stack for further analysis.  
   - Matlab Script: `combineSegfiles4C.m`

2. **Septum Signal Profiling in Individual Cells**  
   Signal intensity profiles across septal regions (S1 and S0) are extracted from each segmented cell.  
   - Matlab Scripts:  
     - `DemoDR_stage35.m`  
     - `DemoDR_stage2.m`  
   - Dependencies:  
     - `bf2rotateRec.m`, `DR4cShow.m`, `lineProfile.m`, `S0S1select.m`, `S0select.m`

3. **Signal Sorting and Trend Visualization**  
   Septum signals across all cells are sorted by septum-specific metrics (e.g., maturity stage) to visualize changes over the cell cycle.  
   - Matlab Script: `Beforedemoprocess.m`

4. **Demograph Smoothing**  
   The generated demograph is smoothed to better reveal population-level trends.  
   - Matlab Script: `demoSmooth.m`  
   - Dependency: `average_colN.m`

---

## 3. S1 Septum Demograph Analysis

- Quantifies the distance between signal peaks in two fluorescence channels across S1 septa.  
- Matlab Script: `DemoS1AnalysisW.m`

---

## 4. Colocalization Analysis (PCC Calculation)

- Calculates Pearson correlation coefficients (PCC) between two channels for each individual cell.  
- Matlab Script: `PCCcaclu.m`

---

## Requirements

- Matlab 2020a or later  
- Fiji (ImageJ) with:
  - `HyperStackReg` plugin  
  - `PureDenoise` plugin  
- Miji interface (ImageJ-Matlab bridge) installed and configured  
- Input: Multichannel Z-stack fluorescence images (e.g., `.tif`) and brightfield images for cell segmentation

---

## How to Use

1. Run `MultiChannelDriftCorrection.m` in Matlab to automatically apply drift correction and Z-projection  
2. Denoise the output images using Fiji's `PureDenoise` plugin  
3. Use `MultiChannelChromaticAberration.m` to correct for chromatic aberration  
4. Stack the corrected images using `combineSegfiles4C.m`  
5. Extract and analyze septal profiles using `DemoDR_stage35.m` and `DemoDR_stage2.m`  
6. Sort and visualize signal progression using `Beforedemoprocess.m` and `demoSmooth.m`  
7. Quantify peak distances using `DemoS1AnalysisW.m`  
8. Calculate PCC for colocalization analysis using `PCCcaclu.m`



