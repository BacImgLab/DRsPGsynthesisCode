# 3D_dSTORM_Septum_Width_Analysis

This module provides tools for quantifying the relative septum width of *Deinococcus radiodurans* cells from 2D projections of 3D-dSTORM data reconstructed at ~200 nm thickness.

The example used throughout this document is a WT sample with one folder of septum profiles per cell-cycle stage (S0 and S1). An automated demo of the FWHM step on the published example data set is available in [`demo/SMLM/`](../demo/SMLM/Readme.md).

---

## 1. ROI Selection

- In Fiji, manually draw a line ROI across the septum of each cell on the 2D reconstructed dSTORM image.
- The line should be perpendicular to the septum and cover the entire signal width.
- Export the intensity profile along the line: this writes one `X,Y` profile per septum (column `X` = distance along the line ROI, column `Y` = grey value), one folder per cell-cycle stage.

---

## 2. Septum Width Calculation

- Use the Matlab script `dSTORMwidthcaclu.m` to compute the relative width of each septum.
- The script batches over every `*.csv` profile in the folder given by `sourceFolder` and computes the **full width at half maximum (FWHM)** of each intensity profile:
  1. `readmatrix` reads the file and drops the leading `X,Y` header row (which arrives as `NaN,NaN`).
  2. The half height is `(max(Y) + min(Y)) / 2`.
  3. Every crossing of the half height is located by linear interpolation between two consecutive samples.
  4. The width is the distance between the first and the last crossing; profiles with fewer than two crossings get `NaN`.
- Output includes:
  - `AllWidths.csv` — one `FileName, HalfPeakWidth` row per profile, in the unit of the `X` column
  - One annotated figure per profile (`<profile name>_HalfPeakWidth.png`): profile, half height, both crossings and the width
- Both outputs are written into a `Plots/` subfolder **inside the profile folder**.
- `sourceFolder` may be defined in the workspace before the script is run (this is how the demo calls it); it starts with `clearvars -except sourceFolder` and falls back to a hard-coded path otherwise.
- The shipped profiles are sampled every 0.005 µm, i.e. X is in micrometre and `HalfPeakWidth` × 1000 gives nanometres.

---

## Requirements

- Matlab 2020a or later; base MATLAB only (`readmatrix`, `writecell`, `plot`, `saveas`)
- Fiji (ImageJ)
- Input: 2D dSTORM images (e.g., reconstructed from ThunderSTORM or equivalent)
- ROI format: ImageJ ROI file (line ROIs saved from Fiji), exported as `X,Y` intensity profiles

---

## How to Use

1. Open the 2D dSTORM image in Fiji
2. Draw line ROIs across septa and export the intensity profile along each line
3. Point `sourceFolder` at the folder of profiles (or edit the default path in the script) and run `dSTORMwidthcaclu.m`
4. Review `Plots/AllWidths.csv` and the per-profile figures

---
