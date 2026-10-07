# 3D_dSTORM_Septum_Width_Analysis

This module provides tools for quantifying the relative septum width of *Deinococcus radiodurans* cells from 2D projections of 3D-dSTORM data reconstructed at ~200 nm thickness.

---

## 1. ROI Selection

- In Fiji, manually draw a line ROI across the septum of each cell on the 2D reconstructed dSTORM image.  
- The line should be perpendicular to the septum and cover the entire signal width.

---

## 2. Septum Width Calculation

- Use the Matlab script `dSTORMwidthcaclu.m` to compute the relative width of each septum.  
- The script calculates the **full width at half maximum (FWHM)** of the intensity profile along the ROI.  
- Output includes:
  - FWHM in pixels (can be converted to nanometers)  
  - Visualization of intensity profile and fitting curve (if applicable)

---

## Requirements

- Matlab 2020a or later  
- Fiji (ImageJ)  
- Input: 2D dSTORM images (e.g., reconstructed from ThunderSTORM or equivalent)  
- ROI format: ImageJ ROI file (line ROIs saved from Fiji)

---

## How to Use

1. Open the 2D dSTORM image in Fiji  
2. Draw line ROIs across septa and export them as `.roi` or `.zip` ROI files  
3. Run `dSTORMwidthcaclu.m` in Matlab to process the ROIs and calculate septum widths  
4. Review output plots and export numerical results


