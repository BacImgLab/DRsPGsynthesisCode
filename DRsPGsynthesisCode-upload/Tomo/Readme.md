# Tomography_Septum_Width_Analysis

This module provides tools for quantifying the relative septum width of *Deinococcus radiodurans* cells from 2D tomography images.

---

## 1. ROI Selection

- Open the 2D tomography image in Fiji.  
- Manually draw line ROIs across the septum of each cell.  


---

## 2. Septum Width Calculation

- Use the Matlab script `TomowidthCaclu.m` to calculate the relative septum width 
- The output includes:
  - Width measurement 


---

## Requirements

- Matlab 2020a or later  
- Fiji (ImageJ)  
- Input: 2D tomography images  
- ROI files exported from Fiji (line ROIs)

---

## How to Use

1. Open tomography images in Fiji  
2. Draw line ROIs across septa and save the ROIs  
3. Run `TomowidthCaclu.m` in Matlab to process ROIs and calculate septum widths  
4. Inspect plots and export numerical data

---
