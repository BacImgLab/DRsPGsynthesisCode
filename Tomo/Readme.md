# Tomography_Septum_Width_Analysis

This module provides tools for quantifying the relative septum width of *Deinococcus radiodurans* cells from 2D tomography images.

The example used throughout this document is a WT sample with one cell in stage S0 and one in stage S1. A documentation package with the published example data set of this module is available in [`demo/Tomo/`](../demo/Tomo/Readme.md).

---

## 1. ROI Selection

- Open the 2D tomography image in Fiji.
- Manually draw line ROIs across the septum of each cell.
- Save each ROI as an ImageJ `.roi` file; its name becomes the base name of every output of step 2.

---

## 2. Septum Width Calculation

- Use the Matlab script `TomowidthCaclu.m` to calculate the relative septum width
- The output includes:
  - Width measurement

The script processes an ImageJ ROI file, asks for two clicks that define the main axis of the septum, then measures the local width perpendicular to that axis:

1. `uigetfile` selects the `.roi` file; `ReadImageJROI.m` returns the polygon vertices, which are converted from pixels to nanometres (`pixelS = 0.86` nm/pixel) and closed into a ring (`CoordXY`).
2. `ginput(2)` collects the two endpoints of the septum main axis.
3. The axis is sampled every `interD = 1` nm; at each sample a perpendicular slice (half-length `M = 200` nm) is intersected with the polygon (`polyxpoly`), and the midpoint of each two-point intersection becomes a point of the middle line.
4. Local tangents of the middle line are estimated by finite differences, a fresh perpendicular slice is cast through every midpoint, and the local width is the distance between its two intersections.
5. The results are saved next to the ROI file:
   - `<ROI base name>.mat` — `Res` = `[distance along the middle line, local width, midpoint_x, midpoint_y]` in nm, and `CoordXY` = the closed polygon in nm;
   - `<ROI base name>.jpg` — the annotated figure (polygon, slices, middle line, width segments).

Because both ends of the workflow are interactive — the ROI is drawn by hand in Fiji and the main axis is picked by hand in the MATLAB figure window — this module has no step that runs unattended, and there is no automated demo for it. The example data and the full walkthrough are in [`demo/Tomo/`](../demo/Tomo/Readme.md).

---

## Requirements

- Matlab 2020a or later, with the **Mapping Toolbox** (`polyxpoly`)
- `ReadImageJROI.m` (shipped in this repository in `SMTanalysis/`)
- Fiji (ImageJ)
- Input: 2D tomography images
- ROI files exported from Fiji (polygon ROIs)

---

## How to Use

1. Open tomography images in Fiji
2. Draw polygon ROIs across septa and save the ROIs
3. Run `TomowidthCaclu.m` in Matlab to process ROIs and calculate septum widths
4. Inspect plots and export numerical data

---
