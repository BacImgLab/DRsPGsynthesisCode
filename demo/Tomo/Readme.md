# Electron tomography septum width analysis (Tomo)

`Tomo/` contains the MATLAB script of the electron-tomography workflow for *Deinococcus radiodurans*.
The workflow measures the **relative width of the division septum** along a manually drawn ROI on a
2D tomographic slice.

`demo/Tomo/` is a **documentation package** of that module: the module script, the reader function
it depends on, the example data set in `demo_data/` and this Readme all live in this single folder,
so the folder can be downloaded and read on its own.

**This package has no runnable entry point, on purpose.** Both ends of the workflow need a human:

* the septum ROI is drawn **by hand in Fiji**, and
* `TomowidthCaclu.m` itself asks for two clicks in a MATLAB figure window to define the main axis
  (`ginput(2)`), after a file dialog (`uigetfile`).

An unattended replay would have to substitute the very step this module exists for, so the package
ships the published example data set and the step-by-step instructions instead. The two
*upstream-preparation* and *downstream-numerics* ends of the other workflows are covered by the
three runnable demo packages (`demo/Colocalization/`, `demo/Lifetime/`, `demo/SMT/`); this module
has no step that runs unattended (see the *Scope of the demo* section of the repository `README.md`).

---

## Workflow

| Step | What happens | Interactive? |
|---|---|---|
| 1 | [ROI selection in Fiji](#1-roi-selection-in-fiji) | **Yes — manual** |
| 2 | [Septum-width calculation in MATLAB](#2-septum-width-calculation-in-matlab) | **Yes — two clicks** |

### 1. ROI selection in Fiji

* Open the 2D tomographic slice in **Fiji/ImageJ**.
* Draw a closed **polygon ROI** around the septum of one dividing cell: the polygon encloses the
  whole septum disc between the two daughter cells.
* Save the ROI (`File > Save As > ROI...`, or the ROI Manager) — this writes an ImageJ `.roi` file
  whose name becomes the base name of every output of step 2.

→ produced: **example data 1** — the 2D tomographic slices and the manually drawn polygon ROIs.

### 2. Septum width calculation in MATLAB

* MATLAB script: `TomowidthCaclu.m`
* **Input**: one ImageJ `.roi` file (example data 1)
* **Output**: one `.mat` result file and one `.jpg` figure per ROI, written **next to the ROI file**
  (example data 2)

What the script does:

1. `uigetfile` asks for the `.roi` file; `ReadImageJROI.m` parses it and returns the polygon
   vertices, which the script converts from pixels to nanometres (`pixelS = 0.86` nm/pixel) and
   closes into a ring (`CoordXY`).
2. The polygon is plotted and **`ginput(2)` asks for two clicks** — the two endpoints of the main
   axis of the septum (roughly the ridge from one outer membrane to the other).
3. The main axis is sampled every `interD = 1` nm; at every sample a perpendicular slice of
   half-length `M = 200` nm is intersected with the polygon (`polyxpoly`). The midpoint of each
   two-point intersection is a point of the **middle line**.
4. Local tangents of the middle line are estimated by finite differences, a fresh perpendicular
   slice is cast through every midpoint, and the **local width** is the distance between the two
   intersections.
5. The results are saved as
   * `Res` — `[distance along the middle line, local width, midpoint_x, midpoint_y]`, in nm,
   * `CoordXY` — the closed polygon in nm,
   into `<ROI base name>.mat`, and the annotated figure into `<ROI base name>.jpg`.

The three parameters at the top of the script are the ones used for the published data:
`pixelS = 0.86` nm/pixel, `interD = 1` nm (sampling spacing along the axis), `M = 200` nm
(half-length of the slice lines, only needs to be long enough to leave the polygon on both sides).

---

## Example data

`demo_data/` holds the two items of the published example data set, under the English name of the
original item; the image, ROI and result files themselves are unchanged. 8 files, ≈ 1.5 MB in total.

| Folder | Content | Role in the workflow |
|---|---|---|
| `01_manual_septum_roi/` | Two 2D tomographic slices of WT cells (`Position_1_Aretomo_bin1-1-23Roi.jpg`, 963 × 1126; `Position_3_Aretomo_bin1-1-17Roi.jpg`, 933 × 1143; 8-bit RGB JPEG, 200 nm scale bar) and two Fiji polygon ROIs: `...S0Roi.roi` (54 vertices) and `...S1Roi_1.roi` (30 vertices). | step 1 (output) / step 2 (input) |
| `02_septum_width_calculation/` | The published script output for the two cells: the annotated figure (`...S0Roi.jpg`, `...S1Roi_2.jpg`, 1609 × 4300) and the result file (`...S0Roi.mat`, `...S1Roi_2.mat`). | step 2 (output) |

The published result files contain:

| Result file | `Res` | Width along the middle line |
|---|---|---|
| `...S0Roi.mat` | 715 × 4, no NaN | 27.9 ± 3.7 nm (range 22.7–51.7 nm) over 714 nm of axis |
| `...S1Roi_2.mat` | 125 × 5, no NaN | 27.7 ± 8.6 nm (range 10.3–76.4 nm) over 126 nm of axis |

### A note on the file names

The two example **inputs** and the two example **outputs** are not two matched pairs:

* `...S0Roi.mat` is the output of the shipped `...S0Roi.roi`: its `CoordXY` (55 × 2, first vertex
  864 / 1334 px at `pixelS = 0.86`) closes exactly the 54-vertex polygon of the shipped ROI.
* `...S1Roi_2.mat` and `...S1Roi_2.jpg` were produced from a **second** ROI of the same cell
  (`...S1Roi_2.roi`, 16 vertices, first vertex 980 / 1659 px) that is not part of the example data.
  The shipped `...S1Roi_1.roi` is a *different* drawing — 30 vertices, at a different location — and
  its coordinates exceed the shipped slice (`y` up to 1659 px against a 933-px image height), so the
  S1 input shown here and the S1 output shown here belong to different drawings.
* The shipped `...S1Roi_2.mat` also carries **five** columns per row of `Res`, where the current
  `TomowidthCaclu.m` writes four; it was saved by an earlier revision of the script.

All of this is fine for reading the workflow, but do not expect to re-derive `...S1Roi_2.mat` from
`...S1Roi_1.roi`.

---

## How to reproduce the example manually

    1. Open demo_data/01_manual_septum_roi/Position_1_Aretomo_bin1-1-23Roi.jpg in Fiji
    2. Draw a polygon around the septum, save it as ...-S0Roi.roi next to the image
    3. In MATLAB, cd into demo_data/01_manual_septum_roi/
    4. Run TomowidthCaclu.m, select ...-S0Roi.roi in the dialog
    5. Click the two endpoints of the septum axis in the figure window
    6. ...-S0Roi.mat and ...-S0Roi.jpg appear next to the ROI;
       compare them with demo_data/02_septum_width_calculation/...-S0Roi.*

---

## Requirements

* MATLAB R2020a or later, with the **Mapping Toolbox** (`polyxpoly`)
* `ReadImageJROI.m` — shipped here and in `SMTanalysis/` (identical copies)
* Fiji (ImageJ) for the ROI drawing
* Input: 2D electron tomography slices, ImageJ polygon ROIs

## Files

| File | Content |
|---|---|
| `TomowidthCaclu.m` | the module script, byte-identical to `Tomo/TomowidthCaclu.m` |
| `ReadImageJROI.m` | the ImageJ-ROI reader, byte-identical to `SMTanalysis/ReadImageJROI.m` |
| `demo_data/01_manual_septum_roi/` | example data 1 — the slices and the drawn ROIs |
| `demo_data/02_septum_width_calculation/` | example data 2 — the published figures and results |
