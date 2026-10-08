# Single-molecule tracking (SMT)

`SMTanalysis/` contains the MATLAB scripts and the Fiji/ImageJ dependencies of the
single-molecule tracking (SMT) workflow for *Deinococcus radiodurans*. The workflow links
single-molecule localisations into trajectories, segments them into unidirectional pieces,
classifies their motion state and quantifies the speed distribution of the directed population.
The example used throughout this document is an **FtsW** tracking experiment (condition label
`FtsW-S1-Cef`, the label `dataprocessSMTWCF.m` writes into its output file name).

`demo/SMT/` is a **self-contained package** of that module: the 31 module scripts, the entry point
`runDemoSMT.m` with its step wrapper `runDemoSMTVelocity.m`, the example data item the demo reads
in `demo_data/` and this Readme all live in this single folder, so the folder can be downloaded
and run on its own.

The workflow is documented against the four items of the published example data set, referred to
below as **example data 1** … **example data 4**. This package ships the one item that the demo
reads — example data 4, the classified trajectory data — in `demo_data/`, under the English name
of the original item; the data file itself is unchanged. The other three items are not needed by
the demo and are available from the corresponding author.

The module consists of three parts:

| Part | Content | Can be run unattended |
|---|---|---|
| 1 | [Trajectory data preprocessing](#1-trajectory-data-preprocessing) | No — Fiji/ThunderSTORM, Cellpose3 and GUI file dialogs |
| 2 | [Trajectory segmentation](#2-trajectory-segmentation) | No — interactive App Designer app |
| 3 | [Trajectory classification and speed fitting](#3-trajectory-classification-and-speed-fitting) | 3.1 No — save dialog and button press; 3.2 **Yes** — covered by `runDemoSMT.m` |

Step 3.2 is the only step of the workflow that consumes a plain, scriptable input, and it is the
one that the automated demo replays (see [Automated demo](#automated-demo)).

---

## 1. Trajectory data preprocessing

### 1.1 ROI selection and cropping

A region of interest (ROI) with uniform illumination across the splitter is identified on
fluorescent-bead images.

* MATLAB script: `SMTdataPrepare.m`

The script organises the raw dual-view (splitter) tracking data — bright-field (`B`),
fluorescence (`F`) and SMT/tracking (`S`) 16-bit TIFF images — crops the predefined ROI from every
image in a batch and writes the organised output into a new folder for downstream analysis. The
cropped images are 743 × 1022 pixel.

### 1.2 Single-molecule localisation

Individual molecules are localised from the raw movies with ThunderSTORM in Fiji.

* Fiji macro: `MacroThurderSTORMDrJ.ijm`

The macro automates the localisation and writes the molecular positions as a ThunderSTORM CSV
table (frame, x and y in nm, localisation precision, photon counts, …).

→ produced: **example data 1**.

### 1.3 Chromatic aberration correction

The chromatic offset between the fluorescence channel and the tracking channel is corrected.

* MATLAB script: `CACorrectionofSMTsplitter.m`

The script matches the localisation points of the ThunderSTORM output between the two channels,
fits a polynomial transformation (moving → fixed channel), applies it to the raw images and saves
the aligned images.

### 1.4 Cell segmentation

The cell contours are segmented from the bright-field images with **Cellpose3** and saved as
binary masks; `roitextRead.m` reads the Cellpose output back into MATLAB, smooths the contours and
fits a rotated bounding box per cell.

### 1.5 Trajectory linking

The single-molecule positions of adjacent frames are linked into continuous trajectories.

* MATLAB GUI: `spotsLinking.m` (requires user interaction)

The GUI takes the localisation table of step 1.2, links the spots across frames and writes one
`.mat` file per field of view with the linked trajectories.

→ produced: **example data 2**.

The trajectories are then assigned to the segmented cells and refined per cell with
`traceAnalysisMain.m` / `traceInROI.m` (which computes the per-trajectory MSD and crops the
bright-field and fluorescence images around each cell), and the per-trajectory summary figures are
produced with `tracePlot.m`.

---

## 2. Trajectory segmentation

Every trajectory is divided into segments of a single direction of motion
(unidirectional segments).

* MATLAB App: `RefineTraceSegDr.mlapp` — interactive
* Dependencies: `MSDsingle2D.m`, `linfitR.m`

The app displays every trajectory together with its bright-field and fluorescence ROI images,
lets the user place the segmentation points and refits every segment linearly. The result is
stored as the `Seg` field of the refined-trajectory structure.

→ produced: **example data 3**.

---

## 3. Trajectory classification and speed fitting

### 3.1 Segment motion-state classification

Every segment is classified as **directional**, **diffusive** or **stationary**.

* MATLAB script: `statesClassifyDr.m`
* Dependencies shipped with the module: `rcdfCal.m`, `mergeSegData.m`, `tracedropoutXY.m`,
  `addVoneFLtrajs.m`, `getProbR.m`, `histLog_xy.m`, `kusumi_xy.m`, `kelsey2.m`,
  `colorCodeTracePlot.m`
* **Requires user interaction**: the script asks for an output file name (`uiputfile`) and waits
  for a button press before it continues.

The script first simulates a library of confined random walks with `rcdfCal.m` (section 0), merges
the segmented trajectories with `mergeSegData.m` (section 1), computes for every segment the
noise-to-displacement ratio `R` by bootstrapping with random dropout (`tracedropoutXY.m`) and the
probability `P` of directed motion against the simulated library (`addVoneFLtrajs.m`,
`getProbR.m`) (section 2), and finally classifies the segments with the thresholds `Rmax1`,
`Rmax2`, `StDmax` and `Pmin` (section 4). The per-segment table is written as the 7-column matrix
`SegVSRPTC` — velocity, standard deviation, `R`, `P`, dwell time, location flag, state — together
with the speed and dwell-time vectors `Vd`/`Vs` and `Td`/`Ts` of the directional and stationary
populations.

→ produced: **example data 4**.

### 3.2 Speed-distribution fitting

The speed distribution of the directed segments is fitted with a single- and a double-log-normal
model, and the mean speed and the population proportions are computed.

* MATLAB script: `dataprocessSMTWCF.m`
* Dependencies: `CDF_logCalc.m`, `logn1cdf.m`, `logn2cdf.m`
* **Input**: example data 4 — the classified trajectory data (`FtsW-all.mat`)
* **Output**: the fitted figures and the fitting results (`FtsW-S1-Cef.mat`)
* This step runs unattended and is covered by the demo.

The script selects the directional (`state = 1`) and stationary (`state = 3`) segments of location
flag 3 from the classification table, computes the empirical CDF of the directional speeds above
the 1 nm/s threshold on 21 log-spaced bins between 1 and 100 nm/s, fits the single- and the
double-log-normal CDF models by least squares, bootstraps both fits over 200 resamples and
reconstructs the two probability density functions for the comparison with the speed histogram.
Everything is stored in the `Result` struct of the output `.mat` file.

Requires the Optimization Toolbox (`lsqcurvefit`) and the Statistics and Machine Learning Toolbox
(`logncdf`, `lognpdf`, `bootstrp`).

---

## Example data

The published example data set of the SMT workflow consists of four items. This package ships
**example data 4** (1 file, 54 KB); the other three items are described here because the workflow
text refers to them, and are available from the corresponding author.

| Supplement item | Folder | Content | Shipped here |
|---|---|---|---|
| 1. ThunderSTORM localisation data | `01_localizations/` | `track0.csv` — 14 673 single-molecule localisations of 200 frames in the ThunderSTORM column format (frame, x/y in nm, sigma, intensity, background, χ², x/y uncertainty). | no |
| 2. Spotslink trajectory data | `02_linked_tracks/` | `long-10-Coord-track0.mat` — `Spots` (200 frames, 14 673 localisations) and `tracksFinal` (419 linked trajectories, 8 845 localisations, 10–185 points each). | no |
| 3. Segmented trajectory data | `03_segmented_traces/` | `TraceRefine-test.mat` — `TrackRefine`, a 372 × 1 cell array with one struct per trajectory: the raw track, the 14-column coordinate matrix, the per-trajectory MSD, the rotated cell ROI and the cropped bright-field and fluorescence images. Pixel size 110 nm, frame time 1 s. | no |
| 4. Classified trajectory data | `04_classified_traces/` | `FtsW-all.mat` — the output of step 3.1: `SegVSRPTC` (1405 × 7 = velocity, StD, R, P, dwell time, location flag, state), the classification thresholds (`Pmin = 0.8`, `Rmax1 = 0.5`, `Rmax2 = 0.3`) and the speed and dwell-time vectors `Vd`/`Vs`/`Td`/`Ts` (487 directional and 917 stationary segments). | **yes** — demo input |

Two notes on this data set:

* `03_segmented_traces/TraceRefine-test.mat` is **31 MB**, above the 25 MiB single-file limit of
  the GitHub web uploader — one more reason to ship only the item the demo needs.
* The `TrackRefine` structures of item 3 do **not** contain the `Seg` field that `mergeSegData.m`
  and `statesClassifyDr.m` expect, so step 3.1 cannot be replayed from the shipped files. The demo
  therefore starts at step 3.2, whose input (`SegVSRPTC`) is example data 4.

---

## Requirements

* MATLAB R2020a or later
* Fiji with ThunderSTORM for the single-molecule localisation
* Cellpose3 for the cell segmentation
* For `dataprocessSMTWCF.m` (step 3.2): the Optimization Toolbox (`lsqcurvefit`) and the Statistics
  and Machine Learning Toolbox (`logncdf`, `lognpdf`, `bootstrp`). The other steps run on base
  MATLAB.
* Windows 10 / 11 (64-bit)

Input data are dual-view single-molecule movies (`.tif`), bright-field images and the ThunderSTORM
localisation tables.

---

## How to use

1. Crop the ROIs with `SMTdataPrepare.m`.
2. Run `MacroThurderSTORMDrJ.ijm` in Fiji to localise the single molecules.
3. Correct the channel offset with `CACorrectionofSMTsplitter.m`.
4. Segment the cells with Cellpose3 and read the masks back with `roitextRead.m`.
5. Link the trajectories with the `spotsLinking.m` GUI, then assign them to the cells with
   `traceAnalysisMain.m` / `traceInROI.m`.
6. Segment the trajectories into unidirectional segments with the `RefineTraceSegDr` app.
7. Classify the motion states with `statesClassifyDr.m`.
8. Fit the speed distributions with `dataprocessSMTWCF.m`.
9. Visualise the segmented trajectories with `tracePlotafterseg.m`, and remove the drifting cells
   with `ReadImageJROI.m` and `RemoveDrift.m`.

## Automated demo

Step 3.2 — the only step with no user interaction — is covered by the demo that ships with this
folder:

    >> cd demo/SMT
    >> runDemoSMT

or, from a terminal: `matlab -batch "cd demo/SMT; runDemoSMT"`. It runs `runDemoSMTVelocity.m`
(part 3.2) on the example data and writes everything into `demo/SMT/demo_output/`. The expected
console output, the expected values and the run time are listed in the [Demo section of the
repository README](../../README.md#demo).

Neither the package nor the example data are modified by the demo: the wrapper assembles an
isolated working copy below `demo_output/`, which is not tracked by git.

## Note on the file and variable names

A few names differ between the manuscript text and the files in this folder:

| In the manuscript text | In this folder | Where |
|---|---|---|
| `MacroThurderSTORMDr.ijm` | `MacroThurderSTORMDrJ.ijm` | Fiji macro of step 1.2 |
| `Spotslink.m` | `spotsLinking.m` | trajectory linking, step 1.5 |
| `mergeSegDrJ.m` | `mergeSegData.m` | dependency of step 3.1 |
| `segUpdate.m`, `segPRplot.m` | local functions inside `statesClassifyDr.m` | dependencies of step 3.1 |
| `MSDcaclulate_2d.m` | `MSDsingle2D.m` | dependency of steps 2 and 3.1 |

Two further notes:

* The module script `dataprocessSMTWCF.m` loads `FtsW-all.mat` and expects the per-segment table
  under the variable name **`DataSMT`**, whereas the published example data set — and the save
  statement of `statesClassifyDr.m` — store exactly the same matrix under the name **`SegVSRPTC`**.
  The published file is left untouched: the demo wrapper adds the alias `DataSMT = SegVSRPTC` in
  its isolated working copy and reports this in the console.
* The GUIDE figure `spotsLinking.fig` that `spotsLinking.m` opens is not part of the repository.

## Contents of this folder

| File | Role |
|---|---|
| `Readme.md` | this document |
| `SMTdataPrepare.m` | 1.1 — batch ROI cropping of the dual-view raw data |
| `MacroThurderSTORMDrJ.ijm` | 1.2 — ThunderSTORM single-molecule localisation in Fiji |
| `CACorrectionofSMTsplitter.m` | 1.3 — chromatic aberration correction of the splitter |
| `roitextRead.m` | 1.4 — reads the Cellpose masks, smooths the contours, rotated bounding boxes |
| `spotsLinking.m` | 1.5 — links the localisations into trajectories (GUI) |
| `traceAnalysisMain.m` | 1.5 — main script: trajectories → cell ROIs → refined traces |
| `traceInROI.m` | 1.5 — assigns trajectories to cells, per-trajectory MSD and ROI crops |
| `tracePlot.m` | 1.5 — per-trajectory summary figures |
| `RefineTraceSegDr.mlapp` | 2 — interactive segmentation into unidirectional segments |
| `linfitR.m` | linear fit of a trace, used by 2 and 3.1 |
| `MSDsingle2D.m` | mean squared displacement of a single trajectory |
| `statesClassifyDr.m` | 3.1 — classification into directional / diffusive / stationary |
| `rcdfCal.m` | 3.1 — simulation of the confined-trajectory library |
| `mergeSegData.m` | 3.1 — merges the segmented trajectories into `mergedAllSeg.mat` |
| `tracedropoutXY.m` | 3.1 — bootstrapping with random dropout |
| `addVoneFLtrajs.m` | 3.1 — R distribution of the velocity-shifted simulations |
| `getProbR.m` | 3.1 — probability of an R interval |
| `histLog_xy.m` | 3.1 — log-binned speed histogram |
| `kusumi_xy.m` | 3.1 (optional) — Kusumi confinement fit of the stationary MSD |
| `kelsey2.m` | 3.1 (optional) — least-squares helper of the Kusumi fit |
| `colorCodeTracePlot.m` | trajectory plot coloured by time |
| `dataprocessSMTWCF.m` | 3.2 — single/double log-normal fit of the speed distribution |
| `CDF_logCalc.m` | 3.2 — empirical CDF on a given bin grid |
| `logn1cdf.m` | 3.2 — single log-normal CDF model |
| `logn2cdf.m` | 3.2 — double log-normal CDF model |
| `tracePlotafterseg.m` | 4 — multi-page TIFF of the segmented trajectories |
| `unwrapTraj.m` | unwraps a track onto the cylindrical cell surface |
| `DrunwrapX.m`, `DrunwrapY.m` | interactive unwrapping GUIs |
| `ReadImageJROI.m` | 5 — reads ImageJ `.roi` files |
| `RemoveDrift.m` | 5 — removes the trajectories of drifting cells |
| `runDemoSMT.m` | demo entry point — runs the unattended step and prints a summary |
| `runDemoSMTVelocity.m` | demo step 1 — assembles an isolated working copy and drives `dataprocessSMTWCF.m` on example data 4 |
| `demo_data/` | the published example data set, item 4 (`04_classified_traces/FtsW-all.mat`, 54 KB) |
