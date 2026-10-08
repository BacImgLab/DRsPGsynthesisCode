# 3D-dSTORM / 3D-SMLM septum width analysis

`3D-dSTORM/` contains the MATLAB script of the 3D-dSTORM workflow for *Deinococcus radiodurans*.
The workflow quantifies the **relative width of the division septum** from the 2D projection image
that is reconstructed from the 3D-dSTORM data (acquisition thickness ≈ 200 nm). The example used
throughout this document is a WT sample, with one folder of septum profiles per cell-cycle stage
(**S0** and **S1**).

`demo/SMLM/` is a **self-contained package** of that module: the module script, the entry point
`runDemoSMLM.m` with its step wrapper, the example data set in `demo_data/` and this Readme all
live in this single folder, so the folder can be downloaded and run on its own.

The workflow consists of two parts:

| Part | Content | Can be run unattended |
|---|---|---|
| 1 | [ROI selection and profile export in Fiji](#1-roi-selection-and-profile-export-in-fiji) | No — manual line ROIs |
| 2 | [Septum width (FWHM) calculation](#2-septum-width-fwhm-calculation) | **Yes** — covered by `runDemoSMLM.m` |

Part 2 is the only scriptable step of the workflow, and it is the one the automated demo replays
(see [Automated demo](#automated-demo)).

---

## 1. ROI selection and profile export in Fiji

* Reconstruct the 2D projection image of the 3D-dSTORM data.
* Open it in **Fiji/ImageJ** and draw a **line ROI** across the septum of every cell, perpendicular
  to the septum and long enough to cover the whole septum signal.
* Export the intensity profile along the line: this writes one `X,Y` profile per septum —
  column `X` is the distance along the line ROI, column `Y` the grey value.

→ produced: **example data 1** — one profile file per septum ROI, one folder per cell-cycle stage.

## 2. Septum width (FWHM) calculation

* MATLAB script: `dSTORMwidthcaclu.m`
* **Input**: the folder of profile files (example data 1)
* **Output**: `Plots/AllWidths.csv` (one `FileName, HalfPeakWidth` row per profile) plus one
  annotated figure `<profile name>_HalfPeakWidth.png` per profile — **both written inside the
  profile folder**

The script batches over every `*.csv` of the folder it is pointed at and computes the **full width
at half maximum (FWHM)** of each profile:

1. `readmatrix` reads the file; the leading `X,Y` header row (which arrives as `NaN,NaN`) is
   dropped.
2. The half height is `(max(Y) + min(Y)) / 2`.
3. Every crossing of the half height is located by linear interpolation between two consecutive
   samples.
4. The width is the distance between the **first** and the **last** crossing; profiles with fewer
   than two crossings get `NaN`.
5. The annotated figure marks the half height, both crossings and the width; `AllWidths.csv`
   collects all widths.

There are no tunable parameters. The width comes out in the unit of the `X` column, so it is
converted to nanometres by the sampling of the profile: the shipped profiles are sampled every
`0.005` µm, i.e. **X is in micrometre and 1 µm = 1000 nm** (the "5nmGauss" of the file names is the
5 nm sampling/Gaussian of the reconstruction).

---

## Example data

`demo_data/` holds the example data set of the workflow, under the English name of the original
item. 116 files, ≈ 330 KB in total.

| Folder | Content | Read by the demo |
|---|---|---|
| `01_wt_s0_s1_profiles/WTS0/` | 51 septum profiles of the S0 cells: `WT_S0_width_roiN_5nmGuass_500mW_roi1_z450_650_Values.csv`, `X,Y` header plus one row per sample. | step (input) |
| `01_wt_s0_s1_profiles/WTS1/` | 65 septum profiles of the S1 cells, same format. | step (input) |

All profiles are sampled on the same 0.005 µm grid; the line ROIs span 0.59–1.21 µm and the grey
values are background-subtracted (baseline ≈ 0).

---

## Automated demo

    >> cd demo/SMLM
    >> runDemoSMLM

or, from a terminal:

    matlab -batch "cd demo/SMLM; runDemoSMLM"

`runDemoSMLM.m` executes the one unattended step of the workflow, driving the module script that
the folder carries alongside the entry point:

1. `runDemoSMLMWidth.m` — **septum width (FWHM) of every profile**. It runs `dSTORMwidthcaclu.m`
   once per stage. The module script resolves its input from the variable `sourceFolder` (its
   header states that the variable may be defined in the workspace before the script is run, and it
   starts with `clearvars -except sourceFolder`), so the wrapper only sets that variable and the
   published script is not modified. Because the script writes its `Plots/` output *inside* the
   profile folder, the wrapper first assembles an isolated working copy of the stage profiles in
   `demo_output/width/<stage>/work/`, runs the script there and collects the products one level up.

`demo_output/width/` then contains:

    demo/SMLM/demo_output/width/
    ├── SMLM_width_summary.csv            all 116 profiles: FileName, HalfPeakWidth,
    │                                    Stage, HalfPeakWidth_nm
    ├── SMLM_width_stats.csv              per stage: nProfiles, nNaN, nWidths,
    │                                    mean/median/sd/min/max in um
    ├── WTS0/
    │   ├── AllWidths.csv                 the summary as written by the module script
    │   ├── *_HalfPeakWidth.png           one annotated FWHM figure per S0 septum (51)
    │   └── work/                         the isolated working copy used for the run
    └── WTS1/
        ├── AllWidths.csv
        ├── *_HalfPeakWidth.png           one annotated FWHM figure per S1 septum (65)
        └── work/

Expected values (MATLAB R2024a):

* **WTS0**: 51 profiles, none of them without an FWHM. Width **0.1097 ± 0.0206 µm** (median
  0.1068, range 0.0681–0.1775 µm) = **109.7 ± 20.6 nm**.
* **WTS1**: 65 profiles, none of them without an FWHM. Width **0.1062 ± 0.0309 µm** (median
  0.0989, range 0.0564–0.2030 µm) = **106.2 ± 30.9 nm**.
* The console ends with `=== SMLM demo finished in ≈ 2 min ===`.

Expected run time: **≈ 2 min** on a normal desktop computer with MATLAB R2024a (measured 106 s and 129 s), almost all of it in
the 116 figures the module script draws and saves. The demo leaves ≈ 11 MB of output and working copies in
`demo/SMLM/demo_output/`; delete that folder afterwards if disk space matters.

### A note on the units

The `HalfPeakWidth` values of `AllWidths.csv` are in the unit of the profile `X` column, not in
pixels as the workflow description states ("以像素为单位的 FWHM"): a profile sampled from 0 to
0.59–1.21 in steps of 0.005 cannot be a pixel axis. Taking X in micrometre — the 5 nm sampling of
the "5nmGauss" reconstruction — gives septum widths of ≈ 100 nm, the order of magnitude expected
for the septal signal of a 2D projection of a ≈ 200 nm thick 3D-dSTORM reconstruction.

---

## Requirements

* MATLAB R2020a or later; **base MATLAB only** (`readmatrix`, `writecell`, `plot`, `saveas`)
* Fiji (ImageJ) for the ROI drawing and the profile export
* Input: one `X,Y` intensity profile per septum ROI

## Files

| File | Content |
|---|---|
| `dSTORMwidthcaclu.m` | the module script, byte-identical to `3D-dSTORM/dSTORMwidthcaclu.m` |
| `runDemoSMLM.m` | the entry point |
| `runDemoSMLMWidth.m` | the step wrapper that drives the module script per stage |
| `demo_data/01_wt_s0_s1_profiles/` | example data 1 — the septum profiles (WTS0, WTS1) |
