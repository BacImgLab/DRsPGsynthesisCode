# Fluorescence_Lifetime_Imaging

Fluorescence_Lifetime_Imaging contains Matlab scripts and Fiji plugins for processing fluorescence lifetime imaging (FLIM) data of *Deinococcus radiodurans*, with the aim of quantifying spatial lifetime distributions and removing background interference.

---

## 1. Data Preparation

- Export intensity and lifetime maps from **LAS X FLIM/FCS** (Leica software).  
- Set export thresholds:
  - **Intensity**: 0–65535  
  - **Lifetime**: 0–10 (e.g., nanoseconds)  
- Save both types of images for all fields of view.

---

## 2. Image Stacking

- Run **Section 1** of the Matlab script `lifeTcalculationbernsen.m`.  
- This step loads all exported intensity and lifetime images and stacks them separately into intensity and lifetime image sequences.

---

## 3. Background Masking

- Use Fiji to remove cytoplasmic and peripheral background signals by applying thresholding.  
- Fiji Plugin: `bersenThtest`  
- This sets background pixels to zero in both the intensity and lifetime stacks.

---

## 4. Lifetime & Intensity Quantification

- Run **Section 3** of the script `lifeTcalculationbernsen.m`.  
- This step reads the background-masked stacks and computes the intensity and lifetime values pixel by pixel.

---

## Requirements

- Matlab 2020a or later  
- Fiji (ImageJ) with:
  - `bersenThtest` plugin  
- LAS X FLIM/FCS (for data acquisition and export)

---

## How to Use

1. Export lifetime and intensity images from LAS X FLIM/FCS (0–65535 and 0–10 scale respectively)  
2. Run Section 1 of `lifeTcalculationbernsen.m` to create stacks of all images  
3. Use Fiji and the `bersenThtest` plugin to mask background regions (set to 0)  
4. Run Section 3 of `lifeTcalculationbernsen.m` to calculate per-pixel lifetime and intensity values


---

## Automated demo

Steps 2 and 4 - the two steps that need no user interaction and no Fiji - are covered by the repository demo:

    >> cd demo/Lifetime
    >> runDemoLifetime

or, from a terminal: `matlab -batch "cd demo/Lifetime; runDemoLifetime"`. It replays Sections 1
and 3 of `lifeTcalculationbernsen.m` on the published example data set and writes everything into
`demo/Lifetime/demo_output/`. See `demo/Lifetime/Readme.md` for the full description of the
workflow, of the example data and of how the demo drives the module script.
