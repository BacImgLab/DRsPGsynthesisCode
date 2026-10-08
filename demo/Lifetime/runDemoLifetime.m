function runDemoLifetime()
%RUNDEMOLIFETIME Run the two unattended steps of the FLIM (lifetime) demo.
%
%   Usage (MATLAB):
%       >> cd demo/Lifetime
%       >> runDemoLifetime
%
%   From a terminal:
%       matlab -batch "cd demo/Lifetime; runDemoLifetime"
%
%   The demo replays the two steps of the fluorescence lifetime imaging workflow
%   that need no human interaction, on the published example data set:
%
%       Step 1  runDemoLifetimeStacking.m   Section 1 of
%                                           lifeTcalculationbernsen.m
%                                           -> intensity and lifetime stacks
%       Step 2  runDemoLifetimeQuantify.m   Section 3 of
%                                           lifeTcalculationbernsen.m
%                                           -> per-pixel lifetime and intensity
%
%   The steps that are not replayed are described in Readme.md: the LAS X
%   FLIM/FCS export happens in the acquisition software, and the background mask
%   is produced by the Fiji plugin `bersenThtest`, which cannot run from MATLAB.
%   The demo therefore starts from the exported images (example data 1) for
%   step 1 and from the published binary mask (example data 3) for step 2.
%
%   Requirements: base MATLAB only. Step 1 uses imread/imwrite, step 2 adds
%   imquantize, histcounts and bar. Neither Fiji/ImageJ nor the raw acquisition
%   data are needed.
%
%   All output goes to demo/Lifetime/demo_output/, which is not tracked by git;
%   the repository tree and the example data are never modified (both wrappers
%   assemble an isolated working copy below demo_output/).
%
%   Expected output and run time are listed in Readme.md.

t0 = tic;

demoDir = fileparts(mfilename('fullpath'));
fprintf('=== FLIM (lifetime) demo ===\n');
fprintf('  module   : %s\n', fullfile(demoDir, 'lifeTcalculationbernsen.m'));
fprintf('  data     : %s\n', fullfile(demoDir, 'demo_data'));
fprintf('  output   : %s\n\n', fullfile(demoDir, 'demo_output'));

fprintf('Step 1 - image stacking\n');
s1 = runDemoLifetimeStacking();

fprintf('\nStep 2 - lifetime and intensity quantification\n');
s2 = runDemoLifetimeQuantify(s1.outDir);

fprintf('\n=== FLIM demo finished in %.1f s ===\n', toc(t0));
fprintf('  step 1 : %d images stacked into 2 stacks\n', ...
        s1.fieldsOfView * s1.imagesPerField * 2);
fprintf('           identical to example data 2 : intensity %d of %d pages, lifetime %d of %d pages\n', ...
        s1.pagesIdentical(1), numel(imfinfo(s1.intensityStack)), ...
        s1.pagesIdentical(2), numel(imfinfo(s1.lifetimeStack)));
fprintf('  step 2 : %d signal pixels, lifetime %.3f +- %.3f ns, intensity %.1f +- %.1f\n', ...
        s2.signalPixels, s2.lifetimeMean, s2.lifetimeSd, s2.intensityMean, s2.intensitySd);
end
