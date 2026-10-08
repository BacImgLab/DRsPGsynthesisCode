function runDemoSMLM()
%RUNDEMOSMLM Run the unattended step of the 3D-dSTORM / 3D-SMLM septum-width module.
%
%   Usage (MATLAB):
%       >> cd demo/SMLM
%       >> runDemoSMLM
%
%   From a terminal:
%       matlab -batch "cd demo/SMLM; runDemoSMLM"
%
%   The 3D-dSTORM septum-width work flow consists of two parts, and only the
%   second one runs without a human:
%
%       Step 1  ROI selection (Fiji/ImageJ) - NOT replayed
%               A line ROI is drawn by hand across the septum of every cell on
%               the 2D projection of the 3D-dSTORM data, perpendicular to the
%               septum and wide enough to cover the whole septum signal. Fiji
%               then exports the intensity profile along that line, which is the
%               example data set shipped in demo_data/01_wt_s0_s1_profiles/.
%
%       Step 2  dSTORMwidthcaclu.m - THIS demo
%               Full width at half maximum (FWHM) of every profile: half height
%               between the maximum and the minimum of the curve, linear
%               interpolation of the two crossings, width = distance between
%               them. One annotated figure per profile plus a summary table.
%
%   This is a self-contained package: the module script, the wrapper, the
%   published example data set (demo_data/) and the module documentation
%   (Readme.md) all live in this folder, so the folder can be downloaded and run
%   on its own.
%
%   Requirements: base MATLAB only (readmatrix, writecell, plot). Neither
%   Fiji/ImageJ nor the raw dSTORM localisation lists are needed.
%
%   All output goes to demo/SMLM/demo_output/, which is not tracked by git; the
%   package and the example data are never modified (the wrapper assembles an
%   isolated working copy below demo_output/).
%
%   Expected output and run time are listed in Readme.md.

t0 = tic;

demoDir = fileparts(mfilename('fullpath'));
fprintf('=== 3D-dSTORM / 3D-SMLM septum-width demo ===\n');
fprintf('  module   : %s\n', fullfile(demoDir, 'dSTORMwidthcaclu.m'));
fprintf('  data     : %s\n', fullfile(demoDir, 'demo_data', '01_wt_s0_s1_profiles'));
fprintf('  output   : %s\n\n', fullfile(demoDir, 'demo_output', 'width'));

fprintf('Step - septum width (FWHM) of every profile\n');
S = runDemoSMLMWidth();

fprintf('\n=== SMLM demo finished in %.1f s ===\n', toc(t0));
for k = 1:height(S.stats)
    fprintf('  %-5s : %2d profiles, FWHM %.1f +- %.1f nm (median %.1f nm, %.1f - %.1f nm)\n', ...
            S.stats.Stage(k), S.stats.nProfiles(k), ...
            1000*S.stats.mean_um(k), 1000*S.stats.sd_um(k), ...
            1000*S.stats.median_um(k), 1000*S.stats.min_um(k), 1000*S.stats.max_um(k));
end
fprintf('  output   : %s\n', S.outDir);

end
