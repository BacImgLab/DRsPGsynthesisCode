function runDemoColoc()
%RUNDEMOCOLOC Run the two unattended steps of the multicolor colocalization module.
%
%   The multicolor work flow (Multicolor Colocalization) consists of two parts.
%
%   The upstream part cannot be replayed unattended, because it needs Fiji/ImageJ
%   or a human: the drift/Z-projection step and the chromatic-aberration step
%   (MultiChannelDriftCorrection.m, MultiChannelChromaticAberration.m), the
%   PureDenoise plugin, the Cellpose3 + DeCNN cell-cycle classification, and the
%   manual drawing of the S0/S1 septum lines in DemoDR_stage2.m/DemoDR_stage35.m.
%
%   Only the two downstream analysis steps run without any user interaction, and
%   those are exactly what this demo executes on the published example data set:
%
%     Step 1 : Colocalization/DemoS1AnalysisW.m
%              S1 septum demograph analysis - fits two Gaussians to every column
%              of the smoothed S1 demograph and measures the distance between the
%              protein ridge (channel 4) and the membrane ridge (channel 2).
%              Input  : example data 6 (smoothed demograph), channels 2 and 4 of S1
%              Output : WmemData.mat and Wmem.tif (named AmemData/Amem.tif in the
%                       manuscript text; see README, "A note on the file names")
%
%     Step 2 : Colocalization/PCCcaclu.m
%              Pairwise Pearson correlation coefficient of the channels of every
%              single cell, plus the merged table.
%              Input  : example data 4 (per-cell septum profiling)
%              Output : <stage>_PCC_Values.csv for every cell-cycle stage and the
%                       combined Merged_PCC_Values.csv
%
%   This is a self-contained package: the two wrapper scripts, the 18 module
%   scripts they drive, the published example data set (demo_data/) and the
%   module documentation (Readme.md) all live in this folder, so the folder can
%   be downloaded and run on its own.
%
%   Usage (MATLAB):
%       >> cd demo/Colocalization
%       >> runDemoColoc
%
%   From a terminal:
%       matlab -batch "cd demo/Colocalization; runDemoColoc"
%
%   Requirements: step 2 runs on base MATLAB. Step 1 additionally needs the
%   Optimization Toolbox (lsqcurvefit) and the Image Processing Toolbox
%   (rgb2ind). Neither Fiji/ImageJ nor the experimental raw data are needed.
%
%   All output goes to demo/Colocalization/demo_output/, which is not tracked by
%   git; the package and the demo data are never modified.
%
%   Expected output and run time are listed in Readme.md (section "Automated
%   demo") and in the Demo section of the repository README.md.
%   The simulated-data demo (../runDemo.m) is independent of this one.

t0 = tic;
status = {};

demoDir  = fileparts(mfilename('fullpath'));
dataDir  = fullfile(demoDir, 'demo_data');
outRoot  = fullfile(demoDir, 'demo_output');
if exist(outRoot, 'dir')
    rmdir(outRoot, 's');
end
mkdir(outRoot);

%% --- Step 1: S1 septum demograph analysis ---
fprintf('\n=== Step 1/2: S1 septum demograph analysis (runDemoColocS1.m) ===\n');
try
    runDemoColocS1(fullfile(dataDir, '06_smoothed_demograph', 'S1'), ...
                   fullfile(outRoot, 'S1_peak_distance'));
    status{end+1} = 'Step 1 (S1 septum demograph analysis) : OK'; %#ok<*AGROW>
catch ME
    status{end+1} = ['Step 1 (S1 septum demograph analysis) : FAILED - ', ME.message];
    rethrow(ME);
end

%% --- Step 2: PCC between the channels ---
fprintf('\n=== Step 2/2: pairwise PCC between the channels (runDemoColocPcc.m) ===\n');
try
    runDemoColocPcc(fullfile(dataDir, '04_septum_profiles'), ...
                    fullfile(outRoot, 'PCC'));
    status{end+1} = 'Step 2 (pairwise PCC)                  : OK';
catch ME
    status{end+1} = ['Step 2 (pairwise PCC)                  : FAILED - ', ME.message];
    rethrow(ME);
end

%% --- Summary ---
fprintf('\n=== Colocalization demo finished in %.1f s ===\n', toc(t0));
for k = 1:numel(status)
    fprintf('  %s\n', status{k});
end
fprintf('Input data : %s\n', dataDir);
fprintf('Output     : %s\n', outRoot);

end
