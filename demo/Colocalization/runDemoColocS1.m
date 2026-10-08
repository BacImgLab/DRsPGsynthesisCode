function runDemoColocS1(dataRoot, outDir)
%RUNDEMOCOLOCS1 Measure the distance between the two signal ridges in S1 septa.
%
%   Step 1 of the multicolor colocalization demo. It drives
%   DemoS1AnalysisW.m, which fits a sum of two Gaussians to every
%   column of the smoothed S1 demograph and records the first peak position for
%   the protein channel (channel 4) and for the membrane channel (channel 2).
%   The distance between the two ridges across the septum is the quantity
%   reported for S1 in the manuscript.
%
%   runDemoColocS1()            uses demo_data/06_smoothed_demograph/S1
%                               and writes into demo_output/S1_peak_distance
%   runDemoColocS1(dataRoot, outDir)
%
%   The module script resolves its two inputs by name in the current folder
%   (demo_S1_smooth_W.mat for the protein, demo_S1_smooth_mem.mat for the
%   membrane), so a working copy is assembled in <outDir>/work.
%
%   Requirements: Optimization Toolbox (lsqcurvefit) and Image Processing
%   Toolbox (rgb2ind). The other demo step runs on base MATLAB only.
%
%   Output:
%       WmemData.mat   mui_poi, mu1_mem (80 x 1, rows = bin along the cell cycle)
%                      and ParameterAll with the fit parameters of every column
%       Wmem.tif       80-page figure stack, one page per demograph column
%                      (uncompressed, roughly 55 MB)
%   The manuscript text refers to these two files as AmemData and Amem.tif.

if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile(fileparts(mfilename('fullpath')), 'demo_data', '06_smoothed_demograph', 'S1');
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'demo_output', 'S1_peak_distance');
end
scriptDir = fileparts(mfilename('fullpath'));

if exist('lsqcurvefit', 'file') ~= 2 || exist('rgb2ind', 'file') ~= 2
    error('runDemoColocS1:missingToolbox', ...
          ['DemoS1AnalysisW.m needs the Optimization Toolbox (lsqcurvefit) and the ' ...
           'Image Processing Toolbox (rgb2ind); both were not found.']);
end

if exist(outDir, 'dir')
    rmdir(outDir, 's');
end
workDir = fullfile(outDir, 'work');
mkdir(workDir);

copyfile(fullfile(dataRoot, 'demo_sm_norm4new.mat'), fullfile(workDir, 'demo_S1_smooth_W.mat'));
copyfile(fullfile(dataRoot, 'demo_sm_norm2new.mat'), fullfile(workDir, 'demo_S1_smooth_mem.mat'));
copyfile(fullfile(scriptDir, 'DemoS1AnalysisW.m'), workDir);
copyfile(fullfile(scriptDir, 'demofitGauss2.m'), workDir);

t0 = tic;
fitStep(workDir);
elapsed = toc(t0);

copyfile(fullfile(workDir, 'WmemData.mat'), outDir);
copyfile(fullfile(workDir, 'Wmem.tif'), outDir);

R = load(fullfile(outDir, 'WmemData.mat'));
dist = abs(R.mui_poi - R.mu1_mem);
info = imfinfo(fullfile(outDir, 'Wmem.tif'));
tifBytes = dir(fullfile(outDir, 'Wmem.tif')).bytes;

fprintf('  columns fitted          : %d (two Gaussians per column)\n', numel(R.mui_poi));
fprintf('  protein peak  (ch4)     : mean %.2f  range %.1f .. %.1f px\n', ...
        mean(R.mui_poi), min(R.mui_poi), max(R.mui_poi));
fprintf('  membrane peak (ch2)     : mean %.2f  range %.1f .. %.1f px\n', ...
        mean(R.mu1_mem), min(R.mu1_mem), max(R.mu1_mem));
fprintf('  ridge separation        : mean %.2f  range %.1f .. %.1f px  (NaN %d)\n', ...
        mean(dist, 'omitnan'), min(dist), max(dist), sum(isnan(dist)));
fprintf('  Wmem.tif                : %d pages, %dx%d, %.1f MB\n', ...
        numel(info), info(1).Width, info(1).Height, tifBytes/1e6);
fprintf('  fitting time            : %.1f s\n', elapsed);
fprintf('  output                  : %s\n', outDir);

end

function fitStep(workDir)
%FITSTEP Run DemoS1AnalysisW.m inside its own workspace.
%   The script starts with "clear; clc; close all;" and closes all figures at
%   the end, so it is kept two levels away from the caller.
here = pwd;
cd(workDir);
fitInner();
cd(here);
end

function fitInner()
DemoS1AnalysisW; %#ok<NASGU>
end
