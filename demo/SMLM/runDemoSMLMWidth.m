function S = runDemoSMLMWidth(dataRoot, outDir)
%RUNDEMOSMLMWIDTH Septum width (FWHM) of every 3D-dSTORM septum profile.
%
%   The only step of the 3D-dSTORM / 3D-SMLM septum-width module that runs
%   unattended. It drives dSTORMwidthcaclu.m on the published example data set,
%   which is a folder of one fluorescence-intensity profile per septum ROI.
%
%   runDemoSMLMWidth()             uses demo_data/01_wt_s0_s1_profiles and
%                                  writes into demo_output/width
%   runDemoSMLMWidth(dataRoot, outDir)
%
%   Expected input layout below dataRoot (one folder per cell-cycle stage):
%       WTS0/*.csv     profile of every S0 septum: comma separated, "X,Y",
%                      one header row, X = distance along the line ROI in
%                      micrometre (sampled every 0.005 um = 5 nm), Y = grey value
%       WTS1/*.csv     the same for the S1 septa
%
%   dSTORMwidthcaclu.m processes every "*.csv" of the folder it is pointed at
%   and writes a "Plots" subfolder *inside that folder* - one annotated figure
%   per profile plus the summary AllWidths.csv. It resolves its input from the
%   variable `sourceFolder` (it starts with "clearvars -except sourceFolder"), so
%   the driver only has to set that variable; the script itself is not modified.
%   An isolated working copy of the profiles is therefore assembled in
%   <outDir>/<stage>/work, so that neither the repository nor the example data
%   are written to.
%
%   Output written into outDir:
%       WTS0/AllWidths.csv            the summary as written by the module script
%       WTS0/*_HalfPeakWidth.png      one annotated FWHM figure per S0 septum
%       WTS0/work/                    the isolated working copy used for the run
%       WTS1/...                      the same for the S1 septa
%       SMLM_width_summary.csv        Stage, FileName, HalfPeakWidth, HalfPeakWidth_nm
%       SMLM_width_stats.csv          per stage: n, nNaN, mean, median, sd, min, max
%
%   S returns the per-stage statistics and the merged table.
%
%   Runs on base MATLAB - no toolbox is required (readmatrix, writecell, plot).

if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile(fileparts(mfilename('fullpath')), ...
                        'demo_data', '01_wt_s0_s1_profiles');
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'demo_output', 'width');
end
scriptDir = fileparts(mfilename('fullpath'));

if exist(outDir, 'dir')
    rmdir(outDir, 's');
end
mkdir(outDir);

stages = {'WTS0', 'WTS1'};
merged = table();
stats  = table();

for k = 1:numel(stages)
    stage = stages{k};
    src   = fullfile(dataRoot, stage);
    if ~exist(src, 'dir')
        error('runDemoSMLMWidth:missingInput', 'Input folder not found: %s', src);
    end

    % --- assemble an isolated working copy of the profiles ---
    stageDir = fullfile(outDir, stage);
    workDir  = fullfile(stageDir, 'work');
    mkdir(workDir);
    copyfile(fullfile(src, '*.csv'), workDir);
    copyfile(fullfile(scriptDir, 'dSTORMwidthcaclu.m'), workDir);

    % --- run the module script in its own workspace ---
    oldDir = pwd;
    widthStep(workDir);
    cd(oldDir);

    % --- collect the products out of the working copy ---
    plotsDir  = fullfile(workDir, 'Plots');
    allWidths = fullfile(plotsDir, 'AllWidths.csv');
    if ~isfile(allWidths)
        error('runDemoSMLMWidth:noOutput', ...
              'dSTORMwidthcaclu.m did not produce AllWidths.csv in %s', plotsDir);
    end
    copyfile(allWidths, stageDir);
    copyfile(fullfile(plotsDir, '*.png'), stageDir);

    T = readtable(allWidths, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
    T.Stage = repmat(string(stage), height(T), 1);
    T.HalfPeakWidth_nm = T.HalfPeakWidth * 1000;      % X is in micrometre
    merged = [merged; T]; %#ok<AGROW>

    w = T.HalfPeakWidth(~isnan(T.HalfPeakWidth));
    stats = [stats; table(string(stage), height(T), sum(isnan(T.HalfPeakWidth)), ...
                          numel(w), mean(w), median(w), std(w), min(w), max(w), ...
             'VariableNames', {'Stage', 'nProfiles', 'nNaN', 'nWidths', ...
                               'mean_um', 'median_um', 'sd_um', 'min_um', 'max_um'})]; %#ok<AGROW>

    fprintf('  %-5s : %2d profiles, %2d without an FWHM (no half-height crossing)\n', ...
            stage, height(T), sum(isnan(T.HalfPeakWidth)));
    fprintf('          FWHM %.4f +- %.4f um  (median %.4f, range %.4f - %.4f) = %.1f +- %.1f nm\n', ...
            mean(w), std(w), median(w), min(w), max(w), 1000*mean(w), 1000*std(w));
    fprintf('          output : %s\n', stageDir);
end

writetable(merged, fullfile(outDir, 'SMLM_width_summary.csv'));
writetable(stats,  fullfile(outDir, 'SMLM_width_stats.csv'));

S.merged  = merged;
S.stats   = stats;
S.outDir  = outDir;
S.stages  = stages;

end

function widthStep(workDir)
%WIDTHSTEP Run dSTORMwidthcaclu.m inside its own workspace.
%   The script starts with "clearvars -except sourceFolder; clc;", which wipes
%   the variables of the function it is called from - so it is called from a
%   dedicated function whose only variable is `sourceFolder`.
here = pwd;
cd(workDir);
widthInner();
cd(here);
end

function widthInner()
%dSTORMwidthcaclu.m reads `sourceFolder` and falls back to a hard-coded path
%when it is absent; `clearvars -except sourceFolder` makes sure it survives.
sourceFolder = pwd; %#ok<NASGU>
dSTORMwidthcaclu; %#ok<NASGU>
end
