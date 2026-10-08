function summary = runDemoSMTVelocity(dataRoot, outDir)
%RUNDEMOSMTVELOCITY Fit the speed distribution of the classified SMT segments.
%
%   Step 1 of the single-molecule-tracking demo. It drives
%   SMTanalysis/dataprocessSMTWCF.m, which fits the empirical cumulative
%   distribution of the directed-segment speeds of FtsW with one and with two
%   log-normal populations, bootstraps both fits (200 resamples) and reconstructs
%   the two probability density functions for the comparison with the histogram.
%
%   runDemoSMTVelocity()                    uses demo_data/04_classified_traces
%                                           and writes into
%                                           demo_output/velocity_distribution
%   runDemoSMTVelocity(dataRoot, outDir)
%
%   Input : demo_data/04_classified_traces/FtsW-all.mat
%           (the classified trajectory data set, example data 4)
%   Output: FtsW-S1-Cef.mat          the complete Result struct of the module script
%           FtsW_speed_CDF.png       CDF of the speeds with the single and the
%                                    double log-normal fit and the residuals
%           FtsW_speed_PDF.png       step histogram with the two reconstructed PDFs
%           SMT_velocity_summary.csv the fitted parameters in a flat two-column table
%
%   Requirements: the Optimization Toolbox (lsqcurvefit) and the Statistics and
%   Machine Learning Toolbox (logncdf, lognpdf, bootstrp).
%
%   Note on the input variable. The module script loads FtsW-all.mat and expects
%   the matrix of per-segment parameters under the name DataSMT, while the
%   published example data set stores exactly the same matrix under the name
%   SegVSRPTC - the name that statesClassifyDr.m (section 4) writes. The
%   published file is left untouched; this wrapper adds the alias DataSMT in its
%   isolated working copy below demo_output/.
%
%   This step is the only step of the SMT workflow that runs without any user
%   interaction, which is why it is the one replayed by the demo. See Readme.md
%   for the interactive steps that come before it.

if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile(fileparts(mfilename('fullpath')), 'demo_data', '04_classified_traces');
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'demo_output', 'velocity_distribution');
end
scriptDir = fileparts(mfilename('fullpath'));

dataFile = fullfile(dataRoot, 'FtsW-all.mat');
if ~isfile(dataFile)
    error('runDemoSMTVelocity:missingInput', 'Input file not found: %s', dataFile);
end

for fn = {'lsqcurvefit', 'logncdf', 'lognpdf', 'bootstrp'}
    if exist(fn{1}, 'file') == 0
        error('runDemoSMTVelocity:missingToolbox', ...
              'This step needs %s (Optimization and Statistics and Machine Learning Toolbox).', fn{1});
    end
end

%% --- isolated working copy -------------------------------------------------
if exist(outDir, 'dir')
    rmdir(outDir, 's');
end
mkdir(outDir);
workDir = fullfile(outDir, 'work');
mkdir(workDir);

copyfile(fullfile(scriptDir, 'dataprocessSMTWCF.m'), workDir);
copyfile(fullfile(scriptDir, 'CDF_logCalc.m'), workDir);
copyfile(fullfile(scriptDir, 'logn1cdf.m'), workDir);
copyfile(fullfile(scriptDir, 'logn2cdf.m'), workDir);

% working copy of the published table, with the variable name the script expects
S = load(dataFile);
S.DataSMT = S.SegVSRPTC;   % alias, see the note in the header
save(fullfile(workDir, 'FtsW-all.mat'), '-struct', 'S');

D       = S.DataSMT;
selMov  = D(:,6) == 3 & D(:,7) == 1;   % directional, region flag 3
selSta  = D(:,6) == 3 & D(:,7) == 3;   % stationary,  region flag 3
nSeg      = size(D, 1);
nMoving   = sum(selMov);
nStation  = sum(selSta);
Vth       = 1;                          % nm/s, as in the module script
Vx        = D(selMov, 1);
nAboveVth = sum(Vx > Vth);

fprintf('  segments in the table     : %d\n', nSeg);
fprintf('  directional segments      : %d (region flag 3)\n', nMoving);
fprintf('  stationary segments       : %d (region flag 3)\n', nStation);
fprintf('  directional, > %g nm/s    : %d  <- used for the CDF\n', Vth, nAboveVth);

%% --- run the module script in the working copy -----------------------------
here   = pwd;
before = findall(groot, 'Type', 'figure');
cd(workDir);
inner();
cd(here);
after  = findall(groot, 'Type', 'figure');

%% --- collect the results ---------------------------------------------------
resFile = fullfile(workDir, 'FtsW-S1-Cef.mat');
if ~isfile(resFile)
    error('runDemoSMTVelocity:noOutput', 'The module script did not write %s.', resFile);
end
copyfile(resFile, outDir);
R = load(resFile);
Result = R.Result;

newFigs = setdiff(after, before);
if ~isempty(newFigs)
    [~, ord] = sort([newFigs.Number]);
    newFigs  = newFigs(ord);
end
figNames = {'FtsW_speed_CDF.png', 'FtsW_speed_PDF.png'};
for k = 1:min(numel(newFigs), numel(figNames))
    exportgraphics(newFigs(k), fullfile(outDir, figNames{k}), 'Resolution', 150);
end
if ~isempty(newFigs)
    close(newFigs);
end

%% --- flat summary table ----------------------------------------------------
q = {}; v = {};
q{end+1} = 'segments_total';              v{end+1} = nSeg;
q{end+1} = 'segments_directional';        v{end+1} = nMoving;
q{end+1} = 'segments_stationary';         v{end+1} = nStation;
q{end+1} = 'speeds_used_for_CDF';         v{end+1} = nAboveVth;
q{end+1} = 'speed_threshold_nm_per_s';    v{end+1} = Vth;

q{end+1} = 'singlefit_P';                 v{end+1} = Result.P_fit1(1);
q{end+1} = 'singlefit_mu';                v{end+1} = Result.P_fit1(2);
q{end+1} = 'singlefit_sigma';             v{end+1} = Result.P_fit1(3);
q{end+1} = 'singlefit_mean_speed_nm_s';   v{end+1} = Result.FIt_Vdirx1(1);
q{end+1} = 'singlefit_sd_speed_nm_s';     v{end+1} = Result.FIt_Vdirx1(2);
q{end+1} = 'singlefit_mean_speed_sem';    v{end+1} = Result.Vdirx1_sem(1);

q{end+1} = 'doublefit_P1';                v{end+1} = Result.P_fit2(1);
q{end+1} = 'doublefit_mu1';               v{end+1} = Result.P_fit2(2);
q{end+1} = 'doublefit_sigma1';            v{end+1} = Result.P_fit2(3);
q{end+1} = 'doublefit_mu2';               v{end+1} = Result.P_fit2(4);
q{end+1} = 'doublefit_sigma2';            v{end+1} = Result.P_fit2(5);
q{end+1} = 'doublefit_frac1_percent';     v{end+1} = Result.Fit_Vdirx2(1);
q{end+1} = 'doublefit_V1_nm_per_s';       v{end+1} = Result.Fit_Vdirx2(2);
q{end+1} = 'doublefit_sd1_nm_per_s';      v{end+1} = Result.Fit_Vdirx2(3);
q{end+1} = 'doublefit_frac2_percent';     v{end+1} = Result.Fit_Vdirx2(4);
q{end+1} = 'doublefit_V2_nm_per_s';       v{end+1} = Result.Fit_Vdirx2(5);
q{end+1} = 'doublefit_sd2_nm_per_s';      v{end+1} = Result.Fit_Vdirx2(6);
q{end+1} = 'doublefit_V1_sem';            v{end+1} = Result.Vdirx2_sem(2);
q{end+1} = 'doublefit_V2_sem';            v{end+1} = Result.Vdirx2_sem(5);

q{end+1} = 'max_abs_residual_single';     v{end+1} = max(abs(Result.Residual1));
q{end+1} = 'max_abs_residual_double';     v{end+1} = max(abs(Result.Residual2));

T = table(q(:), v(:), 'VariableNames', {'Quantity', 'Value'});
writetable(T, fullfile(outDir, 'SMT_velocity_summary.csv'));

fprintf('\n  single log-normal : V = %.2f nm/s (P = %.3f, mu = %.3f, sigma = %.3f)\n', ...
        Result.FIt_Vdirx1(1), Result.P_fit1(1), Result.P_fit1(2), Result.P_fit1(3));
fprintf('  double log-normal : V1 = %.2f nm/s (%.1f %%), V2 = %.2f nm/s (%.1f %%)\n', ...
        Result.Fit_Vdirx2(2), Result.Fit_Vdirx2(1), ...
        Result.Fit_Vdirx2(5), Result.Fit_Vdirx2(4));
fprintf('  residual (max |fit - data|) : %.4f (single) / %.4f (double)\n', ...
        max(abs(Result.Residual1)), max(abs(Result.Residual2)));
fprintf('  figures written   : %d\n', min(numel(newFigs), numel(figNames)));

summary = struct('segments_total', nSeg, ...
                 'segments_directional', nMoving, ...
                 'segments_stationary', nStation, ...
                 'speeds_used_for_CDF', nAboveVth, ...
                 'single_speed', Result.FIt_Vdirx1(1), ...
                 'double_V1', Result.Fit_Vdirx2(2), ...
                 'double_V2', Result.Fit_Vdirx2(5), ...
                 'outDir', outDir);

end

function inner()
%INNER Run the module script in its own workspace.
%   dataprocessSMTWCF.m starts with "clear; clc;", which wipes the workspace of
%   its caller, so it has to run inside a nested function of its own.
dataprocessSMTWCF;
end
