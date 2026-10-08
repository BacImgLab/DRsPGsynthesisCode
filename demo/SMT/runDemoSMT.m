function runDemoSMT()
%RUNDEMOSMT Run the unattended step of the single-molecule-tracking module.
%
%   The single-molecule-tracking work flow (Single_Molecule_Tracking) consists of
%   four parts.
%
%   The upstream part cannot be replayed unattended, because it needs Fiji/ImageJ,
%   Cellpose3 or a human: the ROI cropping (SMTdataPrepare.m), the ThunderSTORM
%   localisation (MacroThurderSTORMDrJ.ijm), the chromatic-aberration correction,
%   the Cellpose3 segmentation, the trajectory linking (spotsLinking.m, which
%   opens file dialogs), the interactive trajectory segmentation in the
%   RefineTraceSegDr App Designer app, and the state classification
%   (statesClassifyDr.m, which opens save dialogs and waits for a button press).
%
%   Only the last analysis step runs without any user interaction, and that is
%   what this demo executes on the published example data set:
%
%     Step 1 : SMTanalysis/dataprocessSMTWCF.m
%              Speed-distribution fitting of the directed segments of FtsW. The
%              empirical CDF of the speeds is fitted with a single and with a
%              double log-normal population, both fits are bootstrapped and the
%              two PDFs are reconstructed for the histogram comparison.
%              Input  : example data 4 (classified trajectory data, FtsW-all.mat)
%              Output : FtsW-S1-Cef.mat, SMT_velocity_summary.csv and the two
%                       figures FtsW_speed_CDF.png / FtsW_speed_PDF.png
%
%   This is a self-contained package: the demo wrapper, the module scripts it
%   drives, the published example data set (demo_data/) and the module
%   documentation (Readme.md) all live in this folder, so the folder can be
%   downloaded and run on its own.
%
%   Usage (MATLAB):
%       >> cd demo/SMT
%       >> runDemoSMT
%
%   From a terminal:
%       matlab -batch "cd demo/SMT; runDemoSMT"
%
%   Requirements: the Optimization Toolbox (lsqcurvefit) and the Statistics and
%   Machine Learning Toolbox (logncdf, lognpdf, bootstrp). Neither Fiji/ImageJ
%   nor the experimental raw data are needed.
%
%   All output goes to demo/SMT/demo_output/, which is not tracked by git; the
%   package and the demo data are never modified.
%
%   Expected output and run time are listed in Readme.md (section "Automated
%   demo") and in the Demo section of the repository README.md.

t0 = tic;
status = {};

demoDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(demoDir, 'demo_data');
outRoot = fullfile(demoDir, 'demo_output');
if exist(outRoot, 'dir')
    rmdir(outRoot, 's');
end
mkdir(outRoot);

%% --- Step 1: speed-distribution fitting ---
fprintf('\n=== Step 1/1: speed-distribution fitting (runDemoSMTVelocity.m) ===\n');
try
    runDemoSMTVelocity(fullfile(dataDir, '04_classified_traces'), ...
                       fullfile(outRoot, 'velocity_distribution'));
    status{end+1} = 'Step 1 (speed-distribution fitting) : OK';
catch ME
    status{end+1} = ['Step 1 (speed-distribution fitting) : FAILED - ', ME.message];
    rethrow(ME);
end

%% --- Summary ---
fprintf('\n=== SMT demo finished in %.1f s ===\n', toc(t0));
for k = 1:numel(status)
    fprintf('  %s\n', status{k});
end
fprintf('Input data : %s\n', dataDir);
fprintf('Output     : %s\n', outRoot);

end
