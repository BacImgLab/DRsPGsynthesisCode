function runDemo()
%RUNDEMO Run the complete, self-contained demo of the analysis code.
%
%   The demo needs no experimental data, no Fiji/ImageJ and no MATLAB toolbox.
%   It performs four steps:
%
%     1. makeDemoData.m        writes the small simulated datasets into
%                              demo/demo_data/
%     2. dSTORMwidthcaclu.m    (3D-dSTORM module) reads the simulated septum
%                              intensity profiles and writes
%                              demo/demo_output/dSTORM_profiles/Plots/
%                              (AllWidths.csv plus one annotated PNG per profile)
%     3. runDemoSMT.m          analyses the simulated tracking dataset with
%                              MSDsingle2D.m, CDF_logCalc.m and linfitR.m and
%                              writes the MSD, step-length CDF and velocity
%                              results into demo/demo_output/
%     4. runDemoPcc.m          (multicolor colocalization module) runs
%                              Colocalization/PCCcaclu.m on the simulated
%                              4-channel stacks and writes
%                              demo/demo_output/PCC/Merged_PCC_Values.csv
%
%   Usage (MATLAB):
%       >> cd demo
%       >> runDemo
%
%   From a terminal:
%       matlab -batch "cd demo; runDemo"
%
%   The expected output and the expected run time are listed in README.md
%   (section "Demo").

t0 = tic;
status = {};

demoDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(demoDir);
dataDir = fullfile(demoDir, 'demo_data');
outDir  = fullfile(demoDir, 'demo_output');

addpath(fullfile(repoDir, '3D-dSTORM'));
addpath(fullfile(repoDir, 'SMTanalysis'));

%% --- Step 1: simulated demo data ---
fprintf('\n=== Step 1/4: generating the simulated demo data ===\n');
try
    makeDemoData(dataDir);
    status{end+1} = 'Step 1 (demo data)             : OK';        %#ok<*AGROW>
catch ME
    status{end+1} = ['Step 1 (demo data)             : FAILED - ', ME.message];
    rethrow(ME);
end

%% --- Step 2: 3D-dSTORM septum width module ---
fprintf('\n=== Step 2/4: 3D-dSTORM septum width (dSTORMwidthcaclu.m) ===\n');
srcDir = fullfile(outDir, 'dSTORM_profiles');
if ~exist(srcDir, 'dir')
    mkdir(srcDir);
end
copyfile(fullfile(dataDir, 'dSTORM_profile_*.csv'), srcDir);
try
    runDemoDstorm(srcDir);
    status{end+1} = 'Step 2 (3D-dSTORM septum width): OK';
catch ME
    status{end+1} = ['Step 2 (3D-dSTORM septum width): FAILED - ', ME.message];
    rethrow(ME);
end

%% --- Step 3: SMT trajectory analysis ---
fprintf('\n=== Step 3/4: SMT trajectory analysis (runDemoSMT.m) ===\n');
try
    runDemoSMT(fullfile(dataDir, 'SMT_tracks_demo.csv'), outDir);
    status{end+1} = 'Step 3 (SMT trajectory analysis): OK';
catch ME
    status{end+1} = ['Step 3 (SMT trajectory analysis): FAILED - ', ME.message];
    rethrow(ME);
end

%% --- Step 4: multicolor colocalization (PCC) module ---
fprintf('\n=== Step 4/4: multicolor colocalization (runDemoPcc.m) ===\n');
try
    runDemoPcc(fullfile(dataDir, 'PCC'), fullfile(outDir, 'PCC'));
    status{end+1} = 'Step 4 (colocalization PCC)    : OK';
catch ME
    status{end+1} = ['Step 4 (colocalization PCC)    : FAILED - ', ME.message];
    rethrow(ME);
end

%% --- Summary ---
fprintf('\n=== Demo finished in %.1f s ===\n', toc(t0));
for k = 1:numel(status)
    fprintf('  %s\n', status{k});
end
fprintf('Input data : %s\n', dataDir);
fprintf('Output     : %s\n', outDir);

end
