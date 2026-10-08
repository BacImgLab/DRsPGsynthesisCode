function runDemoColocPcc(dataRoot, outDir)
%RUNDEMOCOLOCPCC Pairwise Pearson correlation of the channels of every cell.
%
%   Step 2 of the multicolor colocalization demo. It drives
%   PCCcaclu.m on the published per-cell output of the septum
%   profiling (example data 4), which is the input listed for the PCC step.
%
%   runDemoColocPcc()            uses demo_data/04_septum_profiles
%                                and writes into demo_output/PCC
%   runDemoColocPcc(dataRoot, outDir)
%
%   Expected input layout below dataRoot (one folder per cell-cycle stage):
%       stage2/Ch4_DR_*.tif                    4-page stack, pages 1..4 = C1..C4
%       stage2/Processed/Ch4_DR_*_rot.tif      same cell, septum rotated vertical
%       ... stage3, stage4, stage5
%
%   PCCcaclu.m resolves its inputs relative to the current folder: part 1 walks
%   stage2 ... stage5, collects the unrotated originals that belong to a rotated
%   copy in the "Processed" subfolder and puts them into stageN/PCCcaclu; part 2
%   computes the Pearson correlation coefficients between channels 2, 3 and 4 on
%   the background-masked pixels, together with the mean intensity of each
%   channel. An isolated working copy is therefore assembled in <outDir>/work so
%   that neither the repository nor the demo data are modified.
%
%   Output written into outDir:
%       Merged_PCC_Values.csv      FileName, C2C3, C2C4, C3C4, meanC2, meanC3, meanC4
%       stageN_PCC_Values.csv      the same table per cell-cycle stage
%       work/                      the isolated working copy used for the run
%
%   Runs on base MATLAB - no toolbox is required.

if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile(fileparts(mfilename('fullpath')), 'demo_data', '04_septum_profiles');
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'demo_output', 'PCC');
end
scriptDir = fileparts(mfilename('fullpath'));

if exist(outDir, 'dir')
    rmdir(outDir, 's');
end
workDir = fullfile(outDir, 'work');
mkdir(workDir);

%% --- assemble an isolated working copy (PCCcaclu.m is driven by pwd) ---
stages = 2:5;
for s = stages
    name = sprintf('stage%d', s);
    src  = fullfile(dataRoot, name);
    if ~exist(src, 'dir')
        error('runDemoColocPcc:missingInput', 'Input folder not found: %s', src);
    end
    copyfile(src, fullfile(workDir, name));   % copies the folder including Processed
end
copyfile(fullfile(scriptDir, 'PCCcaclu.m'), fullfile(workDir, 'PCCcaclu.m'));

%% --- run the module script in its own workspace ---
oldDir = pwd;
pccStep(workDir);
cd(oldDir);

%% --- collect the result ---
% PCCcaclu.m writes each per-stage table into its own stage folder inside the
% working copy, and the merged table into the root of the working copy.
nCells = 0;
nStageTables = 0;
for s = stages
    perStage = fullfile(workDir, sprintf('stage%d', s), ...
                        sprintf('stage%d_PCC_Values.csv', s));
    if isfile(perStage)
        copyfile(perStage, outDir);
        nCells = nCells + height(readtable(perStage));
        nStageTables = nStageTables + 1;
    end
end

merged = fullfile(workDir, 'Merged_PCC_Values.csv');
if ~isfile(merged)
    error('runDemoColocPcc:noOutput', ...
          'PCCcaclu.m did not produce Merged_PCC_Values.csv in %s', workDir);
end
copyfile(merged, outDir);

T = readtable(merged);
fprintf('  cell-cycle stages      : %d (stage2 ... stage5)\n', numel(stages));
fprintf('  per-stage tables       : %d copied next to the merged one\n', nStageTables);
fprintf('  cells analysed         : %d\n', height(T));
fprintf('  mean PCC ch2-ch3       : %.3f\n', mean(T.C2C3));
fprintf('  mean PCC ch2-ch4       : %.3f\n', mean(T.C2C4));
fprintf('  mean PCC ch3-ch4       : %.3f\n', mean(T.C3C4));
fprintf('  output                 : %s\n', outDir);

end

function pccStep(workDir)
%PCCSTEP Run PCCcaclu.m inside its own workspace.
%   The script starts with "clear; clc;", which would wipe the variables of this
%   function if it ran here, so it is called from a nested function.
here = pwd;
cd(workDir);
pccInner();
cd(here);
end

function pccInner()
PCCcaclu; %#ok<NASGU>
end
