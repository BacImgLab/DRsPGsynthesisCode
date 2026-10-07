function runDemoPcc(dataRoot, outDir)
%RUNDEMOPCC Analyse the simulated multicolor stack with Colocalization/PCCcaclu.m
%
%   runDemoPcc()                 uses demo/demo_data/PCC and writes into
%                                demo/demo_output/PCC
%   runDemoPcc(dataRoot, outDir) uses the given input and output folders
%
%   The wrapper exists for two reasons:
%
%     - Colocalization/PCCcaclu.m resolves the stage folders relative to the
%       current working directory (parentDir = pwd). A self-contained working
%       copy is therefore assembled in <outDir>/work and the script is executed
%       there, so the repository tree itself is never modified.
%     - PCCcaclu.m starts with "clear; clc;", which would erase the variables of
%       the calling function. It is therefore executed inside its own local
%       function, whose workspace is independent of this one.
%
%   The module script has two parts: the first one collects the original images
%   of each stage into a stageX/PCCcaclu folder (it expects rotated copies in a
%   "Processed" folder, which this demo also provides), the second one computes
%   the Pearson correlation coefficients between the channels and writes a CSV
%   per stage plus the merged Merged_PCC_Values.csv.
%
%   Output written into outDir:
%     Merged_PCC_Values.csv   seven columns: FileName, the three pairwise PCC
%                             values (C2C3, C2C4, C3C4) and the mean intensity
%                             of channels 2, 3 and 4 (meanC2, meanC3, meanC4)
%                             -- direct output of PCCcaclu.m
%     work/                   the isolated working copy used for the run

if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile(fileparts(mfilename('fullpath')), 'demo_data', 'PCC');
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'demo_output', 'PCC');
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

repoDir = fileparts(fileparts(mfilename('fullpath')));

%% --- assemble an isolated working copy (PCCcaclu.m is driven by pwd) ---
workDir = fullfile(outDir, 'work');
if exist(workDir, 'dir')
    rmdir(workDir, 's');
end
mkdir(workDir);

stages = 2:5;
for s = stages
    src = fullfile(dataRoot, sprintf('stage%d', s));
    dst = fullfile(workDir, sprintf('stage%d', s));
    if ~exist(src, 'dir')
        error('runDemoPcc:missingInput', 'Input folder not found: %s', src);
    end
    mkdir(dst);
    copyfile(fullfile(src, 'Ch4_DR_*.tif'), dst);

    % Part 1 of PCCcaclu.m copies the original images into stageX/PCCcaclu,
    % taking them from a "Processed" subfolder that holds rotated copies.
    procDir = fullfile(dst, 'Processed');
    mkdir(procDir);
    fl = dir(fullfile(dst, 'Ch4_DR_*.tif'));
    for k = 1:numel(fl)
        [~, b, ~] = fileparts(fl(k).name);
        copyfile(fullfile(dst, fl(k).name), fullfile(procDir, [b, '_rot.tif']));
    end
end

% The module script must sit next to the stage folders (it uses pwd)
copyfile(fullfile(repoDir, 'Colocalization', 'PCCcaclu.m'), ...
         fullfile(workDir, 'PCCcaclu.m'));

%% --- run the module script in its own workspace ---
oldDir = pwd;
restoreDir = onCleanup(@() cd(oldDir));   %#ok<NASGU>
executePcc(workDir);

%% --- collect the result ---
merged = fullfile(workDir, 'Merged_PCC_Values.csv');
if ~isfile(merged)
    error('runDemoPcc:noOutput', ...
          'PCCcaclu.m did not produce Merged_PCC_Values.csv in %s', workDir);
end
copyfile(merged, fullfile(outDir, 'Merged_PCC_Values.csv'));

T = readtable(merged);
fprintf('PCC demo finished.\n');
fprintf('  stages analysed      : %d (stage2 ... stage5)\n', numel(stages));
fprintf('  images analysed      : %d\n', height(T));
fprintf('  mean PCC ch2-ch3     : %.3f\n', mean(T.C2C3));
fprintf('  mean PCC ch2-ch4     : %.3f\n', mean(T.C2C4));
fprintf('  mean PCC ch3-ch4     : %.3f\n', mean(T.C3C4));

end

function executePcc(workDir)
%EXECUTEPCC Run PCCcaclu.m inside workDir.
%   The script clears its own workspace ("clear; clc;"), so it is kept inside a
%   dedicated function to protect the variables of the caller.

cd(workDir);
PCCcaclu;   %#ok<NASGU>   % module script; resolves stage folders relative to pwd

end
