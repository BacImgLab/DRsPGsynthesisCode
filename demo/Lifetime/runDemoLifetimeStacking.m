function summary = runDemoLifetimeStacking(dataRoot, outDir)
%RUNDEMOLIFETIMESTACKING Stack the LAS X FLIM/FCS export into intensity and lifetime stacks.
%
%   Step 1 of the FLIM (fluorescence lifetime imaging) demo. It drives
%   Section 1 of Lifetime/lifeTcalculationbernsen.m, which reads the intensity
%   and the lifetime image of every exported field of view and writes them into a
%   single intensity stack and a single lifetime stack.
%
%   runDemoLifetimeStacking()            uses demo_data/01_lasx_export
%                                        and writes into demo_output/stacking
%   runDemoLifetimeStacking(dataRoot, outDir)
%
%   Input : demo_data/01_lasx_export/<FOV>/<FOV>-50uM-<n>_ch0.ome.tif   intensity
%           demo_data/01_lasx_export/<FOV>/<FOV>-50uM-<n>_ch1.ome.tif   lifetime
%           (example data 1; <FOV> = field of view, <n> = cell in that field)
%   Output: <outDir>/WT50uMintensityStack.tif            one page per image read
%           <outDir>/WT50uMlifetimeStack.tif
%           <outDir>/Lifetime_stacking_summary.csv
%           <outDir>/work/                               isolated working copy
%
%   Requirements: base MATLAB (imread / imwrite only).
%
%   How the module script is driven. lifeTcalculationbernsen.m is shipped
%   unchanged, and the wrapper replays Section 1 of it in an isolated working
%   copy below <outDir>/work. Three lines that describe the acquisition series
%   the script was written for are replaced by the settings of the demo:
%       dirRoot       -> the working copy
%       fnum          -> the number of images per field of view actually shipped
%       the FOV loop  -> the field-of-view folders actually shipped
%   The script keeps its own file names inside the working copy (it uses a `CFX`
%   prefix for the series it was written for), so the wrapper stages the
%   published images under those names and renames the two stacks afterwards.
%   See the note on the file names in Readme.md.
%
%   The published stacks of example data 2 are used as the reference: every page
%   of both stacks is compared with them and the result is printed.

if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile(fileparts(mfilename('fullpath')), 'demo_data', '01_lasx_export');
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'demo_output', 'stacking');
end
scriptDir = fileparts(mfilename('fullpath'));
refDir    = fullfile(scriptDir, 'demo_data', '02_stacked');

% ---- discover the fields of view, the cells and the channels --------------
d = dir(dataRoot);
d = d([d.isdir]);
fovNames = {d.name};
fovNames = fovNames(~ismember(fovNames, {'.', '..'}));
if isempty(fovNames)
    error('runDemoLifetimeStacking:missingData', 'No field-of-view folder in %s.', dataRoot);
end
if ~all(~cellfun(@isempty, regexp(fovNames, '^\d+$', 'once')))
    error('runDemoLifetimeStacking:dataLayout', ...
          'Field-of-view folders must have numeric names; found: %s', strjoin(fovNames, ', '));
end
fovVals = cellfun(@str2double, fovNames);
[fovVals, ord] = sort(fovVals);
fovNames = fovNames(ord);

cellsPerFov = cell(numel(fovNames), 1);
cellsFirst  = [];
for k = 1:numel(fovNames)
    f     = fovNames{k};
    files = dir(fullfile(dataRoot, f, '*.ome.tif'));
    tok      = regexp({files.name}, '^(\d+)-50uM-(\d+)_ch(\d+)\.ome\.tif$', 'tokens', 'once');
    if numel(files) == 0 || any(cellfun(@isempty, tok))
        error('runDemoLifetimeStacking:dataLayout', ...
              'Expected <FOV>-50uM-<n>_ch<c>.ome.tif in %s.', fullfile(dataRoot, f));
    end
    cellNums = cellfun(@(t) str2double(t{2}), tok);
    chanNums = cellfun(@(t) str2double(t{3}), tok);
    chans    = unique(chanNums(:))';
    if ~isequal(chans, [0 1])
        error('runDemoLifetimeStacking:channels', ...
              'Expected channels 0 (intensity) and 1 (lifetime) in %s, found %s.', ...
              f, mat2str(chans));
    end
    cellsPerFov{k} = unique(cellNums(:))';
    if isempty(cellsFirst)
        cellsFirst = cellsPerFov{k};
    elseif ~isequal(cellsFirst, cellsPerFov{k})
        error('runDemoLifetimeStacking:dataLayout', ...
              'The fields of view do not contain the same cells.');
    end
end
cellsPerFov = cellsPerFov{:};
nCells      = numel(cellsFirst);

fprintf('  field-of-view folders    : %d  %s\n', numel(fovNames), mat2str(fovVals));
fprintf('  images per field of view : %d\n', nCells);
fprintf('  channels                 : 0 = intensity, 1 = lifetime\n');
fprintf('  images to stack          : %d\n', numel(fovNames) * nCells * 2);

% ---- isolated working copy -----------------------------------------------
if exist(outDir, 'dir')
    rmdir(outDir, 's');
end
mkdir(outDir);
workDir = fullfile(outDir, 'work');
mkdir(workDir);

for k = 1:numel(fovNames)
    f = fovNames{k};
    mkdir(fullfile(workDir, f));
    for n = cellsFirst
        for c = [0 1]
            src = fullfile(dataRoot, f, sprintf('%s-50uM-%d_ch%d.ome.tif', f, n, c));
            copyfile(src, fullfile(workDir, f, sprintf('CFX-%d_ch%d.ome.tif', n, c)));
        end
    end
end

% ---- build and run the driver for Section 1 ------------------------------
fprintf('  adapting Section 1 of lifeTcalculationbernsen.m:\n');
lifetimeSectionDriver( ...
    fullfile(scriptDir, 'lifeTcalculationbernsen.m'), ...
    '^%%\s*section\s*1', ...
    workDir, 'run_section1.m', ...
    { '^dirRoot\s*=', sprintf('dirRoot = ''%s'';', workDir), 1; ...
      '^fnum\s*=',    sprintf('fnum = %d; %% number of images per field of view (set by the demo)', nCells), 1; ...
      '^for jj = 1 : 3', sprintf('for jj = %s %% field-of-view folders to stack (set by the demo)', mat2str(fovVals)), 1 });

here = pwd;
cd(workDir);
localRunSection1();       % the driver starts with clear; clc; -> own workspace
cd(here);

% ---- collect and check the two stacks ------------------------------------
pairs = { 'CFXintensityStack.tif', 'WT50uMintensityStack.tif'; ...
          'CFXlifetimeStack.tif',  'WT50uMlifetimeStack.tif'  };
pages     = zeros(1, size(pairs, 1));
pagesSame = zeros(1, size(pairs, 1));
for k = 1:size(pairs, 1)
    srcF = fullfile(workDir, pairs{k, 1});
    if ~isfile(srcF)
        error('runDemoLifetimeStacking:noOutput', ...
              'Section 1 did not write %s.', srcF);
    end
    dstF = fullfile(outDir, pairs{k, 2});
    copyfile(srcF, dstF);

    info      = imfinfo(dstF);
    pages(k)  = numel(info);
    refF      = fullfile(refDir, pairs{k, 2});
    if isfile(refF)
        nRef  = numel(imfinfo(refF));
        same  = 0;
        for p = 1:min(pages(k), nRef)
            if isequal(imread(dstF, p), imread(refF, p))
                same = same + 1;
            end
        end
        pagesSame(k) = same;
    end
    fprintf('  %-24s : %d pages, %d x %d %s\n', pairs{k, 2}, pages(k), ...
            info(1).Height, info(1).Width, class(imread(dstF, 1)));
end
fprintf('  identical to example data 2 : intensity %d of %d pages, lifetime %d of %d pages\n', ...
        pagesSame(1), pages(1), pagesSame(2), pages(2));

% ---- summary table -------------------------------------------------------
q = {}; v = {};
q{end+1} = 'fields_of_view';            v{end+1} = numel(fovNames);
q{end+1} = 'images_per_field_of_view';  v{end+1} = nCells;
q{end+1} = 'images_stacked';            v{end+1} = numel(fovNames) * nCells * 2;
q{end+1} = 'intensity_stack_pages';     v{end+1} = pages(1);
q{end+1} = 'lifetime_stack_pages';      v{end+1} = pages(2);
q{end+1} = 'intensity_pages_identical_to_reference'; v{end+1} = pagesSame(1);
q{end+1} = 'lifetime_pages_identical_to_reference';  v{end+1} = pagesSame(2);
T = table(q(:), v(:), 'VariableNames', {'Quantity', 'Value'});
writetable(T, fullfile(outDir, 'Lifetime_stacking_summary.csv'));

summary = struct('outDir', outDir, ...
                 'fieldsOfView', numel(fovNames), ...
                 'imagesPerField', nCells, ...
                 'intensityStack', fullfile(outDir, pairs{1, 2}), ...
                 'lifetimeStack',  fullfile(outDir, pairs{2, 2}), ...
                 'pagesIdentical', pagesSame);
end

% ---------------------------------------------------------------------------

function localRunSection1()
%LOCALRUNSECTION1 Run the generated driver in a workspace of its own.
%   The driver is Section 1 of the module script, which starts with
%   "clear; clc;" and would therefore wipe the workspace of its caller, so it has
%   to run inside a function of its own.
run_section1;
end
