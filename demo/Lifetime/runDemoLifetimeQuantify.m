function summary = runDemoLifetimeQuantify(stackDir, maskDir, outDir)
%RUNDEMOLIFETIMEQUANTIFY Per-pixel lifetime and intensity quantification of the masked stacks.
%
%   Step 2 of the FLIM (fluorescence lifetime imaging) demo. It drives
%   Section 3 of Lifetime/lifeTcalculationbernsen.m, which reads the intensity
%   stack, the lifetime stack and the Bernsen binary mask, sets every background
%   pixel to zero in both stacks and computes the lifetime and the intensity of
%   every remaining pixel, together with the lifetime distribution.
%
%   runDemoLifetimeQuantify()                     uses demo_data/02_stacked and
%                                                 demo_data/03_background_masked
%                                                 and writes into
%                                                 demo_output/quantification
%   runDemoLifetimeQuantify(stackDir, maskDir, outDir)
%
%   Input : <stackDir>/WT50uMintensityStack.tif            6 pages, uint16
%           <stackDir>/WT50uMlifetimeStack.tif             6 pages, uint16
%           <maskDir>/WT50uMintensityStack-binary.tif      6 pages, uint8, 0/255
%           (the stacks are example data 2, the mask is example data 3 - the mask
%           is the output of the Fiji plugin `bersenThtest`, which cannot be run
%           from MATLAB; see Readme.md)
%   Output: <outDir>/WT50uMintensityStack-filter.tif   background set to zero
%           <outDir>/WT50uMlifetimeStack-filter.tif
%           <outDir>/LifetimeResults.mat               hLT, Result, LT_all,
%                                                      LT_mean, In_mean, In_all
%           <outDir>/Lifetime_histogram.png            lifetime histogram
%           <outDir>/Lifetime_summary.csv              the numbers printed below
%           <outDir>/Lifetime_per_page.csv             pixel count and means per page
%           <outDir>/work/                             isolated working copy
%
%   Requirements: base MATLAB (imread / imwrite / imquantize / histcounts / bar).
%
%   How the module script is driven. lifeTcalculationbernsen.m is shipped
%   unchanged, and the wrapper replays Section 3 of it in an isolated working
%   copy below <outDir>/work. Two lines are replaced:
%       dirRoot   -> the working copy
%       filenameS -> a fixed file name (the shipped script asks for the output
%                    `.mat` with uiputfile, i.e. it waits for a click on Save)
%   The script keeps its own file names inside the working copy (it uses a `CFX`
%   prefix for the series it was written for), so the wrapper stages the
%   published stacks and the mask under those names and renames the results
%   afterwards. See the note on the file names in Readme.md.

if nargin < 1 || isempty(stackDir)
    stackDir = fullfile(fileparts(mfilename('fullpath')), 'demo_data', '02_stacked');
end
if nargin < 2 || isempty(maskDir)
    maskDir = fullfile(fileparts(mfilename('fullpath')), 'demo_data', '03_background_masked');
end
if nargin < 3 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'demo_output', 'quantification');
end
scriptDir = fileparts(mfilename('fullpath'));

% ---- the file names of the published example data -------------------------
fIntensity = 'WT50uMintensityStack.tif';
fLifetime  = 'WT50uMlifetimeStack.tif';
fBinary    = 'WT50uMintensityStack-binary.tif';
fInFilter  = 'WT50uMintensityStack-filter.tif';
fLtFilter  = 'WT50uMlifetimeStack-filter.tif';

inI = fullfile(stackDir, fIntensity);
inL = fullfile(stackDir, fLifetime);
inB = fullfile(maskDir,  fBinary);
for f = {inI, inL, inB}
    if ~isfile(f{1})
        error('runDemoLifetimeQuantify:missingInput', 'Input file not found: %s', f{1});
    end
end
nStack = numel(imfinfo(inI));
nMask  = numel(imfinfo(inB));
if nStack ~= nMask
    error('runDemoLifetimeQuantify:sizeMismatch', ...
          'The stacks have %d pages and the mask has %d.', nStack, nMask);
end
fprintf('  intensity stack          : %s (%d pages)\n', fIntensity, nStack);
fprintf('  lifetime stack           : %s (%d pages)\n', fLifetime, numel(imfinfo(inL)));
fprintf('  binary mask              : %s (%d pages)\n', fBinary, nMask);

% ---- isolated working copy -----------------------------------------------
if exist(outDir, 'dir')
    rmdir(outDir, 's');
end
mkdir(outDir);
workDir = fullfile(outDir, 'work');
mkdir(workDir);

% the module script resolves its inputs by name, with the prefix of the series
% it was written for
copyfile(inI, fullfile(workDir, 'CFXintensityStack.tif'));
copyfile(inL, fullfile(workDir, 'CFXlifetimeStack.tif'));
copyfile(inB, fullfile(workDir, 'CFXintensityStack-binary.tif'));

% ---- build and run the driver for Section 3 ------------------------------
fprintf('  adapting Section 3 of lifeTcalculationbernsen.m:\n');
lifetimeSectionDriver( ...
    fullfile(scriptDir, 'lifeTcalculationbernsen.m'), ...
    '^%%\s*section\s*3', ...
    workDir, 'run_section3.m', ...
    { '^dirRoot\s*=', sprintf('dirRoot = ''%s'';', workDir), 1; ...
      '^filenameS\s*=', 'filenameS = ''LifetimeResults.mat''; % fixed name: the demo must not wait for a save dialog', 1 });

here   = pwd;
before = findall(groot, 'Type', 'figure');
cd(workDir);
localRunSection3();
cd(here);
after  = findall(groot, 'Type', 'figure');

% ---- collect the results -------------------------------------------------
matFile = fullfile(workDir, 'LifetimeResults.mat');
if ~isfile(matFile)
    error('runDemoLifetimeQuantify:noOutput', 'Section 3 did not write %s.', matFile);
end
copyfile(matFile, fullfile(outDir, 'LifetimeResults.mat'));
R = load(matFile);

newFigs = setdiff(after, before);
if ~isempty(newFigs)
    [~, ord] = sort([newFigs.Number]);
    newFigs  = newFigs(ord);
end
figWritten = 0;
if ~isempty(newFigs)
    exportgraphics(newFigs(1), fullfile(outDir, 'Lifetime_histogram.png'), 'Resolution', 150);
    figWritten = 1;
end
if ~isempty(newFigs)
    close(newFigs);
end

pairs = { 'CFXMintensityStack-filter.tif',  fInFilter; ...
          'CFXXMlifetimeStack-filter.tif',  fLtFilter };
pagesSame = zeros(1, 2);
for k = 1:2
    srcF = fullfile(workDir, pairs{k, 1});
    if ~isfile(srcF)
        error('runDemoLifetimeQuantify:noOutput', 'Section 3 did not write %s.', srcF);
    end
    dstF = fullfile(outDir, pairs{k, 2});
    copyfile(srcF, dstF);

    refF = fullfile(maskDir, pairs{k, 2});
    if isfile(refF)
        nRef = numel(imfinfo(refF));
        same = 0;
        for p = 1:min(numel(imfinfo(dstF)), nRef)
            if isequal(imread(dstF, p), imread(refF, p))
                same = same + 1;
            end
        end
        pagesSame(k) = same;
    end
end

LT_all  = R.LT_all;
In_all  = R.In_all;
LT_mean = R.LT_mean;
In_mean = R.In_mean;
hLT     = R.hLT;
Result  = R.Result;

nPix       = numel(LT_all);
perPage    = arrayfun(@(r) numel(r.LT_sample), Result);
ltPerPage  = arrayfun(@(r) mean(r.LT_sample), Result);
inPerPage  = arrayfun(@(r) mean(r.In_sample), Result);
[pk, ipk]  = max(hLT(:, 2));

fprintf('  images analysed          : %d\n', numel(Result));
fprintf('  signal pixels            : %d\n', nPix);
fprintf('  lifetime  (ns)           : mean %.4f  sd %.4f  range %.4f - %.4f\n', ...
        LT_mean(1), LT_mean(2), min(LT_all), max(LT_all));
fprintf('  intensity (a.u.)         : mean %.2f  sd %.2f  range %g - %g\n', ...
        In_mean(1), In_mean(2), min(In_all), max(In_all));
fprintf('  lifetime histogram       : %d bins, %.2f - %.2f ns, peak %.4f at %.2f ns\n', ...
        size(hLT, 1), hLT(1, 1), hLT(end, 1), pk, hLT(ipk, 1));
fprintf('  pixels with NaN lifetime : %d\n', sum(isnan(LT_all)));
fprintf('  identical to example data 3 : intensity %d of %d pages, lifetime %d of %d pages\n', ...
        pagesSame(1), nStack, pagesSame(2), nStack);
fprintf('  figures written          : %d\n', figWritten);

% ---- summary tables ------------------------------------------------------
q = {}; v = {};
q{end+1} = 'images_analysed';        v{end+1} = numel(Result);
q{end+1} = 'signal_pixels';          v{end+1} = nPix;
q{end+1} = 'lifetime_mean_ns';       v{end+1} = LT_mean(1);
q{end+1} = 'lifetime_sd_ns';         v{end+1} = LT_mean(2);
q{end+1} = 'lifetime_min_ns';        v{end+1} = min(LT_all);
q{end+1} = 'lifetime_max_ns';        v{end+1} = max(LT_all);
q{end+1} = 'intensity_mean';         v{end+1} = In_mean(1);
q{end+1} = 'intensity_sd';           v{end+1} = In_mean(2);
q{end+1} = 'intensity_min';          v{end+1} = min(In_all);
q{end+1} = 'intensity_max';          v{end+1} = max(In_all);
q{end+1} = 'histogram_bins';         v{end+1} = size(hLT, 1);
q{end+1} = 'histogram_first_ns';     v{end+1} = hLT(1, 1);
q{end+1} = 'histogram_last_ns';      v{end+1} = hLT(end, 1);
q{end+1} = 'histogram_peak_ns';      v{end+1} = hLT(ipk, 1);
q{end+1} = 'histogram_peak_prob';    v{end+1} = pk;
q{end+1} = 'pixels_with_nan';        v{end+1} = sum(isnan(LT_all));
q{end+1} = 'intensity_pages_identical_to_reference'; v{end+1} = pagesSame(1);
q{end+1} = 'lifetime_pages_identical_to_reference';  v{end+1} = pagesSame(2);
T = table(q(:), v(:), 'VariableNames', {'Quantity', 'Value'});
writetable(T, fullfile(outDir, 'Lifetime_summary.csv'));

Tp = table((1:numel(Result))', perPage(:), ltPerPage(:), inPerPage(:), ...
           'VariableNames', {'Page', 'SignalPixels', 'MeanLifetime_ns', 'MeanIntensity'});
writetable(Tp, fullfile(outDir, 'Lifetime_per_page.csv'));

summary = struct('outDir', outDir, ...
                 'imagesAnalysed', numel(Result), ...
                 'signalPixels', nPix, ...
                 'lifetimeMean', LT_mean(1), ...
                 'lifetimeSd', LT_mean(2), ...
                 'intensityMean', In_mean(1), ...
                 'intensitySd', In_mean(2), ...
                 'histogramPeak', pk, ...
                 'pagesIdentical', pagesSame);
end

% ---------------------------------------------------------------------------

function localRunSection3()
%LOCALRUNSECTION3 Run the generated driver in a workspace of its own.
%   Section 3 also opens the lifetime histogram figure and assigns the variables
%   the wrapper reads back from LifetimeResults.mat, so it is kept out of the
%   wrapper's workspace.
run_section3;
end
