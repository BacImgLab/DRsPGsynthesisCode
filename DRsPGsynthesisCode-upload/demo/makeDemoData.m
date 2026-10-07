function makeDemoData(outDir)
%MAKEDEMODATA Generate the small simulated datasets used by the demo.
%
%   makeDemoData()        writes the demo data into <repo>/demo/demo_data
%   makeDemoData(outDir)  writes the demo data into the folder outDir
%
%   Dataset 1 - 3D-dSTORM septum width module
%     dSTORM_profile_01.csv ... dSTORM_profile_05.csv
%     Two columns (Distance_nm, GrayValue): the simulated intensity profile of
%     one septum, sampled along a line ROI drawn across the septum on a 2D
%     projection of 3D-dSTORM data. The files are direct input for
%     3D-dSTORM/dSTORMwidthcaclu.m.
%
%   Dataset 2 - Single-molecule tracking (SMT) module
%     SMT_tracks_demo.csv
%     Columns (track, frame, x, y, intensity): a simulated 2D single-molecule
%     tracking dataset, 12 molecules x 200 frames, 0.16 um per pixel and
%     110 ms per frame. The file is direct input for demo/runDemoSMT.m, i.e.
%     for MSDsingle2D.m, CDF_logCalc.m and linfitR.m.
%
%   Dataset 3 - Multicolor colocalization module (Pearson correlation)
%     PCC/stage2/Ch4_DR_001.tif ... PCC/stage5/Ch4_DR_002.tif
%     Eight simulated 4-channel (4-page) 64 x 64 uint16 TIFF stacks. Channels 2,
%     3 and 4 share a common spatial structure with decreasing correlation
%     (ch2-ch3 > ch2-ch4 > ch3-ch4), mimicking a multicolor colocalization
%     dataset. The files are direct input for Colocalization/PCCcaclu.m.
%
%   All datasets are simulated. They contain no experimental measurements and
%   no experimental information. A fixed random seed (rng(2026,'twister')) is
%   used, so the files are reproducible; the script does not require any
%   MATLAB toolbox.

if nargin < 1 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'demo_data');
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ---------------- 1. Simulated dSTORM septum intensity profiles ----------------
% x: distance along the line ROI (nm). y: gray value of the reconstructed
% localization density. Each profile is a broadened peak (septum) on a
% constant background level, with a small deterministic ripple added to mimic
% detector/background structure.
x    = (0:1:150)';          % distance along the line ROI (nm)
bg   = 100;                 % background gray value
wid  = [10 12 14 16 18];    % simulated septum half-widths sigma (nm)

for k = 1:numel(wid)
    amp = 600 + 50*k;                                           % peak amplitude
    xc  = 75 + 2*(k-1);                                         % peak position (nm)
    y   = bg + amp*exp(-0.5*((x-xc)/wid(k)).^2) ...
              + 6*sin(2*pi*x/37) + 3*cos(2*pi*x/13);            % + deterministic ripple
    T = table(x, y, 'VariableNames', {'Distance_nm', 'GrayValue'});
    writetable(T, fullfile(outDir, sprintf('dSTORM_profile_%02d.csv', k)));
end

%% ---------------- 2. Simulated SMT dataset ----------------
rng(2026, 'twister');           % fixed seed -> reproducible demo data
nTr = 12;                       % number of simulated molecules
nFr = 200;                      % number of frames per molecule
px  = 0.16;                     % pixel size (um)
dt  = 0.11;                     % frame interval (s)
D   = [0.005 0.02 0.05];        % diffusion coefficients used (um^2/s)
loc = 0.03;                     % localization uncertainty (um)

track = zeros(nTr*nFr,1); frame = track; X = track; Y = track; I = track;
n = 0;
for t = 1:nTr
    d = D(mod(t-1, numel(D)) + 1);
    s = sqrt(2*d*dt)/px;                          % per-frame step size (pixels)
    x = zeros(nFr,1); y = zeros(nFr,1);
    x(1) = 20 + 20*rand; y(1) = 20 + 20*rand;     % random start position
    for f = 2:nFr
        x(f) = x(f-1) + s*randn;
        y(f) = y(f-1) + s*randn;
    end
    x = x + (loc/px)*randn(nFr,1);                % localization uncertainty
    y = y + (loc/px)*randn(nFr,1);
    idx = n + (1:nFr); n = n + nFr;
    track(idx) = t;
    frame(idx) = (1:nFr)';
    X(idx) = x;
    Y(idx) = y;
    I(idx) = 1000 + 200*randn(nFr,1);
end

T = table(track, frame, X, Y, I, ...
          'VariableNames', {'track', 'frame', 'x', 'y', 'intensity'});
writetable(T, fullfile(outDir, 'SMT_tracks_demo.csv'));

%% ---------------- 3. Simulated multicolor colocalization stacks ----------------
% One 4-channel (4-page) uint16 stack per field of view, in the folder layout
% expected by Colocalization/PCCcaclu.m (stage2 ... stage5, file name
% Ch4_DR_00X.tif). Channel 2 is the reference signal; channels 3 and 4 are
% progressively weaker, noisier copies of it, so that the Pearson correlation
% coefficients come out as ch2-ch3 > ch2-ch4 > ch3-ch4.
rng(2026, 'twister');           % fixed seed -> reproducible demo data
stages = 2:5;                   % stage folders expected by PCCcaclu.m
nImg   = 2;                     % number of fields of view per stage
side   = 64;                    % image side (pixels)

for s = stages
    stageDir = fullfile(outDir, 'PCC', sprintf('stage%d', s));
    if ~exist(stageDir, 'dir')
        mkdir(stageDir);
    end
    for k = 1:nImg
        base = 150 + 30*s;                                              % offset per stage
        Ch1 = uint16(250*rand(side, side));                             % channel 1 (mostly background)
        Ch2 = uint16(base + 800*rand(side, side));                      % channel 2 (reference signal)
        Ch3 = uint16(min(65535, double(Ch2)*0.70 + 150*rand(side, side)));   % channel 3 (strongly correlated)
        Ch4 = uint16(min(65535, double(Ch2)*0.45 + 300*rand(side, side)));   % channel 4 (partly correlated)
        fn = fullfile(stageDir, sprintf('Ch4_DR_%03d.tif', k));
        imwrite(Ch1, fn, 'tif', 'WriteMode', 'overwrite');
        imwrite(Ch2, fn, 'tif', 'WriteMode', 'append');
        imwrite(Ch3, fn, 'tif', 'WriteMode', 'append');
        imwrite(Ch4, fn, 'tif', 'WriteMode', 'append');
    end
end

fprintf('Simulated demo data written to %s\n', outDir);
fprintf('  dSTORM_profile_01..%02d.csv : %d simulated septum profiles\n', ...
        numel(wid), numel(wid));
fprintf('  SMT_tracks_demo.csv         : %d simulated trajectories x %d frames\n', ...
        nTr, nFr);
fprintf('  PCC/stage2..stage5/         : %d simulated 4-channel stacks\n', ...
        numel(stages)*nImg);

end
