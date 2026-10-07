function runDemoSMT(dataFile, outDir)
%RUNDEMOSMT Analyse the simulated SMT dataset with the functions of SMTanalysis/.
%
%   runDemoSMT()                  uses demo/demo_data/SMT_tracks_demo.csv and
%                                 writes the results into demo/demo_output
%   runDemoSMT(dataFile, outDir)  uses the given input file and output folder
%
%   The function demonstrates three numerical analysis functions of this
%   repository on the small simulated tracking dataset:
%
%     MSDsingle2D.m   mean squared displacement of every trajectory
%     CDF_logCalc.m   cumulative distribution of the frame-to-frame step length
%     linfitR.m       linear fit of the position vs. time of every trajectory
%
%   Output files written into outDir:
%     SMT_MSD_demo.csv        lag time (s), mean MSD (um^2), SEM (um^2)
%     SMT_stepCDF_demo.csv    step length (um), cumulative probability
%     SMT_velocity_demo.csv   per-trajectory fitted velocity (um/s)
%     SMT_MSD_demo.png        log-log MSD curve
%     SMT_stepCDF_demo.png    step-length CDF
%
%   The acquisition parameters of the simulated data are 0.16 um per pixel and
%   110 ms per frame. Only base MATLAB is required.

if nargin < 1 || isempty(dataFile)
    dataFile = fullfile(fileparts(mfilename('fullpath')), 'demo_data', 'SMT_tracks_demo.csv');
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'demo_output');
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

px   = 0.16;    % pixel size (um)
dt   = 0.11;    % frame interval (s)
dark = 0;       % dark time between two frames (s)
Nlag = 50;      % number of MSD lag times used by MSDsingle2D

T = readtable(dataFile);
trk = unique(T.track);

MSDall = [];  steps = [];
vx = zeros(numel(trk),1); vy = vx; speed = vx;
displ = vx; resid = vx;

for k = 1:numel(trk)
    t = sortrows(T(T.track == trk(k), :), 'frame');

    % --- MSD of this trajectory (repository function) ---
    trace = [t.frame, t.x, t.y];                 % [frame, x(px), y(px)]
    m = MSDsingle2D(trace, px, dt, dark, Nlag);  % col 1: lag time, col 2: MSD (um^2)
    if isempty(MSDall)
        MSDall = nan(numel(trk), size(m,1));
        lagT   = m(:,1);
    end
    MSDall(k,:) = m(:,2)';

    % --- frame-to-frame step lengths (um) ---
    steps = [steps; px*sqrt(diff(t.x).^2 + diff(t.y).^2)];

    % --- linear fit of position vs. time (repository function) ---
    tm = t.frame*dt;
    [~, px_fit, dx, stx] = linfitR(tm, t.x*px);
    [~, py_fit, dy, sty] = linfitR(tm, t.y*px);
    vx(k) = px_fit(1);  vy(k) = py_fit(1);
    speed(k) = sqrt(vx(k)^2 + vy(k)^2);
    displ(k) = sqrt(dx^2 + dy^2);
    resid(k) = (stx + sty)/2;
end

meanMSD = mean(MSDall, 1, 'omitnan')';
semMSD  = std(MSDall, 0, 1, 'omitnan')' / sqrt(size(MSDall,1));

xbin = 0:0.005:0.5;                 % step-length bins (um)
CDF  = CDF_logCalc(steps, xbin);    % repository function

MSDtab = table(lagT(:), meanMSD, semMSD, ...
               'VariableNames', {'LagTime_s', 'MeanMSD_um2', 'SEM_um2'});
writetable(MSDtab, fullfile(outDir, 'SMT_MSD_demo.csv'));

CDFtab = table(xbin(:), CDF(:), ...
               'VariableNames', {'StepLength_um', 'CumulativeProbability'});
writetable(CDFtab, fullfile(outDir, 'SMT_stepCDF_demo.csv'));

Vtab = table(trk(:), vx, vy, speed, displ, resid, ...
             'VariableNames', {'track', 'vx_um_per_s', 'vy_um_per_s', ...
                               'speed_um_per_s', 'displacement_um', 'fitResidual_um'});
writetable(Vtab, fullfile(outDir, 'SMT_velocity_demo.csv'));

f1 = figure('Visible','off','Position',[100 100 560 460]);
loglog(lagT, meanMSD, 'o-', 'LineWidth', 1.4, 'MarkerSize', 4); hold on;
loglog(lagT, lagT*mean(meanMSD(1:10)./lagT(1:10)), 'k--', 'LineWidth', 1);
xlabel('Lag time (s)'); ylabel('MSD (\mum^2)');
title(sprintf('Simulated SMT data: MSD (mean of %d trajectories)', numel(trk)));
legend({'mean MSD \pm SEM', 'linear reference'}, 'Location','northwest');
grid on; box on;
saveas(f1, fullfile(outDir, 'SMT_MSD_demo.png')); close(f1);

f2 = figure('Visible','off','Position',[100 100 560 460]);
plot(xbin, CDF, 'LineWidth', 1.6);
xlabel('Step length (\mum)'); ylabel('Cumulative probability');
title('Simulated SMT data: cumulative distribution of step lengths');
grid on; box on;
saveas(f2, fullfile(outDir, 'SMT_stepCDF_demo.png')); close(f2);

fprintf('SMT demo finished.\n');
fprintf('  trajectories analysed       : %d\n', numel(trk));
fprintf('  median step length          : %.3f um\n', median(steps));
fprintf('  steps with length <= 0.5 um : %.1f %%\n', 100*CDF(end));
fprintf('  mean fitted speed           : %.4f um/s\n', mean(speed));
fprintf('  MSD at the last lag (%4.2f s): %.4f um^2\n', lagT(end), meanMSD(end));

end
