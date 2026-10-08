%% tracePlotaftersrg
%
% Description:
% This script visualizes and exports trajectory segment data along with 
% associated bright-field (BF) and fluorescence (FL) ROI images.
% The final output is a multi-page TIFF file showing various trace-specific 
% information for each segmented trajectory.
%
% Input:
% - A .mat file containing a variable `TrackRefine`, which is a cell array
%   of structures. Each structure represents a refined trajectory and includes:
%   * .Seg: segmented trace information
%   * .Coord: full trace coordinates
%   * .ROIBFrot: rotated BF image
%   * .ROIFLrot: rotated FL image
%   * .ROIRotate, .ROIRecRotate: ROI boundary coordinates
%   * .s_MSD: MSD curve data
%   * .TraceId: unique trace identifier
%   * .Pixelsize: pixel-to-nanometer conversion factor
%   * .ExpT: exposure time per frame
%
% Output:
% - A multi-page TIFF file named `test.tif`, where each page contains:
%   1. BF image with overlaid trajectory and ROI
%   2. XY trajectory in real space (nm)
%   3. X vs time plot with velocity fit
%   4. Y vs time plot with velocity fit
%   5. Mean square displacement (MSD) curve
%   6. Fluorescence intensity vs time
%
% Requirements:
% - The external function `colorCodeTracePlot` must be available in the path

clear;
clc;

% Load the trace data structure saved by the traceInROI function
load('TraceRefine-test.mat');

kk = 1; % Counter to track the first figure when saving TIFF

% Loop over all trajectories in the dataset
for ii = 1 : size(TrackRefine, 1)

    % Skip if no segmented data exists for the current trace
    if ~isfield(TrackRefine{ii,1}, 'Seg')
        continue;
    end

    % Extract trajectory data
    Seg = TrackRefine{ii,1}.Seg;
    Ls = length(Seg);
    Trace = TrackRefine{ii,1}.Coord;
    BFimage = TrackRefine{ii,1}.ROIBFrot;
    FLimage = TrackRefine{ii,1}.ROIFLrot;
    ROIRotate = TrackRefine{ii,1}.ROIRotate;
    ROIRecRotate = TrackRefine{ii,1}.ROIRecRotate;
    s_MSD = TrackRefine{ii,1}.s_MSD;
    TraceID = TrackRefine{ii,1}.TraceId;
    PixelS = TrackRefine{ii,1}.Pixelsize;
    ExpT = TrackRefine{ii,1}.ExpT;

    % Loop through all segments of the current trajectory
    for jj = 1 : Ls

        % Initialize a new figure (hidden for batch processing)
        h = figure('visible', 'off');
        set(h, 'position', [100, 100, 1200, 700]);

        % Extract segment-level trace data
        T = Seg(jj).TimeT_TXY;      % Time points
        X = Seg(jj).TraceTx;        % X coordinates
        Y = Seg(jj).TraceTy;        % Y coordinates
        R = Seg(jj).Rxy;            % Correlation coefficient
        Vx = Seg(jj).Vx;            % Velocity in X
        Vy = Seg(jj).Vy;            % Velocity in Y
        V = sqrt(Vx^2 + Vy^2);      % Net velocity
        s_MSD = Seg(jj).MSD;        % MSD curve data
        I = Seg(jj).IntensityValues;% Fluorescence intensity over time

        % --- Plot 1: BF image with trajectory and ROI overlaid ---
        hs1 = subplot('position', [0.03, 0.5, 0.19, 0.4]);
        set(gca, 'XTick', [], 'YTick', []);
        hold on;
        plot(ROIRotate(:,1), ROIRotate(:,2), ':', 'linewidth', 2, 'color', 'y');
        plot(ROIRecRotate(:,1), ROIRecRotate(:,2), '-', 'linewidth', 1.5, 'color', 'c');
        colorCodeTracePlot(T, [X, Y] / PixelS, 0.8);
        set(gca, 'FontSize', 8);
        title(TraceID, 'FontSize', 14);
        axis equal;

        % --- Plot 2: XY trace in real units (nm) ---
        hs2 = subplot('position', [0.03, 0.05, 0.19, 0.4]);
        set(gca, 'XTick', [], 'YTick', []);
        plot(Trace(:,6) * PixelS, Trace(:,7) * PixelS, '-ok', 'markersize', 3);
        hold on;
        colorCodeTracePlot(T, [X, Y], 0.8);
        axis equal;
        title('Positions in nm', 'FontSize', 10);

        % --- Plot 3: X position vs time ---
        hs3 = subplot('position', [0.26, 0.5, 0.51, 0.41]);
        set(gca, 'XTick', [], 'YTick', []);
        plot(Trace(:,1), Trace(:,8), '-ok', 'markersize', 3);
        hold on;
        plot(T, X, '-ro', 'linewidth', 2);
        px = polyfit(T, X, 1);
        plot(T, polyval(px, T), '-b', 'linewidth', 3);
        title(['Trace#' num2str(ii) '    Segment#' num2str(jj) ...
            '         V = ' num2str(V, '%.2f') ' nm/s     R = ' num2str(R, '%.2f')], ...
            'FontSize', 14);
        xlim([min(Trace(:,1)), max(Trace(:,1))]);

        % --- Plot 4: Y position vs time ---
        hs4 = subplot('position', [0.26, 0.05, 0.51, 0.41]);
        set(gca, 'XTick', [], 'YTick', []);
        plot(Trace(:,1), Trace(:,9), '-ok', 'markersize', 3);
        hold on;
        plot(T, Y, '-ro', 'linewidth', 2);
        py = polyfit(T, Y, 1);
        plot(T, polyval(py, T), '-b', 'linewidth', 3);
        title(['Vx = ' num2str(Vx, '%.2f') ' nm/s     Vy = ' num2str(Vy, '%.2f') ' nm/s'], ...
            'FontSize', 10);
        xlim([min(Trace(:,1)), max(Trace(:,1))]);

        % --- Plot 5: MSD curve ---
        hs5 = subplot('position', [0.80, 0.5, 0.17, 0.41]);
        hold on;
        errorbar(s_MSD(:,1), s_MSD(:,2), s_MSD(:,4), 'linewidth', 2);
        set(gca, 'FontSize', 8);
        xlabel('Lag Time (s)', 'FontSize', 9);
        ylabel('MSD', 'FontSize', 9);

        % --- Plot 6: Intensity over time ---
        hs6 = subplot('position', [0.80, 0.05, 0.17, 0.41]);
        hold on;
        plot(T, I, '-ok', 'markersize', 3);
        set(gca, 'FontSize', 8);
        xlabel('Time', 'FontSize', 9);
        ylabel('I', 'FontSize', 9);

        % --- Export current figure as a TIFF page ---
        frame = getframe(h);
        imgplot = frame2im(frame); % Convert to image

        % Save as indexed image for TIFF
        [imind, cm] = rgb2ind(imgplot, 256);

        if kk == 1
            % First image: create new TIFF file
            imwrite(imgplot, 'test.tif', 'tif', 'WriteMode', 'overwrite', 'Compression', 'none');
            kk = 2;
        else
            % Append subsequent images to TIFF
            imwrite(imgplot, 'test.tif', 'tif', 'WriteMode', 'append', 'Compression', 'none');
        end

        % Close the figure to free memory
        close(h);
    end

    ii % Display current trajectory index
end

% Final message
display('All traces saved!');

% Open a blank figure and close it to reset GUI state
figure('visible', 'on');
close all;
