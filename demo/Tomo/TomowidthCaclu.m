%% TomowidthCaclu.m
%
% Purpose:
%   This script processes a Region of Interest (ROI) polygon from an ImageJ ROI file,
%   allows the user to manually select the main axis by clicking two endpoints,
%   then calculates and visualizes:
%     - The middle (center) line of the ROI along the main axis,
%     - The local width of the ROI perpendicular to that middle line at sampled intervals.
%
% Workflow:
%   1. Load ROI file and extract polygon coordinates.
%   2. User selects two endpoints defining the main axis line.
%   3. Sample points along this main axis line at intervals specified by `interD`.
%   4. For each sample point, find intersections of a perpendicular slice with the ROI polygon,
%      then compute the midpoint of the slice intersection (middle line point).
%   5. Compute local tangent directions along the middle line,
%      and recalculate perpendicular slices to accurately measure local ROI widths.
%   6. Visualize ROI, main axis, slices, intersection points, middle line, and widths.
%   7. Save results and figures for further analysis.
%
% Inputs:
%   pixelS - pixel size in nm (scaling factor)
%   interD - sampling spacing in nm
%   M      - half-length of slice lines in nm (used to ensure slices are long enough)
%
% Outputs:
%   A .mat file containing:
%     - Res: Nx4 array of [distance_along_main_axis, local_width, midpoint_x, midpoint_y]
%     - CoordXY: ROI polygon coordinates in nm
%   A .jpg image file of the figure showing all overlays.
%
% Dependencies:
%   - ReadImageJROI (user-supplied function to read ImageJ ROI files)
%   - polyxpoly (MATLAB built-in for polygon-line intersections)
%
% Usage:
%   Run the script, select the ROI file when prompted,
%   click exactly two points on the displayed ROI to define the main axis,
%   then the script performs the analysis automatically.

clear
clc

% Parameters
pixelS = 0.86;   % pixel size in nm
interD = 1;      % sampling spacing in nm
M = 200;         % half-width for slicing lines (in nm)

% --- Load ROI file ---
[filename, pathname] = uigetfile('*.roi', 'Select ROI file');
if isequal(filename,0)
    disp('No file selected. Exiting.');
    return;
end
fullROIpath = fullfile(pathname, filename);
sROI = ReadImageJROI(fullROIpath);

% Convert ROI coordinates to physical units (nm) and close polygon
CoordXY = sROI.mnCoordinates * pixelS;
CoordXY = [CoordXY; CoordXY(1,:)];  % ensure closed polygon
X = CoordXY(:,1);
Y = CoordXY(:,2);

% Plot ROI and prompt user to select main axis endpoints
h = figure;
h.WindowState = 'maximized';
plot(X, Y, 'b-', 'LineWidth', 1); hold on;
title('Click two endpoints on the ROI');
axis equal

% User selects two points defining the main axis
[xp, yp] = ginput(2);
P1 = [xp(1), yp(1)];
P2 = [xp(2), yp(2)];
plot(xp, yp, 'ro-', 'MarkerSize', 8, 'LineWidth', 1.5);

% Draw main axis line
plot([P1(1) P2(1)], [P1(2) P2(2)], 'r-', 'LineWidth', 2);

% Compute sampling points along main axis
D = interD;
v = P2 - P1;      % vector along main axis
L = norm(v);      % length of main axis
dve = v / L;      % unit vector along main axis
pve = [-dve(2), dve(1)];  % perpendicular unit vector

samps = 0 : D : L;  % sample positions along main axis

% Preallocate arrays for midpoints of intersections
midX = nan(numel(samps),1);
midY = nan(numel(samps),1);

for k = 1:numel(samps)
    % Origin of current slice line
    O = P1 + samps(k) * dve;
    A = O + M * pve;
    B = O - M * pve;
    plot([A(1) B(1)], [A(2) B(2)], 'b-');
    
    % Find intersections between slice line and ROI polygon
    [xi, yi] = polyxpoly(X, Y, [A(1) B(1)], [A(2) B(2)]);
    plot(xi, yi, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
    
    % If exactly two intersection points, compute midpoint
    if numel(xi) == 2
        midX(k) = mean(xi);
        midY(k) = mean(yi);
    else
        warning('Slice %d found %d intersections', k, numel(xi));
    end
end

plot(midX, midY, 'md-', 'MarkerSize', 5, 'LineWidth', 1.5);
xlabel('X (nm)');
ylabel('Y (nm)');
title('Midpoints of Perpendicular Slice Crossings');

% Collect midpoint coordinates
midPoints = [midX, midY];

% Compute tangent vectors and widths at each midpoint
numMP = size(midPoints,1);
tangents = nan(numMP, 2);
widths = nan(numMP, 1);

for k = 1:numMP
    % Compute local tangent using finite differences
    if k == 1
        v = midPoints(2,:) - midPoints(1,:);
    elseif k == numMP
        v = midPoints(end,:) - midPoints(end-1,:);
    else
        v = midPoints(k+1,:) - midPoints(k-1,:);
    end
    t = v / norm(v);  % unit tangent vector
    tangents(k,:) = t;
    
    % Perpendicular direction
    p = [-t(2), t(1)];
    
    % Build long perpendicular slice line through midpoint
    O = midPoints(k,:);
    A = O + M * p;
    B = O - M * p;
    plot([A(1) B(1)], [A(2) B(2)], 'c-', 'LineWidth', 1);
    
    % Find intersections of perpendicular slice with ROI polygon
    [xi, yi] = polyxpoly(X, Y, [A(1) B(1)], [A(2) B(2)]);
    if numel(xi) == 2
        plot(xi, yi, 'ko', 'MarkerSize', 6, 'LineWidth', 1.5);
        % Width is distance between intersection points
        widths(k) = hypot(xi(1)-xi(2), yi(1)-yi(2));
        % Draw width segment
        plot([xi(1) xi(2)], [yi(1) yi(2)], 'g-', 'LineWidth', 2);
    else
        widths(k) = NaN;
        warning('Slice %d gave %d intersections', k, numel(xi));
    end
end

% Overlay center line
plot(midPoints(:,1), midPoints(:,2), 'm-', 'LineWidth', 2);

xlabel('X (nm)');
ylabel('Y (nm)');

% Filter out NaNs in widths and corresponding midpoints
midPointsNonnan = midPoints(~isnan(widths), :);
widthsNonnan = widths(~isnan(widths));

% Compute distances along middle line from first midpoint
ref = midPointsNonnan(1, :);
diffs = midPointsNonnan - ref;
distances = sqrt(sum(diffs.^2, 2));

% Organize results as [distance_along_main_axis, width, midpoint_x, midpoint_y]
Res = [distances, widthsNonnan, midPointsNonnan];

% Save results to .mat file
[~, baseName, ~] = fileparts(filename);
matFilename = fullfile(pathname, [baseName, '.mat']);
save(matFilename, 'Res', 'CoordXY');

% Save figure as jpg
figFilename = fullfile(pathname, [baseName, '.jpg']);
saveas(h, figFilename);

%% Example: calculate tangent and perpendicular between two midpoints
n = 5; % example index
v_pair = midPoints(n+1,:) - midPoints(n,:);
dir_pair = v_pair / norm(v_pair);
perp_pair = [-dir_pair(2), dir_pair(1)];
fprintf('Tangent at midpoint %d→%d = [%.3f, %.3f]\n', n, n+1, dir_pair(1), dir_pair(2));
fprintf('Perpendicular direction = [%.3f, %.3f]\n', perp_pair(1), perp_pair(2));

%% Optional: Interpolate contour (not currently used)
% [Xq, Yq] = interpContour(CoordXY(:,1), CoordXY(:,2), 10);
% figure; hold on; axis equal;
% plot(CoordXY(:,1), CoordXY(:,2), 'bo-', 'DisplayName', 'Original');
% plot(Xq, Yq, 'r.-', 'DisplayName', 'Interpolated');
% legend;

function [Xq, Yq] = interpContour(X, Y, N)
% interpContour Upsample a 2D contour by factor N
%
% Inputs:
%   X, Y - column vectors of original contour points
%   N    - integer upsampling factor (>=1)
%
% Outputs:
%   Xq, Yq - upsampled contour points
%
% Note: Treats contour as open by default.

    % Validate inputs
    X = X(:);
    Y = Y(:);
    if N < 1 || round(N) ~= N
        error('N must be a positive integer.');
    end

    M = numel(X);
    totPts = (M-1)*N + 1;
    Xq = zeros(totPts,1);
    Yq = zeros(totPts,1);

    idx = 1;
    for i = 1:(M-1)
        t = linspace(0, 1, N+1).'; % parameter vector
        Xi = (1-t)*X(i) + t*X(i+1);
        Yi = (1-t)*Y(i) + t*Y(i+1);
        Xq(idx:idx+N-1) = Xi(1:end-1);
        Yq(idx:idx+N-1) = Yi(1:end-1);
        idx = idx + N;
    end

    % Add last original point
    Xq(end) = X(end);
    Yq(end) = Y(end);
end
