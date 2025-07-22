function ROIs = roitextRead(filename, smoothsize, splinesize)
%
% PURPOSE:
%   Reads ROI segmentation results from a Cellpose-generated .txt file,
%   converts them to MATLAB format, and performs:
%     1. Smoothing of the ROI contour points using a moving average.
%     2. Rectangular bounding box fitting via PCA (rotated bounding box).
%     3. Spline interpolation to achieve subpixel resolution of the contour.
%
%   The output includes raw and smoothed contours, bounding box coordinates,
%   PCA coefficients, centroid coordinates, and cell length & width.
%
% INPUTS:
%   filename   - Full path to the Cellpose .txt file containing ROI contours.
%   smoothsize - Window size for moving average smoothing (default: 5).
%   splinesize - Number of interpolation points per original outline point
%                for spline fitting (default: 10).
%
% OUTPUT:
%   ROIs - A cell array where each element contains a structure with fields:
%          ROIraw      : Original contour points (Nx2).
%          ROIsubpixel : Smoothed, spline-interpolated contour points.
%          ROIrectan   : Coordinates of rotated bounding rectangle (5x2).
%          PCAcoeff    : PCA rotation matrix (principal axes).
%          CenterXY    : Centroid of the ROI.
%          Cellsize    : [Length, Width] of the bounding box.
%
% AUTHOR: Xinxing Yang
% DATE  : 2024-08-22

%% Prepare: Check the number of input arguments
if nargin < 3
    splinesize = 10;  % Default value for splinesize
end
if nargin < 2
    smoothsize = 5;   % Default value for smoothsize
end

%% Calculation: Read the text file and convert all ROIs to MATLAB format
fid = fopen(filename, 'r'); % Open the file for reading

ROIs = {}; % Initialize cell array to store ROIs

lineIndex = 1;
while ~feof(fid)
    currentLine = fgetl(fid); % Read one line from the file
    
    rawCoords = strsplit(currentLine, ','); % Split line by commas
    
    % Convert strings to numeric and adjust coordinate system (+1)
    numCoords = str2double(rawCoords) + 1; % ImageJ is zero-based, MATLAB is 1-based
    coords = reshape(numCoords, 2, [])'; % Reshape into Nx2 matrix (x,y pairs)
    
    % Define smoothing function: moving average filter
    smoothCoords = @(coords, windowSize) filter(ones(1, windowSize) / windowSize, 1, coords);
    
    % Duplicate coordinates to connect start & end for smoothing/interpolation
    x = [coords(:,1); coords(:,1)]; 
    y = [coords(:,2); coords(:,2)];
    
    % Create parameter vector for interpolation (not directly used here)
    t = 1:length(x);
    
    % Smooth coordinates with moving average
    x_smooth = smoothCoords(x, smoothsize);
    y_smooth = smoothCoords(y, smoothsize);
    
    % Extract middle half of smoothed coordinates to avoid edge artifacts
    x_smooth_final = x_smooth(floor(length(x)/4):floor(length(x)/4)+length(x)/2);
    y_smooth_final = y_smooth(floor(length(x)/4):floor(length(x)/4)+length(x)/2);
    
    % Interpolate the outline with spline to increase resolution (subpixel)
    t_interp = linspace(1, length(x_smooth_final), splinesize * length(x_smooth_final));
    x_interp = spline(1:length(x_smooth_final), x_smooth_final, t_interp);
    y_interp = spline(1:length(x_smooth_final), y_smooth_final, t_interp);
    
    % Final subpixel-smoothed outline coordinates
    coords_sp = [x_interp', y_interp'];
    
    % Calculate rotated bounding rectangle and other properties via PCA
    [coords_rec, coeff, coordCenter, Length, Width] = calcRotatedBoundingBox(coords);
    
    % Store data in output structure
    ROIs{lineIndex, 1}.ROIraw      = coords;
    ROIs{lineIndex, 1}.ROIsubpixel = coords_sp;
    ROIs{lineIndex, 1}.ROIrectan   = coords_rec;
    ROIs{lineIndex, 1}.PCAcoeff    = coeff;
    ROIs{lineIndex, 1}.CenterXY    = coordCenter;
    ROIs{lineIndex, 1}.Cellsize    = [Length, Width];
    
    lineIndex = lineIndex + 1;
    
    % Optional plotting for checking (commented out)
    % plot(x_smooth_final, y_smooth_final, '-o', 'LineWidth', 1, 'Marker', 'o');
    % hold on;
    % plot(x_interp, y_interp, 'b', 'LineWidth', 1);
    % plot(coords_rec(:,1), coords_rec(:,2), 'g--', 'LineWidth', 1);
end

fclose(fid); % Close the file

%% Subfunction to calculate rotated bounding box using PCA
function [ROIrec, coeff, coordCenter, Length, Width] = calcRotatedBoundingBox(coords)
    % Center the coordinates
    meanCoords = mean(coords);
    coordCenter = meanCoords; % centroid
    
    centeredCoords = coords - meanCoords;
    
    % Perform PCA to find principal axes
    [coeff, ~, ~] = pca(centeredCoords);
    
    % Rotate points to align with principal axes
    rotatedCoords = centeredCoords * coeff;
    
    % Find min/max along principal axes
    minX = min(rotatedCoords(:,1));
    maxX = max(rotatedCoords(:,1));
    minY = min(rotatedCoords(:,2));
    maxY = max(rotatedCoords(:,2));
    
    Length = maxX - minX;
    Width = maxY - minY;
    
    % Define the rectangle in rotated frame
    rect = [
        minX, minY;
        maxX, minY;
        maxX, maxY;
        minX, maxY;
        minX, minY
    ];
    
    % Transform rectangle back to original coordinate system
    ROIrec = rect * coeff' + meanCoords;
end

end
