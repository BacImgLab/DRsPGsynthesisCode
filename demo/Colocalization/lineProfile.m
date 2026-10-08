function [grayValues, lineX, lineY, Len] = lineProfile(img, x, y, lineWidth)
% lineProfile - Extract intensity profile along a user-defined line in an image
%
% This function mimics ImageJ's "Plot Profile" feature, allowing users to
% sample intensity values along a line between two selected points (x, y).
% It supports line averaging (smoothing) with a user-defined line width.
%
% INPUTS:
%   img       - 2D grayscale input image
%   x, y      - Vectors of length 2 specifying line endpoints: [x1, x2], [y1, y2]
%   lineWidth - Odd positive integer specifying the width (in pixels) of the sampling region
%
% OUTPUTS:
%   grayValues - Vector of sampled intensity values along the line (averaged if lineWidth > 1)
%   lineX      - x-coordinates along the line (subpixel)
%   lineY      - y-coordinates along the line (subpixel)
%   Len        - Length of the line (Euclidean distance between the two points)
%
% Author: Xinxing Yang
% Date: 2025-02-10

%% Validate lineWidth: must be a positive odd integer
if mod(lineWidth, 2) == 0 || lineWidth <= 0 || mod(lineWidth, 1) ~= 0
    error('Invalid input! The line width must be a positive odd integer.');
end

%% Step 1: Calculate line length (Euclidean distance)
Len = sqrt((diff(x))^2 + (diff(y))^2);

%% Step 2: Round the input coordinates to nearest pixels
x1 = round(x(1)); y1 = round(y(1)); % Starting point
x2 = round(x(2)); y2 = round(y(2)); % Ending point

%% Step 3: Generate evenly spaced subpixel sampling coordinates along the line
numPoints = round(sqrt((x2 - x1)^2 + (y2 - y1)^2));  % Number of samples
lineX = linspace(x1, x2, numPoints);  % X-coordinates along the line
lineY = linspace(y1, y2, numPoints);  % Y-coordinates along the line

%% Step 4: Sample pixel intensities using interpolation
% If lineWidth is 1, directly sample along the center line
if lineWidth == 1
    grayValues = zeros(1, numPoints);
    for i = 1:numPoints
        grayValues(i) = interp2(img, lineX(i), lineY(i), 'linear', 0); % 0 = fill value for out-of-bounds
    end
else
    %% Step 5: Sample a stripe region perpendicular to the line (for averaging)
    
    % Compute unit vector along the main line
    lineVec = [x2 - x1, y2 - y1];
    lineVecNorm = lineVec / norm(lineVec);
    
    % Compute perpendicular unit vector
    perpVec = [-lineVecNorm(2), lineVecNorm(1)];
    
    % Initialize matrices for sampling region coordinates
    halfWidth = (lineWidth - 1) / 2;
    grayRegion = zeros(lineWidth, numPoints); % Store sampled values

    % Loop through each offset line within the line width
    for offset = -halfWidth : halfWidth
        % Shift coordinates perpendicular to the line
        offsetX = lineX + offset * perpVec(1);
        offsetY = lineY + offset * perpVec(2);
        
        % Interpolate pixel intensities along this offset line
        for j = 1:numPoints
            grayRegion(offset + halfWidth + 1, j) = interp2(img, offsetX(j), offsetY(j), 'linear', 0);
        end
    end

    % Average over all offset lines to get the final profile
    grayValues = mean(grayRegion, 1);
end

end
