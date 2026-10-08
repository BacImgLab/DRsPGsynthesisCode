function [angle, rotatedBox] = bf2rotateRec(bfimg)
% bf2rotateRec - Estimate the orientation of the largest object in a BF image 
%                and compute its rotated bounding box.
%
% This function performs the following operations on a bright-field (BF) image:
%   1. Threshold the input image (implicitly by treating 0 as background).
%   2. Find all connected components (ROIs).
%   3. Select the largest component based on area.
%   4. Measure its orientation (angle between the major axis and the x-axis).
%   5. Rotate the image to align the major axis horizontally.
%   6. Compute the bounding box of the rotated ROI.
%
% INPUT:
%   bfimg - A binary or segmented bright-field image containing one or more ROIs.
% OUTPUT:
%   angle - The orientation angle (in degrees) of the major axis with respect to x-axis.
%   rotatedBox - The bounding box [x, y, width, height] of the largest object after rotation.
% NOTE:
%   The input image is expected to be a single segmented DR-stage cell image (e.g., stage 3–5).
%   Assumes the long axis of the cell represents the main orientation.
%
% By Xinxing Yang, 2025-02-02

% Step 1: Label connected components (8-connected neighborhood)
cc = bwconncomp(bfimg, 8);

% Step 2: Compute region properties: orientation, bounding box, area, centroid
stats = regionprops(cc, 'Orientation', 'BoundingBox', 'Area', 'Centroid');

% Step 3: Find the largest connected component based on area
areas = [stats.Area];
[~, idx] = max(areas);
largestRegion = stats(idx);

% Step 4: Get the orientation of the largest region (degrees, counterclockwise from x-axis)
angle = largestRegion.Orientation;

% Step 5: Rotate the input image to align the major axis horizontally
% imrotate uses counterclockwise rotation, so we rotate by -angle
rotatedMask = imrotate(bfimg, -angle, 'bilinear', 'crop');

% Step 6: Find the bounding box of the rotated image
cc_rotated = bwconncomp(rotatedMask);
stats_rotated = regionprops(cc_rotated, 'BoundingBox');

% Get the bounding box of the first (usually largest) component after rotation
% This assumes the largest component still appears first after rotation
rotatedBox = stats_rotated(1).BoundingBox;

% % Optional Step 7: Visualization (uncomment to view result)
% figure;
% imshow(rotatedMask);
% hold on;
% rectangle('Position', rotatedBox, 'EdgeColor', 'r', 'LineWidth', 2);
% title('Rotated Mask with Bounding Box');
% hold off;

end
