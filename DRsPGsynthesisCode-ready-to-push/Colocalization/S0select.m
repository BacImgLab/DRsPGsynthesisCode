function [xS0, yS0] = S0select(ImShowing, rotatedBox, filename)
% S0select  ─ Interactive selection of a single septal line (S0)
%
%   [xS0, yS0] = S0select(ImShowing, rotatedBox, filename)
%
%   ImShowing   : 2‑D image to display (typically the membrane channel, already rotated)
%   rotatedBox  : Bounding box of the rotated cell [x, y, width, height]
%   filename    : (optional) name for saving the annotated figure (default: 'temp.fig')
%
%   The function displays the image with contrast stretching and its bounding box,
%   then asks the user to pick two end‑points for the main septum (S0):
%     • S0 (red) – central septum of the dividing cell
%
%   The line is redrawn until the user confirms it. The coordinates of the two
%   selected points are returned. The final annotated figure is saved and closed.
%
%   Last modified: 2025‑02‑17  (Xinxing Yang)

% -------------------------------------------------------------------------
% 0. Default filename if not provided
% -------------------------------------------------------------------------
if nargin < 3
    filename = 'temp.fig';  % Default value for filename if not provided
end

% -------------------------------------------------------------------------
% 1. Display parameters
% -------------------------------------------------------------------------
% scale factor: this factor means how big the image will be showed on the
% screen: 100 is the original size
MagScale = 600;

% the lowP and highP determines how the image will be adjusted in contrast.
% the user can adjust this if the image looks dim or over exposed.
LowP = 5;
HighP = 99;

% -------------------------------------------------------------------------
% 2. Show image with adjusted contrast and bounding box
% -------------------------------------------------------------------------
% the function will output the poisitions of the lines
h = figure;

% Find the minimum and maximum non-zero pixel values
minVal = prctile(double(ImShowing(ImShowing > 0)), LowP);   % 5th percentile
maxVal = prctile(double(ImShowing(ImShowing > 0)), HighP);  % 95th percentile

% Display the image with adjusted contrast using imshow
imshow(ImShowing, [minVal, maxVal], 'InitialMagnification', MagScale);
hold on;

% Draw the bounding box on top
rectangle('Position', rotatedBox, 'EdgeColor', 'y', 'LineWidth', 2);

% -------------------------------------------------------------------------
% 3. Select S0 line (central septum)
% -------------------------------------------------------------------------
% Add title for S0 selection
title('Select the two end points of Septum 0!');
%draw a line to find the profile of S0,S1 left and right
selectionS0Good = false;
xS0 = []; yS0 = [];

% Loop until the user is satisfied with the selection of S0
while ~selectionS0Good
    % initiate some handles for line plot
    S0handle = [];
    S1lefthandle = [];
    S1righthandle = [];

    % Get the two endpoints of the line using ginput
    [xS0, yS0] = ginput(2);  % Get two points from the user

    % Display the selected points on the image
    S0handle = line(xS0, yS0, 'Color', 'r', 'LineWidth', 3);  % Draw a red line between the points

    % Step 4: Ask the user if they are satisfied with the selection
    answer = questdlg('Are you satisfied with the selected points?', ...
                      'Line Selection', ...
                      'Yes', 'No', 'No');

    % If the user is satisfied, break out of the loop
    if strcmp(answer, 'Yes')
        selectionS0Good = true;
    else
        % If not satisfied, clear the plot and prompt to select again
        delete(S0handle)
        disp('Please select the end points of S0 again.');
    end
end

% -------------------------------------------------------------------------
% 4. Save the annotated figure and close
% -------------------------------------------------------------------------
saveas(h, filename);  % Save the figure to file
close(h);             % Close the figure

end
