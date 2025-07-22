%% RemoveDrift
%
% Purpose:
%   This script identifies and removes trajectories that belong to drifting 
%   cells (i.e., cells that move significantly, typically due to image drift 
%   or experimental instability). The process is split into three parts:
%
%   1. Visualize all trajectories in the original brightfield (BF) image.
%      These will help the user manually identify drifting cells.
%
%   2. Manually select drifting cells in FIJI/ImageJ using ROIs and save 
%      them as a .roi or .zip file.
%
%   3. Load the saved ROIs and remove all trajectories whose centroid falls 
%      within these regions.
%
% Requirements:
%   - The 'TrackRefine' structure file (MAT-file)
%   - The corresponding BF image (TIFF)
%   - The function `ReadImageJROI.m` to read ImageJ ROI files
%
% Output:
%   - A new 'TrackRefine' file excluding trajectories in the drifting cells
%   - A figure with all original trajectories overlaid for manual inspection
%

%% Section 1: Load trajectory data and plot over BF image
clear
clc

% Load trajectory structure from a .mat file
[filenameTrack, pathnameTrack] = uigetfile('.mat', 'Select the TrackRefine file');
load([pathnameTrack filenameTrack]);

% Load the corresponding brightfield (BF) image
[filenameImg, pathnameImg] = uigetfile('.tif', 'Select the BF image');
img = imread([pathnameImg filenameImg]);
filenameSave = ['traj-' filenameImg];

% Plot all trajectories on top of the BF image
h = imshow(img, []);
hold on
qq = 0; % Counter for trajectories

for ii = 1 : length(TrackRefine)
    if isfield(TrackRefine{ii,1}, 'Seg') % Check if valid trace
        Traj = TrackRefine{ii,1}.TrackOriginal;
        plot(Traj(:,2), Traj(:,3), '-r', 'LineWidth', 3); % x=column 2, y=column 3
        qq = qq + 1;
    end
end

% Save the image for manual inspection in ImageJ/FIJI
saveas(h, [pathnameImg filenameSave]);
close all

%% Section 2: Manual annotation in FIJI/ImageJ
% 1. Open the saved image with trajectories (from Section 1) and the original
%    BF image stack in ImageJ or FIJI.
% 2. Use ROI Manager to identify drifting cells and mark them using ROI tools.
% 3. Save the ROIs as a .roi or .zip file (to be used in Section 3).

%% Section 3: Remove trajectories in drifted (bad) cells
% Load the ROIs created in ImageJ/FIJI
[filenameROI, pathnameROI] = uigetfile({'*.roi;*.zip'}, 'Select the ROI file');
sROI = ReadImageJROI([pathnameROI filenameROI]);

% Ensure sROI is in cell format
if ~iscell(sROI)
    sROI = {sROI};
end

% Initialize new TrackRefine structure without drifted cells
TrackRefineNew = {};
kk = 1; % Index for saving filtered tracks

% Check each trajectory
for ii = 1 : length(TrackRefine)
    Xtraj = TrackRefine{ii,1}; % Get current trajectory structure
    
    if isfield(Xtraj, 'Seg') % Skip empty entries
        Traj = Xtraj.TrackOriginal;

        % Compute trajectory center
        Ycenter = mean(Traj(:,2)); % x-position
        Xcenter = mean(Traj(:,3)); % y-position

        % Check if the center falls inside any ROI
        flagIn = 0; % 1 if inside an ROI
        for jj = 1 : length(sROI)
            Box = sROI{jj}.vnRectBounds; % [Top Left X, Top Left Y, Bottom Right X, Bottom Right Y]
            if Ycenter >= Box(1) && Ycenter <= Box(3) && ...
               Xcenter >= Box(2) && Xcenter <= Box(4)
                flagIn = 1; % Mark as inside ROI
                break;
            end
        end

        % If not in any ROI, retain the trajectory
        if flagIn == 0
            TrackRefineNew{kk,1} = Xtraj;
            kk = kk + 1;
        end
    end
end

% Save the cleaned trajectory set
TrackRefine = TrackRefineNew;
filenameTrackS = ['Driftout-' filenameTrack];
save([pathnameTrack filenameTrackS], 'TrackRefine');

disp('Drifted trajectories removed and saved.');
