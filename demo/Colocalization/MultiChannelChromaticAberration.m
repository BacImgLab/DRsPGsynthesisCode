%% Multichannel Chromatic Aberration Correction after Drift Correction
%
% This script performs chromatic aberration correction for multichannel images.
% It assumes the images have already undergone drift correction.
% channel1: BF; channel2:647; channel3:561; channel4:488. 
%
% [Processing Steps]:
%   1. Load beads images from multiple positions for 3 fluorescence channels (488, 561, 647 nm);
%   2. Select region of interest (ROI) on maximum projection of the 647 nm channel from beads images;
%   3. Crop and stack the ROI from all bead images across positions;
%   4. Compute spatial transformations to align 488 nm and 647 nm channels to the 561 nm reference channel;
%   5. Apply transformations to experimental datasets for chromatic shift correction;
%   6. Save corrected image stacks and transformation parameters.
%
% [Requirements]:
%   - User interaction for ROI selection.
%
% [Output]:
%   - ROI stacks of beads: roi488stack.tif, roi561stack.tif, roi647stack.tif;
%   - Maximum projections of beads: roi488max.tif, roi561max.tif, roi647max.tif;
%   - Aligned images: roi488aligned.tif, roi647aligned.tif;
%   - Transformation matrices saved in tforms.mat;
%   - Corrected experimental images saved as multi-page TIFFs.
% Cite:
% MATLAB version: R2024a (recommended)

clear; clc;

%% Section0: Set up parameters for bead image processing
clear; clc;
Position = 100; % Number of bead positions (folders) to process
ROIsize = 1600; % Size of the square ROI in pixels
Edgesize = 10;  % Pixels to trim from edges to remove calibration edge effects
beadFolder = 'E:\Dr strain\20250518\beads_0\';  % Folder containing bead images, Replace with the actual path to your dataset

%% Section1: Find the ROI from the 647 channel of beads
% This section loads all the 647 channel images across different positions of beads,
% creates a maximum intensity projection, and allows the user to select a square region of interest (ROI) interactively.
% The selected ROI coordinates are saved for later cropping.

for ii = 1 : Position
    foldname = [beadFolder 'Pos' num2str(ii) '\img_000000000_Channel2_000.tif'];
    Im(:,:,ii) = imread(foldname);
end

% Maximum intensity projection across the Z-stack (positions)
ImMax = max(Im, [], 3);

roijudge = 0;  % Flag for user confirmation of ROI selection
while roijudge == 0
    % Display the max projection with adjusted contrast
    h = imshow(ImMax, [0.1*max(ImMax(:)), 0.4*max(ImMax(:))]);
    colormap(cool)
    
    % User selects a point which will be the center of ROI
    [x, y] = ginput(1);
    
    % Define the square ROI centered at the selected point
    ROI = [round(x) - round(ROIsize/2), round(y) - round(ROIsize/2), ROIsize, ROIsize]; % [x, y, width, height]
    
    hold on
    % Draw the ROI rectangle on the image
    plot([ROI(1), ROI(1) + ROI(3), ROI(1) + ROI(3), ROI(1), ROI(1)], ...
         [ROI(2), ROI(2), ROI(2) + ROI(4), ROI(2) + ROI(4), ROI(2)], '-y', 'LineWidth', 3);
    
    % Ask the user to confirm if the ROI is good
    choice = questdlg('Is the ROI good?', 'User Feedback', 'Yes', 'No', 'Yes');
    
    switch choice
        case 'Yes'
            roijudge = 1;  % Exit the loop and proceed
        case 'No'
            close all;  % Close the figure and ask again
    end
end

% Save the selected ROI coordinates to a .mat file
save([beadFolder 'ROI.mat'], 'ROI');

% Save the maximum intensity projection image for reference
imwrite(ImMax, [beadFolder 'MaxofBeads.tif']);

%% Section2: Save the cropped stacks of the three fluorescence channels of beads using the ROI
% For each position, this section reads images from three channels (647, 561, 488),
% crops them according to the previously selected ROI, and appends the cropped images to multi-page TIFF stacks for each channel.
% Output:
%   - roi647stack.tif  : Cropped stack for 647 channel
%   - roi561stack.tif  : Cropped stack for 561 channel
%   - roi488stack.tif  : Cropped stack for 488 channel

for ii = 1 : Position
    % Define file paths for each channel's image at the current position
    foldname647 = [beadFolder 'Pos' num2str(ii) '\img_000000000_Channel2_000.tif'];
    foldname561 = [beadFolder 'Pos' num2str(ii) '\img_000000000_Channel3_000.tif'];
    foldname488 = [beadFolder 'Pos' num2str(ii) '\img_000000000_Channel4_000.tif'];
    
    % Read and crop the 647 channel image
    Im = imread(foldname647);
    Im647 = Im(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3));
    imwrite(Im647, [beadFolder 'roi647stack.tif'], 'WriteMode', 'append');
    
    % Read and crop the 561 channel image
    Im = imread(foldname561);
    Im561 = Im(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3));
    imwrite(Im561, [beadFolder 'roi561stack.tif'], 'WriteMode', 'append');
    
    % Read and crop the 488 channel image
    Im = imread(foldname488);
    Im488 = Im(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3));
    imwrite(Im488, [beadFolder 'roi488stack.tif'], 'WriteMode', 'append');
    
    % Display progress
    disp(['Processing position ' num2str(ii)]);
end

disp('All images were written!');

%% Section3: Transformation calculation for chromatic aberration correction
% This section reads the saved stacks for each fluorescence channel of the beads, calculates the maximum intensity projection (MIP) for each channel,
% and performs Chromatic aberration to align the 488 and 647 channels to the 561 channel.
% The resulting transformation matrices and aligned images are saved for later use.
% Outputs:
%   - roi488max.tif       : Max projection image of 488 channel
%   - roi561max.tif       : Max projection image of 561 channel (reference)
%   - roi647max.tif       : Max projection image of 647 channel
%   - roi488aligned.tif   : Aligned 488 channel max projection
%   - roi647aligned.tif   : Aligned 647 channel max projection
%   - tforms.mat          : Transformation objects for 488 and 647 channels

% Pre-allocate 3D arrays to store stacks
for ii = 1 : Position
    Im488temp(:,:,ii) = imread([beadFolder 'roi488stack.tif'], ii);
    Im561temp(:,:,ii) = imread([beadFolder 'roi561stack.tif'], ii);
    Im647temp(:,:,ii) = imread([beadFolder 'roi647stack.tif'], ii);
end

% Calculate maximum intensity projections (MIPs)
Im488Max = max(Im488temp, [], 3);
Im561Max = max(Im561temp, [], 3);
Im647Max = max(Im647temp, [], 3);

% Perform registration: align 488 and 647 channels to 561 channel
[alignedIm488, tform488] = imreg2Dr(Im561Max, Im488Max);
[alignedIm647, tform647] = imreg2Dr(Im561Max, Im647Max);

% Save the maximum projection images and aligned results
imwrite(Im488Max, [beadFolder 'roi488max.tif']);
imwrite(Im561Max, [beadFolder 'roi561max.tif']);
imwrite(Im647Max, [beadFolder 'roi647max.tif']);
imwrite(alignedIm488, [beadFolder 'roi488aligned.tif']);
imwrite(alignedIm647, [beadFolder 'roi647aligned.tif']);

% Save the transformation objects for future use
save([beadFolder 'tforms.mat'], 'tform488', 'tform647');

disp('Transformation matrices calculated and saved!');

%% Section4: Processing Drift-Corrected Average Intensity Projection Images for Chromatic Aberration Correction
% This section processes drift-corrected average intensity projection (AVG) images for Chromatic Aberration Correction.

[dirExpName] = uigetdir; % Select the root folder containing multiple experiments
[~, Cfolder, ~] = fileparts(dirExpName);

% Get all experiment subfolders
subfolders = dir(dirExpName);
subfolders = subfolders([subfolders.isdir]);
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));

% Loop over each experiment folder
for ii = 1:length(subfolders)
    subfolderPath = fullfile(dirExpName, subfolders(ii).name);
    
    % Get all fields of view (sub-subfolders)
    subsubfolders = dir(subfolderPath);
    subsubfolders = subsubfolders([subsubfolders.isdir]);
    subsubfolders = subsubfolders(~ismember({subsubfolders.name}, {'.', '..'}));
    
    for j = 1:length(subsubfolders)
        subsubfolderPath = fullfile(subfolderPath, subsubfolders(j).name);
        
        % Locate drift-corrected AVG images for each channel
        tiffFileBF = dir(fullfile(subsubfolderPath, '*0000_Channel1_000*.tif'));
        tiffFile647 = dir(fullfile(subsubfolderPath, '*AVG_channel2*.tif'));
        tiffFile561 = dir(fullfile(subsubfolderPath, '*AVG_channel3*.tif'));
        tiffFile488 = dir(fullfile(subsubfolderPath, '*AVG_channel4*.tif'));
        
        % Read images
        imgBF  = imread(fullfile(subsubfolderPath, tiffFileBF.name));
        img647 = imread(fullfile(subsubfolderPath, tiffFile647.name));
        img561 = imread(fullfile(subsubfolderPath, tiffFile561.name));
        img488 = imread(fullfile(subsubfolderPath, tiffFile488.name));
        
        % Crop images using predefined ROI
        imgBF  = imgBF(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3));
        img647 = img647(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3));
        img561 = img561(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3));
        img488 = img488(ROI(2):ROI(2)+ROI(4), ROI(1):ROI(1)+ROI(3));
        
        % Apply chromatic aberration correction transforms to 647 and 488 channels
        img647align = imwarp(img647, tform647, 'OutputView', imref2d(size(img647)));
        img488align = imwarp(img488, tform488, 'OutputView', imref2d(size(img488)));
        
        % Trim edges to remove potential artifacts from transformation and calibration
        imgBF        = imgBF(Edgesize:end-Edgesize, Edgesize:end-Edgesize);
        img647       = img647(Edgesize:end-Edgesize, Edgesize:end-Edgesize);
        img561       = img561(Edgesize:end-Edgesize, Edgesize:end-Edgesize);
        img488       = img488(Edgesize:end-Edgesize, Edgesize:end-Edgesize);
        img647align  = img647align(Edgesize:end-Edgesize, Edgesize:end-Edgesize);
        img488align  = img488align(Edgesize:end-Edgesize, Edgesize:end-Edgesize);
        
        % Construct output filenames for saving corrected and uncorrected images
        FilenameBF       = fullfile(dirExpName, [Cfolder '_BF_' subfolders(ii).name '.tif']);
        Filename647      = fullfile(dirExpName, [Cfolder '_647_' subfolders(ii).name '.tif']);
        Filename647align = fullfile(dirExpName, [Cfolder '_647align_' subfolders(ii).name '.tif']);
        Filename561      = fullfile(dirExpName, [Cfolder '_561_' subfolders(ii).name '.tif']);
        Filename488      = fullfile(dirExpName, [Cfolder '_488_' subfolders(ii).name '.tif']);
        Filename488align = fullfile(dirExpName, [Cfolder '_488align_' subfolders(ii).name '.tif']);
        
        % Save images as multi-frame TIFF stacks, overwrite for first frame, append thereafter
        if j == 1
            imwrite(imgBF, FilenameBF, 'WriteMode', 'overwrite');
            imwrite(img647, Filename647, 'WriteMode', 'overwrite');
            imwrite(img647align, Filename647align, 'WriteMode', 'overwrite');
            imwrite(img561, Filename561, 'WriteMode', 'overwrite');
            imwrite(img488, Filename488, 'WriteMode', 'overwrite');
            imwrite(img488align, Filename488align, 'WriteMode', 'overwrite');
        else
            imwrite(imgBF, FilenameBF, 'WriteMode', 'append');
            imwrite(img647, Filename647, 'WriteMode', 'append');
            imwrite(img647align, Filename647align, 'WriteMode', 'append');
            imwrite(img561, Filename561, 'WriteMode', 'append');
            imwrite(img488, Filename488, 'WriteMode', 'append');
            imwrite(img488align, Filename488align, 'WriteMode', 'append');
        end
    end
    disp(['Folder ' subfolders(ii).name ' processed']);
end
