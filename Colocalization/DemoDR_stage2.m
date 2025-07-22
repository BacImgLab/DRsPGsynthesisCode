%% Construct the septal intensity profiles of different channels for D. radiodurans stage 2 cells
% This script processes fluorescence microscopy images of D. radiodurans cells at division stage 2.
% It performs the following steps:
% 1. User selects a folder containing segmented multi-channel TIFF images.
% 2. For each image:
%    - Reads the brightfield channel to compute the rotation angle and bounding box of the cell.
%    - Rotates all fluorescence channels to horizontally align the cell.
%    - Displays the rotated multi-channel image stack for user quality control and acceptance.
%    - If accepted, uses the membrane channel to interactively select the septum S0 line position.
%    - Extracts the intensity profile along septum S0 for all fluorescence channels,
%      using a specified linewidth for averaging.
%    - Calculates septal geometric parameters including:
%        * D (cell width from bounding box long axis),
%        * d (septum length from intensity profile),
%        * d/D ratio indicating septal constriction.
%    - Saves rotated images, septum position data, profiles, and parameters in a 'Processed' subfolder.
%
% Parameters:
%   Wline        - Width in pixels of the line for averaging intensity profiles.
%   ChannelNum   - Number of image channels in the TIFF files.
%   Channelbf    - Index of the brightfield channel for rotation angle calculation.
%   ChannelMem   - Index of the membrane channel for septum selection and profile calculation.
%
% Outputs:
%   - Rotated multi-channel TIFF images saved with suffix '_rot.tif'.
%   - Data files (.mat) containing septal positions, rotation info, geometric parameters, and profiles.
%   - CSV files for septal intensity profiles (S0 profile saved; S1L and S1R placeholders present).
%
% Dependencies:
%   - bf2rotateRec: calculates rotation angle and bounding box from brightfield image.
%   - DR4cShow: interactive quality control GUI for selecting usable cells.
%   - S0select: interactive tool to select septum S0 position on membrane channel image.
%   - lineProfile: extracts intensity profile along a specified line with averaging width. 
%
% This script supports quantitative analysis of septal morphology and fluorescence intensity 
% distributions in dividing D. radiodurans cells at early division stages.

clear; clc;

% 1. setup the parameters
Wline = 7; % The linewidth used to construct the intensity profile (in pixels)
ChannelNum = 4; % Total number of image channels

Channelbf = 1; % Index of brightfield channel used to determine rotation angle
ChannelMem = 2; % Index of membrane channel used to find septum position

% 2. select folder(or folders) for processing
selpath = uigetdir('', 'Select the folder contains segmented cells'); % user selects folder with images
resultFolder = fullfile(selpath, 'Processed'); % folder to save processed results
if ~exist(resultFolder, 'dir')
    mkdir(resultFolder);  % Create the folder if it doesn't exist
end

% Find all TIFF image files in the selected folder
tifFiles = dir(fullfile(selpath, '*.tif'));

% 3. walk through all the images in the folder
for idf = 1 : numel(tifFiles)
    imageFilename = tifFiles(idf).name;              % get the image filename
    [~, name, ext] = fileparts(imageFilename);        % get the filename without extension
    imageFilepath = fullfile(selpath, imageFilename); % full path for reading the original image

    % Define output filenames for saving processed data and profiles
    newImgName = strcat(name, '_rot', '.tif');  
    newDataName = strcat(name, '_data', '.mat');
    newprofName1 = strcat(name, '_S0', '.csv');
    newprofName2 = strcat(name, '_S1L', '.csv'); % unused here, but kept for consistency
    newprofName3 = strcat(name, '_S1R', '.csv'); % unused here, but kept for consistency

    imageSavename = fullfile(resultFolder, newImgName);   % save path for rotated images
    dataSavename = fullfile(resultFolder, newDataName);   % save path for septum positions and metadata
    profileSavename1 = fullfile(resultFolder, newprofName1); % save path for S0 intensity profile
    profileSavename2 = fullfile(resultFolder, newprofName2); % placeholder for S1L profile path
    profileSavename3 = fullfile(resultFolder, newprofName3); % placeholder for S1R profile path

    % Read the brightfield image to determine cell rotation angle and bounding box
    imgbf = imread(imageFilepath, Channelbf);
    [angle, rotatedBox] = bf2rotateRec(imgbf);

    % 4. read all the channels and rotate images to align cells horizontally
    for ii = 1 : ChannelNum
        imgTemp = imread(imageFilepath, ii);           % read each channel
        imgRot = imrotate(imgTemp, -angle, 'bilinear', 'crop'); % rotate image by negative angle to horizontal
        ImageAllChannel(:,:,ii) = imgRot;              % store rotated image in 3D array
    end

    % 5. determine whether the current cell image is usable or not
    % By default imgFlag=1 means good image, but here we use an interactive GUI to check quality
    imgFlag = DR4cShow(ImageAllChannel); % display images and get user input for acceptance

    % 6. if the cell is approved, find septum position and extract profiles
    if imgFlag == 1
        % Use the membrane channel image to select S0 septum position interactively
        ImShowing = ImageAllChannel(:,:,ChannelMem);
        [xS0, yS0] = S0select(ImShowing, rotatedBox); % user selects two points defining septum S0

        % Initialize matrix to store S0 intensity profiles from all channels
        S0profile = [];

        % Extract intensity profile along septum S0 for each channel
        for ii = 1 : ChannelNum
            imgC = ImageAllChannel(:,:,ii); % current channel image
            % Extract intensity profile along line defined by (xS0,yS0), with width Wline
            [S0p, lineX0, lineY0, Len0] = lineProfile(imgC, xS0, yS0, Wline);
            S0profile = [S0profile, S0p']; % append profile as a new column

            % For membrane channel, calculate septum geometry d, D and ratio d/D
            if ii == ChannelMem
                Ds0 = rotatedBox(3); % Cell width (long axis length) from bounding box
                ds0 = Len0;          % Length of the septum profile
                DratioS0 = ds0 / Ds0; % Ratio of septum length to cell width
            end
        end

        % 7. save all the results into a .mat file for further analysis
        data.LACell = rotatedBox(3); % Long axis length of the cell
        data.SACell = rotatedBox(4); % Short axis length of the cell
        data.LoS0 = sqrt(diff(xS0)^2 + diff(yS0)^2); % Euclidean length of septum S0

        % Save septum coordinates and full line coordinates for profile plotting
        data.xS0 = xS0; 
        data.yS0 = yS0;
        data.XYS0 = [lineX0', lineY0'];

        % Save rotated bounding box and rotation angle
        data.rotatedBox = rotatedBox;
        data.angle = angle;

        % Save septum geometric parameters
        data.ds0 = ds0;       % septum length
        data.Ds0 = Ds0;       % cell width
        data.dDratioS0 = DratioS0; % ratio

        % Save the intensity profiles of septum S0 across all channels
        data.S0profile = S0profile;

        % Save the data struct to the .mat file
        save(dataSavename, 'data');

        % 8. write the rotated images for all channels into a multi-page TIFF file
        imwrite(ImageAllChannel(:,:,1), imageSavename, "tif", "WriteMode", "overwrite", "Compression", "none");
        for jj = 2 : ChannelNum
            imwrite(ImageAllChannel(:,:,jj), imageSavename, "tif", "WriteMode", "append", "Compression", "none");
        end
    end
end
