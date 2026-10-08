% This script merges four single-channel TIFF images (C1, C2, C3, C4),
% assumed to be saved under separate folders, into a single multi-frame TIFF file per field. 
% The input folders should follow the structure:
%   basePath\C1\stageName
%   basePath\C2\stageName
%   basePath\C3\stageName
%   basePath\C4\stageName
% The merged files will be saved under:
%   basePath\stageName\Ch4_xxx.tif
clear;clc;
%% Section0: Define the base directory and stage folder name
basePath = 'D:\ImageData\SA\satest\FtsW-PBP2\mem\123';  % Replace with your actual root folder path
stageName = 'stage1';  % Replace with your actual stage folder name
outputFolder = fullfile(basePath, stageName);  
% Define the folder format to read images from each channel's subfolder
readStr = ['%s\\' stageName];  % Used to build full paths for each channel (e.g., C1\stage1)
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
%% Section 1: Collect .tif filenames from all four channels (C1–C4)
% This section scans subdirectories C1, C2, C3, and C4 under the base path,
% specifically in each one's subfolder named after the stage (e.g., 'stage1').
% It extracts and stores the list of .tif files for each channel in a cell array.

channelFolders = {'C1', 'C2', 'C3', 'C4'};  % Names of subfolders for the four image channels
fileList = cell(1, 4);  % Preallocate a cell array to store file name lists for each channel

for ch = 1:4
    % Construct the full path to the subfolder of the current channel
    channelFolderPath = fullfile(basePath, sprintf(readStr, channelFolders{ch}));
    
    % Retrieve a list of all .tif files in the current channel folder
    tifFiles = dir(fullfile(channelFolderPath, '*.tif'));
    
    % Store the file names (without paths) in the corresponding cell of fileList
    fileList{ch} = {tifFiles.name};
end
%% Section 2: Combine corresponding images from all four channels into a single multi-channel TIFF
% This section loops through all image files in the first channel folder (C1),
% assumes the other channels have matching filenames, reads corresponding images
% from each channel, merges them into a 3D array, and saves them as a multi-page TIFF.

for i = 1:length(fileList{1})
    % Get the current image file name from channel 1 (used as reference)
    imgName = fileList{1}{i};

    % Preallocate a cell array to hold the 2D image data from each channel
    channels = cell(1, 4);

    % Loop through each of the four channels (C1 to C4)
    for ch = 1:4
        % Construct the full file path for the image in the current channel folder
        channelFolderPath = fullfile(basePath, sprintf(readStr, channelFolders{ch}));
        imgPath = fullfile(channelFolderPath, imgName);

        % Check if the image exists at the constructed path
        if exist(imgPath, 'file')
            % Read the image and store it in the corresponding channel slot
            channels{ch} = imread(imgPath);
        else
            % If the file is missing, show a warning and leave the channel empty
            warning('Image %s not found in channel %d folder: %s', imgName, ch, imgPath);
        end
    end

    % Combine the four channel images into a single 3D matrix: height × width × 4
    multiChannelImage = cat(3, channels{1}, channels{2}, channels{3}, channels{4});

    % Extract the base file name (without extension) from the original file
    [~, name, ~] = fileparts(imgName);

    % Build the full path for the output file in the target folder
    outputFileName = fullfile(outputFolder, sprintf('Ch4_%s.tif', name));

    % Save the first channel as the base page of a multi-page TIFF
    imwrite(multiChannelImage(:,:,1), outputFileName, ...
        'WriteMode', 'overwrite', 'Compression', 'none');

    % Append the remaining three channels as additional pages in the TIFF
    for ii = 2:4
        imwrite(multiChannelImage(:,:,ii), outputFileName, ...
            'WriteMode', 'append', 'Compression', 'none');
    end

    % Display current progress
    display(['Combining file ' num2str(i)]);
end
