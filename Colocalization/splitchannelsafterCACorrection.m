%% splitmultichannelsafterCACorrection.m
% This script splits multi-frame multi-channel TIFF images into individual
% single-frame TIFF files, saving them into separate folders for each channel, after chromatic aberration correction.
% Folder structure assumption:
% rootFolder/
%   subfolder1/
%     subsubfolder1/
%       multi-frame TIFF files for BF, 647, 561, 488 channels
% Output structure:
% rootFolder/subfolderX/subsubfolderY/position/C[channel]/frame_channel.tif
% Channels:
%   C1 - BF (brightfield)
%   C2 - 647 nm fluorescence
%   C3 - 561 nm fluorescence
%   C4 - 488 nm fluorescence
%Cite
clear; clc; 

% Select the root folder containing all sample subfolders
rootFolder = uigetdir;
if rootFolder == 0
    error('No folder selected. Exiting.'); 
end

% Get all subfolders (experiments) inside the root folder
subfolders = dir(rootFolder);
subfolders = subfolders([subfolders.isdir]); 
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); 

% Loop through each experiment subfolder
for ii = 1:length(subfolders)
    subfolderPath = fullfile(rootFolder, subfolders(ii).name); 
    
    % Get all position sub-subfolders inside this experiment folder
    subsubfolders = dir(subfolderPath);
    subsubfolders = subsubfolders([subsubfolders.isdir]);
    subsubfolders = subsubfolders(~ismember({subsubfolders.name}, {'.', '..'}));
    
    % Loop through each position folder
    for jj = 1:length(subsubfolders)
        subsubfolderPath = fullfile(subfolderPath, subsubfolders(jj).name); 
        
        % Find TIFF files for each channel based on filename patterns
        tiffFilesBF  = dir(fullfile(subsubfolderPath, '*BF*.tif'));  % Brightfield channel
        tiffFiles647 = dir(fullfile(subsubfolderPath, '*647*.tif')); % 647 nm fluorescence channel
        tiffFiles561 = dir(fullfile(subsubfolderPath, '*561*.tif')); % 561 nm fluorescence channel
        tiffFiles488 = dir(fullfile(subsubfolderPath, '*488*.tif')); % 488 nm fluorescence channel
        
        % --- Process BF channel TIFF files ---
        for kk = 1:length(tiffFilesBF)
            filename = tiffFilesBF(kk).name;
            % Check if filename contains the pattern '_kk' indicating the correct file
            if contains(filename, ['_' num2str(kk)])
                filePath = fullfile(subsubfolderPath, filename); 
                info = imfinfo(filePath); 
                numFrames = length(info); 
                
                % Create output folder: position/kk/C1 for BF channel
                outFolder = fullfile(subsubfolderPath, num2str(kk), 'C1');
                if ~isfolder(outFolder) 
                    mkdir(outFolder);
                end
                
                % Loop through each frame and save as individual TIFF
                for frameIdx = 1:numFrames
                    img = uint16(imread(filePath, frameIdx)); 
                    outName = fullfile(outFolder, [num2str(frameIdx) 'C1.tif']); 
                    imwrite(img, outName); 
                end
            end
        end
        
        % --- Process 488 nm channel TIFF files (similar process) ---
        for kk = 1:length(tiffFiles488)
            filename = tiffFiles488(kk).name;
            if contains(filename, ['_' num2str(kk)])
                filePath = fullfile(subsubfolderPath, filename);
                info = imfinfo(filePath);
                numFrames = length(info);
                
                outFolder = fullfile(subsubfolderPath, num2str(kk), 'C4'); 
                if ~isfolder(outFolder)
                    mkdir(outFolder);
                end
                
                for frameIdx = 1:numFrames
                    img = uint16(imread(filePath, frameIdx));
                    outName = fullfile(outFolder, [num2str(frameIdx) 'C4.tif']);
                    imwrite(img, outName);
                end
            end
        end
        
        % --- Process 561 nm channel TIFF files ---
        for kk = 1:length(tiffFiles561)
            filename = tiffFiles561(kk).name;
            if contains(filename, ['_' num2str(kk)])
                filePath = fullfile(subsubfolderPath, filename);
                info = imfinfo(filePath);
                numFrames = length(info);
                
                outFolder = fullfile(subsubfolderPath, num2str(kk), 'C3'); 
                if ~isfolder(outFolder)
                    mkdir(outFolder);
                end
                
                for frameIdx = 1:numFrames
                    img = uint16(imread(filePath, frameIdx));
                    outName = fullfile(outFolder, [num2str(frameIdx) 'C3.tif']);
                    imwrite(img, outName);
                end
            end
        end
        
        % --- Process 647 nm channel TIFF files ---
        for kk = 1:length(tiffFiles647)
            filename = tiffFiles647(kk).name;
            if contains(filename, ['_' num2str(kk)])
                filePath = fullfile(subsubfolderPath, filename);
                info = imfinfo(filePath);
                numFrames = length(info);
                
                outFolder = fullfile(subsubfolderPath, num2str(kk), 'C2'); 
                if ~isfolder(outFolder)
                    mkdir(outFolder);
                end
                
                for frameIdx = 1:numFrames
                    img = uint16(imread(filePath, frameIdx));
                    outName = fullfile(outFolder, [num2str(frameIdx) 'C2.tif']);
                    imwrite(img, outName);
                end
            end
        end
        
        % Display progress info for current subfolder and subsubfolder
        fprintf('Processed subfolder %s, subsubfolder %s\n', subfolders(ii).name, subsubfolders(jj).name);
    end
end
