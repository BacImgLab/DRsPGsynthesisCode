%% PreparePCCcaclu
%
% Description:
% This script processes multiple folders (stage2 to stage5), each of which
% contains a subfolder named "Processed". Within "Processed", it searches
% for all rotated image files matching the pattern 'Ch4_DR_*_rot.tif'.
% 
% For every such rotated image, the script extracts the image number,
% locates the corresponding original (non-rotated) image from the parent
% stage folder, and copies it to a new subfolder named "PCCcaclu"
% (created if it doesn't already exist).
%
% Purpose:
% Useful for preparing datasets (e.g., for PCC or colocalization analysis)
% where original unrotated images are needed alongside processed versions.
%
clear; clc;

% 1. Define the list of stage folders to process
stageFolders = {'stage2', 'stage3', 'stage4', 'stage5'};

% 2. Loop through each stage folder
for i = 1:length(stageFolders)
    stageFolder = stageFolders{i};  % Current stage folder name
    
    % Construct path to the 'Processed' subfolder
    processedFolder = fullfile(stageFolder, 'Processed');
    
    % Check if 'Processed' folder exists
    if ~isfolder(processedFolder)
        disp([' Processed folder not found in ', stageFolder]);
        continue;  % Skip to next stage if folder not found
    end
    
    % Search for all rotated images in 'Processed' folder
    rotFiles = dir(fullfile(processedFolder, 'Ch4_DR_*_rot.tif'));
    
    if isempty(rotFiles)
        disp(['️ No rotated images found in ', processedFolder]);
        continue;
    end
    
    % 3. Create PCCcaclu folder if it doesn't exist
    pccFolder = fullfile(stageFolder, 'PCCcaclu');
    if ~isfolder(pccFolder)
        mkdir(pccFolder);
        disp([' Created folder: ', pccFolder]);
    end
    
    % 4. Process each rotated image
    for j = 1:length(rotFiles)
        rotFileName = rotFiles(j).name;  % e.g., 'Ch4_DR_001_rot.tif'
        
        % Extract image number from file name
        prefix = 'Ch4_DR_';
        suffix = '_rot.tif';
        imageNumber = rotFileName(length(prefix)+1:end-length(suffix));
        
        % Construct path to the original image
        originalFile = fullfile(stageFolder, ['Ch4_DR_', imageNumber, '.tif']);
        
        % Check if the original file exists
        if isfile(originalFile)
            % Copy original image to PCCcaclu folder
            copyfile(originalFile, fullfile(pccFolder, ['Ch4_DR_', imageNumber, '.tif']));
            disp([' Copied: Ch4_DR_', imageNumber, '.tif → ', pccFolder]);
        else
            disp([' Original image not found: Ch4_DR_', imageNumber, '.tif in ', stageFolder]);
        end
    end
end

disp(' Done processing all stage folders.')

%% PCCcaclu
%
% Description:
% This script processes multi-channel TIF images across multiple stages
% (stage2, stage3, stage4, stage5), calculating the Pearson Correlation 
% Coefficients (PCC) between specific channels. It also calculates the 
% mean intensity for each channel, excluding zero-value (background) pixels.
%
% The script performs the following:
% 1. For each stage:
%    - Reads each 4-layer stacked TIF image (Ch4_DR_*.tif).
%    - Calculates the PCC values between channels:
%         - Channel 2 vs Channel 3 (C2C3)
%         - Channel 2 vs Channel 4 (C2C4)
%         - Channel 3 vs Channel 4 (C3C4)
%    - Computes the mean intensity of each channel (after removing background and subtracting 100).
%    - Saves the results as a CSV file (e.g., stage2_PCC_Values.csv).
% 2. Merges all per-stage CSV files into a single consolidated file named 'Merged_PCC_Values.csv'.
%

clear;
clc;

%% ---------------------- Part 1: Calculate PCC for Each Stage ----------------------

% 1. Define working directory and stage folders
parentDir = pwd;  % Current working directory
stages = {'stage2', 'stage3', 'stage4', 'stage5'};  % Folders representing each stage

% 2. Loop through each stage
for stageIdx = 1:length(stages)
    stageName = stages{stageIdx};
    stageFolderPath = fullfile(parentDir, stageName);  % Path to current stage folder
    
    % Get all 4-channel TIF image files (assumed format: Ch4_DR_*.tif)
    tifFiles = dir(fullfile(stageFolderPath, 'Ch4_DR_*.tif'));
    
    % Preallocate results: columns for filename, PCCs, and mean intensities
    pccResults = cell(length(tifFiles), 7);  
    
    % 3. Process each TIF stack
    for i = 1:length(tifFiles)
        tifFilePath = fullfile(stageFolderPath, tifFiles(i).name);
        
        % Load the first layer of the 4-channel image stack
        stack = imread(tifFilePath, 'Index', 1);
        
        % Read remaining layers (channels 2 to 4)
        for j = 2:4
            stack(:,:,j) = imread(tifFilePath, 'Index', j);
        end
        
        % Extract individual channels and convert to double
        A2 = double(stack(:,:,2));  % Channel 2
        A3 = double(stack(:,:,3));  % Channel 3
        A4 = double(stack(:,:,4));  % Channel 4
        
        % 4. Create a mask to exclude background (zero-value pixels)
        mask = (A2 > 0) & (A3 > 0) & (A4 > 0);
        
        % 5. Calculate Pearson Correlation Coefficients (PCC)
        PCC_C2C3 = corr(A2(mask), A3(mask));  % Between Channel 2 and 3
        PCC_C2C4 = corr(A2(mask), A4(mask));  % Between Channel 2 and 4
        PCC_C3C4 = corr(A3(mask), A4(mask));  % Between Channel 3 and 4
        
        % 6. Compute mean intensity of each channel (exclude 0s and subtract 100)
        meanC2 = mean(A2(A2 > 0), 'all') - 100;
        meanC3 = mean(A3(A3 > 0), 'all') - 100;
        meanC4 = mean(A4(A4 > 0), 'all') - 100;
        
        % 7. Store results
        pccResults{i, 1} = tifFiles(i).name;
        pccResults{i, 2} = PCC_C2C3;
        pccResults{i, 3} = PCC_C2C4;
        pccResults{i, 4} = PCC_C3C4;
        pccResults{i, 5} = meanC2;
        pccResults{i, 6} = meanC3;
        pccResults{i, 7} = meanC4;
    end
    
    % 8. Save per-stage PCC results to a CSV file
    pccTable = cell2table(pccResults, ...
        'VariableNames', {'FileName', 'C2C3', 'C2C4', 'C3C4', 'meanC2', 'meanC3', 'meanC4'});
    
    csvFileName = sprintf('%s_PCC_Values.csv', stageName);
    csvFilePath = fullfile(stageFolderPath, csvFileName);
    writetable(pccTable, csvFilePath);
    
    fprintf('PCC values for %s saved to %s\n', stageName, csvFilePath);
end

%% ---------------------- Part 2: Merge All CSVs into One File ----------------------

% Initialize merged table
mergedData = [];

% Loop through each stage again to collect results
for stageIdx = 1:length(stages)
    stageName = stages{stageIdx};
    stageFolderPath = fullfile(parentDir, stageName);
    
    % CSV file name for current stage
    csvFileName = sprintf('%s_PCC_Values.csv', stageName);
    csvFilePath = fullfile(stageFolderPath, csvFileName);
    
    % Check if CSV file exists
    if isfile(csvFilePath)
        % Read the CSV data
        currentData = readtable(csvFilePath);
        
        % Append to merged data
        if isempty(mergedData)
            mergedData = currentData;
        else
            mergedData = [mergedData; currentData];
        end
    else
        warning('File not found: %s\n', csvFilePath);
    end
end

% Save merged results
outputFileName = 'Merged_PCC_Values.csv';
outputFilePath = fullfile(parentDir, outputFileName);
writetable(mergedData, outputFilePath);

disp(['Merged PCC CSV saved as: ', outputFilePath]);
