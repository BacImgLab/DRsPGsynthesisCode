%% BeforeDemoProcess.m
%
% Script Overview:
% This script performs multiple processing steps on multi-channel TIFF images and corresponding .mat data files 
% obtained from experiments. It extracts and calculates intensity metrics, organizes data, merges left and right side data,
% renames fields for consistency, sorts data based on features, and finally extracts aligned and normalized channel intensity profiles 
% for subsequent analysis.
%
% Directory structure (assuming current working directory):
% - Original data files (.mat and corresponding .rot.tif) are located in the current directory (Processed folder)
% - Processed output files are saved into subfolders: extracted_files, extracted_files/S1, extracted_files/S1/SortedData, etc.
%
% Main Functional Sections:
%
% Section0: Extract intensity metrics from multi-channel TIFF files
%   - Iterate over all .mat files in the current folder, matching their corresponding .rot.tif files
%   - TIFF files contain multiple pages, each page corresponds to a channel (up to 4 channels)
%   - Compute non-zero pixel count (X), total intensity sum (Y), background-corrected intensity (E = Y - 100 * X), 
%     and average signal intensity per pixel (T = E / X) for each channel
%   - Write these metrics back into the 'data' structure inside the .mat files and save to the extracted_files folder
%
% Section1: Sort S1 data based on features
%   - Call the custom function sortMatFilesByFeatures to sort data by features (ds1, Ds1, dDratioS1) and save results
%
% Section2: Extract aligned normalized channel profiles for S0
%   - Load sorted S0 feature data and use extractAlignedChannels function to extract aligned and normalized intensity profiles
%
% Section3: Merge left (S1l) and right (S1r) S1 data
%   - Load .mat files from extracted_files, extract left and right side features (ds1, Ds1, dDratioS1), profiles and channel info
%   - Save these to extracted_files/S1 folder with filenames suffixed by 'L' or 'R' indicating side
%
% Section4: Rename left/right fields in S1 folder to unified field names
%   - Rename fields such as ds1l/ds1r to ds1, etc., in all .mat files in S1 folder for consistent processing
%
% Section5: Sort data in S1 folder by features
%   - Similar to Section1, sort S1 data by features and save results
%
% Section6: Extract aligned normalized channel profiles for S1
%   - Load sorted S1 feature data and extract aligned normalized intensity profiles using extractAlignedChannels
%
% Important Notes:
% - This script assumes a fixed background intensity of 100 units; modify the background subtraction formula in Section0 if different
% - TIFF filenames must strictly correspond to .mat files by replacing 'data.mat' with 'rot.tif'
% - Each section depends on the output of previous sections; run sequentially to ensure data integrity
% - Auxiliary functions 'sortMatFilesByFeatures' and 'extractAlignedChannels' must be defined and available in the path
% - Maintain folder structure and file naming conventions to avoid path errors
%
% Usage Recommendations:
% - Set the current directory to the folder containing the .mat and .tif files (Processed folder)
% - Run sections in order and monitor command window for errors
% - Verify field correspondences with data structure documentation
%
% Author:
% Date:
%
% -------------------- Begin Code --------------------
%% Section0 : Extract intensity metrics from multi-channel TIFF and update .mat files
clear;
clc;

% Get all .mat files in the current directory
ProcessedFiles = dir(fullfile(pwd, '*.mat'));

% Get all .tif files in the current directory
tifFiles = dir(fullfile(pwd, '*.tif'));

% Loop through all .mat files to find and process matching .tif files
for i = 1:length(ProcessedFiles)
    % Get current .mat file name and construct corresponding .tif filename
    ProcessedFileName = ProcessedFiles(i).name;
    tifFileName = strrep(ProcessedFileName, 'data.mat', 'rot.tif'); % Replace 'data.mat' with 'rot.tif'
    
    % Check if the corresponding .tif file exists
    if exist(fullfile(pwd, tifFileName), 'file') == 2
        % Read TIFF file info to determine number of channels (pages)
        tiffInfo = imfinfo(fullfile(pwd, tifFileName));
        numChannels = length(tiffInfo); % Number of pages/channels in TIFF stack

        % Initialize 3D array to store TIFF stack data (height x width x channels)
        tifStack = zeros(tiffInfo(1).Height, tiffInfo(1).Width, numChannels);

        % Read each channel/page from the TIFF file
        for ch = 1:numChannels
            tifStack(:,:,ch) = imread(fullfile(pwd, tifFileName), ch);
        end

        % Initialize intensity metric variables for up to 4 channels
        Channel1E = 0; Channel2E = 0; Channel3E = 0; Channel4E = 0;
        Channel1T = 0; Channel2T = 0; Channel3T = 0; Channel4T = 0;

        % Limit the number of channels to 4 maximum
        numChannels = min(numChannels, 4);

        % Calculate intensity metrics for each channel
        for ch = 1:numChannels
            % Extract the channel image data
            channelData = tifStack(:,:,ch);

            % Count non-zero pixels indicating signal presence (X)
            X = sum(channelData ~= 0, 'all');

            % Sum of all pixel intensities in the channel (Y)
            Y = sum(channelData, 'all');

            % Display intermediate metrics
            fprintf('Channel %d - X (Non-zero pixel count): %d\n', ch, X);
            fprintf('Channel %d - Y (Total pixel intensity): %f\n', ch, Y);

            % Calculate background-corrected total intensity (E) and average intensity (T)
            E = Y - 100 * X; % assuming background intensity of 100 per pixel
            T = E / X;

            % Display calculated intensity metrics
            fprintf('Channel %d - E (Background-corrected total intensity): %f\n', ch, E);
            fprintf('Channel %d - T (Average intensity per signal pixel): %f\n', ch, T);

            % Store metrics into corresponding variables based on channel number
            if ch == 1
                Channel1E = E;
                Channel1T = T;
            elseif ch == 2
                Channel2E = E;
                Channel2T = T;
            elseif ch == 3
                Channel3E = E;
                Channel3T = T;
            elseif ch == 4
                Channel4E = E;
                Channel4T = T;
            end
        end

        % Load the corresponding .mat file and extract 'data' structure
        matData = load(fullfile(pwd, ProcessedFileName));
        data = matData.data;

        % Update the 'data' structure with new intensity metrics
        data.Channel1E = Channel1E;
        data.Channel2E = Channel2E;
        data.Channel3E = Channel3E;
        data.Channel4E = Channel4E;
        data.Channel1T = Channel1T;
        data.Channel2T = Channel2T;
        data.Channel3T = Channel3T;
        data.Channel4T = Channel4T;

        % Create subfolder 'extracted_files' to save updated .mat files
        subfolderName = fullfile(pwd, 'extracted_files');
        if ~exist(subfolderName, 'dir')
            mkdir(subfolderName);
        end

        % Save the updated 'data' structure into the new subfolder
        save(fullfile(subfolderName, ProcessedFileName), 'data');

        disp(['Processed: ', ProcessedFileName]);
    else
        disp(['No matching .tif file found for: ', ProcessedFileName]);
    end
end
%% Section1 : Sort all .mat files in the subfolder by their feature values
featureList1 = {'ds1l', 'Ds1l', 'dDratioS1l', 'ds1r', 'Ds1r', 'dDratioS1r', 'ds0', 'Ds0', 'dDratioS0'};
saveNameList1 = {'ds1l', 'es1l', 'dDratioS1l', 'ds1r', 'es1r', 'dDratioS1r', 'ds0', 'es0', 'dDratioS0'};
sortedFolderPath = fullfile(subfolderName, 'Sorted');
sortMatFilesByFeatures(subfolderName, featureList1, saveNameList1, sortedFolderPath);

%% Section2: Extract Ds0-S0profile/ds0-S0profile/dDratioS0-S0profile as corresponding matrices for Ds0/ds0/dDratioS0 sorted Demo（S0）
sortedFolderPath = fullfile(subfolderName, 'sorted');
data = load(fullfile(sortedFolderPath, 'es0_Sorted.mat'));
cellData_Ds0 = data.sortedCells;
extractAlignedChannels(cellData_Ds0, 'S0profile', 60, 'e0_S0', sortedFolderPath);

sortedFolderPath = fullfile(subfolderName, 'sorted');
data = load(fullfile(sortedFolderPath, 'ds0_Sorted.mat'));
cellData_ds0 = data.sortedCells;
extractAlignedChannels(cellData_ds0, 'S0profile', 60, 'd0_S0', sortedFolderPath);

sortedFolderPath = fullfile(subfolderName, 'sorted');
data = load(fullfile(sortedFolderPath, 'dDratioS0_Sorted.mat'));
cellData_dDratioS0 = data.sortedCells;
extractAlignedChannels(cellData_dDratioS0, 'S0profile', 60, 'dDratioS0', sortedFolderPath);

%% Section3: Merge S1l and S1r into S1

% Define the paths
sourceFolder = fullfile('extracted_files');          % Folder containing raw .mat files
targetFolder = fullfile(sourceFolder, 'S1');          % Folder to save merged output files

% Create the target folder if it doesn't exist
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

% Retrieve all .mat files from the source folder
matFiles = dir(fullfile(sourceFolder, '*.mat'));

% Process each .mat file
for i = 1:length(matFiles)
    % Load the current file
    currentFile = matFiles(i).name;
    currentFilePath = fullfile(sourceFolder, currentFile);
    loadedData = load(currentFilePath);
    
    % For access convenience
    data = loadedData.data;  
    disp(['Processing file: ', currentFile]);

    % Extract and save S1r (right) data if fields exist
    if all(isfield(data, {'ds1r', 'Ds1r', 'dDratioS1r', 'S1rprofile', ...
                          'Channel1E', 'Channel2E', 'Channel3E', 'Channel4E', ...
                          'Channel1T', 'Channel2T', 'Channel3T', 'Channel4T'}))
                      
        % Assign variables
        ds1r = data.ds1r;
        Ds1r = data.Ds1r;
        dDratioS1r = data.dDratioS1r;
        S1rprofile = data.S1rprofile;
        Channel1E = data.Channel1E;
        Channel2E = data.Channel2E;
        Channel3E = data.Channel3E;
        Channel4E = data.Channel4E;
        Channel1T = data.Channel1T;
        Channel2T = data.Channel2T;
        Channel3T = data.Channel3T;
        Channel4T = data.Channel4T;
        
        % Generate filename and path
        [~, fileName, ~] = fileparts(currentFile);  % Remove extension
        newFileNameR = [fileName, 'R.mat'];         % Add suffix R
        newFilePathR = fullfile(targetFolder, newFileNameR);
        
        % Save S1r data
        save(newFilePathR, 'ds1r', 'Ds1r', 'dDratioS1r', 'S1rprofile', ...
             'Channel1E', 'Channel2E', 'Channel3E', 'Channel4E', ...
             'Channel1T', 'Channel2T', 'Channel3T', 'Channel4T');
        disp(['Saved: ', newFileNameR]);
    end

    % Extract and save S1l (left) data if fields exist
    if all(isfield(data, {'ds1l', 'Ds1l', 'dDratioS1l', 'S1lprofile', ...
                          'Channel1E', 'Channel2E', 'Channel3E', 'Channel4E', ...
                          'Channel1T', 'Channel2T', 'Channel3T', 'Channel4T'}))
                      
        % Assign variables
        ds1l = data.ds1l;
        Ds1l = data.Ds1l;
        dDratioS1l = data.dDratioS1l;
        S1lprofile = data.S1lprofile;
        Channel1E = data.Channel1E;
        Channel2E = data.Channel2E;
        Channel3E = data.Channel3E;
        Channel4E = data.Channel4E;
        Channel1T = data.Channel1T;
        Channel2T = data.Channel2T;
        Channel3T = data.Channel3T;
        Channel4T = data.Channel4T;
        
        % Generate filename and path
        newFileNameL = [fileName, 'L.mat'];         % Add suffix L
        newFilePathL = fullfile(targetFolder, newFileNameL);
        
        % Save S1l data
        save(newFilePathL, 'ds1l', 'Ds1l', 'dDratioS1l', 'S1lprofile', ...
             'Channel1E', 'Channel2E', 'Channel3E', 'Channel4E', ...
             'Channel1T', 'Channel2T', 'Channel3T', 'Channel4T');
        disp(['Saved: ', newFileNameL]);
    end
end

disp('All S1r and S1l components have been processed and saved to the S1 folder.');

%% Section4: Rename the fields in S1 folder files

% Define the path of .mat files in the S1 folder
s1MatFiles = dir(fullfile(targetFolder, '*.mat'));
% Step 1: Process all files with 'L' (left) in the filename
for i = 1:length(s1MatFiles)
    currentFile = s1MatFiles(i).name;
    currentFilePath = fullfile(targetFolder, currentFile);
    
    % Only process files with 'L' (case-insensitive)
    if contains(upper(currentFile), 'L')
        % Load current file
        loadedData = load(currentFilePath);
        disp(['Renaming fields in file (L): ', currentFile]);

        % Rename fields: ds1l → ds1, etc.
        if isfield(loadedData, 'ds1l')
            loadedData.ds1 = loadedData.ds1l;
            loadedData = rmfield(loadedData, 'ds1l');
        end
        if isfield(loadedData, 'Ds1l')
            loadedData.Ds1 = loadedData.Ds1l;
            loadedData = rmfield(loadedData, 'Ds1l');
        end
        if isfield(loadedData, 'dDratioS1l')
            loadedData.dDratioS1 = loadedData.dDratioS1l;
            loadedData = rmfield(loadedData, 'dDratioS1l');
        end
        if isfield(loadedData, 'S1lprofile')
            loadedData.S1profile = loadedData.S1lprofile;
            loadedData = rmfield(loadedData, 'S1lprofile');
        end

        % Save updated data
        save(currentFilePath, '-struct', 'loadedData');
        disp(['Saved renamed file: ', currentFile]);
    end
end
disp('✔ All L-files have completed field renaming.');
% Step 2: Process all files with 'R' (right) in the filename
for i = 1:length(s1MatFiles)
    currentFile = s1MatFiles(i).name;
    currentFilePath = fullfile(targetFolder, currentFile);
    
    % Only process files with 'R' (case-insensitive)
    if contains(upper(currentFile), 'R')
        % Load current file
        loadedData = load(currentFilePath);
        disp(['Renaming fields in file (R): ', currentFile]);

        % Rename fields: ds1r → ds1, etc.
        if isfield(loadedData, 'ds1r')
            loadedData.ds1 = loadedData.ds1r;
            loadedData = rmfield(loadedData, 'ds1r');
        end
        if isfield(loadedData, 'Ds1r')
            loadedData.Ds1 = loadedData.Ds1r;
            loadedData = rmfield(loadedData, 'Ds1r');
        end
        if isfield(loadedData, 'dDratioS1r')
            loadedData.dDratioS1 = loadedData.dDratioS1r;
            loadedData = rmfield(loadedData, 'dDratioS1r');
        end
        if isfield(loadedData, 'S1rprofile')
            loadedData.S1profile = loadedData.S1rprofile;
            loadedData = rmfield(loadedData, 'S1rprofile');
        end

        % Save updated data
        save(currentFilePath, '-struct', 'loadedData');
        disp(['Saved renamed file: ', currentFile]);
    end
end
disp('✔ All R-files have completed field renaming.');

%% Section5 : Sort all .mat files in the S1 folder by their feature values
featureList2 = {'ds1', 'Ds1', 'dDratioS1'};
saveNameList2 = {'ds1', 'es1', 'dDratioS1'};
targetFolder = fullfile('S1');
targetedFolder = fullfile(targetFolder, 'SortedData');

sortMatFilesByFeatures(targetFolder, featureList2, saveNameList2, targetedFolder);
disp('Check if cellData_ds1 exists:')
whos cellData_ds1

%% Section6: Extract ds1-S1profile/Ds1-S1profile/dDratioS1-S1profile as corresponding matrices for ds1/Ds1/dDratioS1 sorted Demo（S1）
data = load(fullfile(targetedFolder, 'ds1_Sorted.mat'));  % 加载文件
cellData_ds1 = data.sortedCells;  % 从结构体里提取变量
extractAlignedChannels(cellData_ds1, 'S1profile', 60, 'd1_S1', targetedFolder);

data = load(fullfile(targetedFolder, 'es1_Sorted.mat'));
cellData_Ds1 = data.sortedCells;
extractAlignedChannels(cellData_Ds1, 'S1profile', 60, 'e1_S1', targetedFolder);

data = load(fullfile(targetedFolder, 'dDratioS1_Sorted.mat'));
cellData_dDratioS1 = data.sortedCells;
extractAlignedChannels(cellData_dDratioS1, 'S1profile', 60, 'dDratioS1', targetedFolder);


%% Function Definitions: Extract and Align Normalized S0 and S1 Channel Profiles
function extractAlignedChannels(cellData, profileField, maxRows, outputPrefix, saveFolder)
% extractAlignedChannels
% General function to extract, normalize, and center-align profile data
% across 4 channels from a cell array of structs.
%
% INPUTS:
% - cellData:      cell array of structs, each with fields including profileField,
%                  and 'Channel1T', ..., 'Channel4T'
% - profileField:  string, name of the profile field to process ('S0profile' or 'S1profile')
% - maxRows:       fixed row size for padded vectors (e.g., 60)
% - outputPrefix:  prefix string for saving output .mat files (e.g., 'd0_S0' or 'd1_S1')
% - saveFolder:    full path where .mat files will be saved
%
% OUTPUTS:
% - Saves four .mat files into 'saveFolder' as:
%       [outputPrefix]Channel1.mat
%       [outputPrefix]Channel2.mat
%       [outputPrefix]Channel3.mat
%       [outputPrefix]Channel4.mat
%
% Each file contains a matrix where:
% - Rows correspond to aligned and normalized profile data (length = maxRows)
% - Columns correspond to individual cells

    numCells = length(cellData);
    channelDataCell = cell(4, numCells);

    for i = 1:numCells
        profileData = cellData{i}.(profileField);
        ChannelT = [cellData{i}.Channel1T, cellData{i}.Channel2T, ...
                    cellData{i}.Channel3T, cellData{i}.Channel4T];
        for col = 1:4
            normData = (profileData(:, col) - 100) / ChannelT(col);
            nonZeroData = normData(normData ~= 0);
            len = length(nonZeroData);
            if len < maxRows
                padTop = floor((maxRows - len) / 2);
                padBottom = maxRows - len - padTop;
                paddedData = [NaN(padTop, 1); nonZeroData; NaN(padBottom, 1)];
            else
                paddedData = nonZeroData(1:maxRows);
            end
            channelDataCell{col, i} = paddedData;
        end
    end

    for col = 1:4
        channelMatrix = cell2mat(channelDataCell(col, :)');
        save(fullfile(saveFolder, sprintf('%sChannel%d.mat', outputPrefix, col)), ...
             'channelMatrix');
    end

    fprintf('Channels saved: %sChannel1-4.mat in %s\n', outputPrefix, saveFolder);
end
function sortMatFilesByFeatures(sourceFolder, featureNames, saveNames, saveFolder)
% sortMatFilesByFeatures:
% Sort .mat files in sourceFolder by multiple features and save sorted cell arrays.
%
% Inputs:
% - sourceFolder: string, path of folder containing .mat files to process
% - featureNames: cell array of strings, feature field names to sort by (must exist in each file)
% - saveNames: cell array of strings, prefixes for saved sorted .mat files (one per feature)
% - saveFolder: string, path to folder to save sorted .mat files
%
% Outputs:
% - Saves one .mat file per feature in saveFolder, each containing a cell array "sortedCells"
%   with the data structs sorted by that feature in descending order.

    % Create save folder if it does not exist
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end

    % Get list of all .mat files in sourceFolder
    matFiles = dir(fullfile(sourceFolder, '*.mat'));
    nFiles = length(matFiles);
    nFeatures = length(featureNames);

    % Initialize structure array to hold file info and feature values
    fileData = struct('fileName', {}, 'data', {});
    for f = 1:nFeatures
        fileData(1).(featureNames{f}) = [];
    end

    % Load each file and extract features
    for i = 1:nFiles
        fullPath = fullfile(sourceFolder, matFiles(i).name);
        loaded = load(fullPath);

        % Support both data struct inside 'data' or fields at root
        if isfield(loaded, 'data')
            d = loaded.data;
        else
            d = loaded;
        end

        fileData(i).fileName = matFiles(i).name;
        fileData(i).data = d;

        % Extract feature values, error if missing
        for f = 1:nFeatures
            featName = featureNames{f};
            if isfield(d, featName)
                fileData(i).(featName) = d.(featName);
            else
                error('Feature field "%s" not found in file %s', featName, matFiles(i).name);
            end
        end
    end

    % Sort files by each feature in descending order and save results
    for f = 1:nFeatures
        featName = featureNames{f};
        saveName = saveNames{f};

        [~, sortIdx] = sort([fileData.(featName)], 'descend');
        sortedCells = cell(nFiles, 1);
        for i = 1:nFiles
            sortedCells{i} = fileData(sortIdx(i)).data;
        end

        save(fullfile(saveFolder, [saveName '_Sorted.mat']), 'sortedCells');
    end

    fprintf('Sorted data saved for features: %s\n', strjoin(featureNames, ', '));
end
