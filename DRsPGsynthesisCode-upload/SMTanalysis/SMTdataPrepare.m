%% SMTdataPrepare
%
% Description:
% This script organizes single-molecule tracking (SMT) image data from dual-view microscopy(split tracking). 
% It handles three types of 16-bit TIFF images:
%   - Bright field (B): contains label 'B' in the filename
%   - Fluorescence (F): contains label 'F' (e.g., Fluo_001.tif)
%   - SMT/Tracking (S): contains label 'S' (e.g., STORM_Cycle0_001.tif)
%
% It crops predefined ROIs, stacks or averages the images, and saves
% organized output into a new folder for downstream analysis.
%
clear;              
clc;               

% Define ROI for long-wavelength channel (usually SMT or BF)
ROI1 = [32, 184, 743, 1022];     % Format: [x, y, width, height]
% Define ROI for short-wavelength channel (usually Fluorescence)
ROI2 = [828, 177, 743, 1022];    % Format: [x, y, width, height]

n_BF = 1;            % Number of bright field images to include in stack
n_FL = 50;           % Number of fluorescence slices to average

% Let user select the parent directory containing subfolders of images
rootDir = uigetdir([], 'Select the Parent Directory Containing Image Folders');  

% Create output folder name using current date
outputRootDir = fullfile(rootDir, ['OrganizedIm-', datestr(now, 'mm-dd-yyyy')]);

% Check if output folder exists; if not, create it
if ~exist(outputRootDir, 'dir')
    mkdir(outputRootDir);       % Create directory if it doesn't exist
    fprintf('Output folder is created: %s\n', outputRootDir);
else
    fprintf('Output folder already exists: %s\n', outputRootDir);
end

% List all entries in root directory
subfolders = dir(rootDir);

% Keep only directories
subfolders = subfolders([subfolders.isdir]);

% Remove system entries and output folder from list
subfolders = subfolders(~ismember({subfolders.name}, ...
    {['OrganizedIm-', datestr(now, 'mm-dd-yyyy')], '.', '..'}));

% Count number of valid subfolders
numSubfolders = length(subfolders);

% Display how many folders will be processed
fprintf('Starting on %d folders ...\n', numSubfolders);

for k = 1:numSubfolders
    subfolderName = subfolders(k).name;              % Get subfolder name
    subfolderPath = fullfile(rootDir, subfolderName);% Full path to subfolder

    % Get all .tif image files in the subfolder
    files = dir(fullfile(subfolderPath, '*.tif'));

    % Initialize image lists for each type
    images_B = {};    % Bright field
    images_F = {};    % Fluorescence
    images_S = {};    % SMT / Tracking

    % Classify each image file based on filename pattern
    for i = 1:length(files)
        fileName = files(i).name;                      % Get filename
        fullPath = fullfile(subfolderPath, fileName);  % Full path to image
        if contains(fileName, 'B')                     % Bright field
            images_B{end+1,1} = fullPath;
        elseif contains(fileName, 'F')                 % Fluorescence
            images_F{end+1,1} = fullPath;
        elseif contains(fileName, 'S')                 % SMT / Tracking
            images_S{end+1,1} = fullPath;
        end
    end

    % Sort fluorescence images based on number in filename (Fluo_XXX)
    if ~isempty(images_F)
        numbers = regexp(images_F, '(?<=Fluo_)\d+', 'match');     % Extract numbers
        numbers = cellfun(@(x) str2double(x{1}), numbers);        % Convert to numeric
        [~, sortIdx] = sort(numbers);                             % Get sort indices
        images_F = images_F(sortIdx);                             % Reorder filenames
    end

    % Sort SMT images based on number in filename (STORM_Cycle0_XXX)
    if ~isempty(images_S)
        numbers = regexp(images_S, '(?<=STORM_Cycle0_)\d+', 'match'); % Extract numbers
        numbers = cellfun(@(x) str2double(x{1}), numbers);            % Convert to numeric
        [~, sortIdx] = sort(numbers);                                 % Sort
        images_S = images_S(sortIdx);                                 % Reorder
    end

    % Display current subfolder being processed
    fprintf('Processing %d: %s\n', k, subfolderName);

    % Stack and crop bright field images
    if ~isempty(images_B)
        BF_stackALL = stackImages(images_B, ROI1, n_BF);
    else
        BF_stackALL = [];
    end

    % Stack and average fluorescence images
    if ~isempty(images_F)
        FL_stackALL = stackImages(images_F, ROI2, n_FL);
    else
        FL_stackALL = [];
    end

    % Stack SMT tracking images
    if ~isempty(images_S)
        trackALL = stackImages(images_S, ROI1);
    else
        trackALL = [];
    end

    % Create output subfolder for current dataset
    outputSubfolder = fullfile(outputRootDir, subfolderName);
    if ~exist(outputSubfolder, 'dir')
        mkdir(outputSubfolder);  % Make subfolder if not present
        fprintf('Subfolder created: %s\n', outputSubfolder);
    end

    % Save bright field stack and first frame
    if ~isempty(BF_stackALL)
        saveStack(BF_stackALL, fullfile(outputSubfolder, ['BF' num2str(k-1) '_stack.tif']));  % Full stack
        imwrite(BF_stackALL(:, :, 1), fullfile(outputSubfolder, ['BF' num2str(k-1) '.tif']), 'WriteMode', 'overwrite');  % First frame
    end

    % Save fluorescence stack and averaged projection
    if ~isempty(FL_stackALL)
        saveStack(FL_stackALL, fullfile(outputSubfolder, ['FL' num2str(k-1) '_stack.tif']));  % Full stack

        % Perform average projection (Z-mean)
        FL_allread = im2double(FL_stackALL);        % Convert to double precision
        FL_meanreal = mean(FL_allread, 3);          % Compute mean across slices
        FL_mean16 = uint16(FL_meanreal * 65535);    % Convert to 16-bit scale
        resultImagePath = fullfile(outputSubfolder, ['FL' num2str(k-1) '.tif']); % Output file path
        imwrite(FL_mean16, resultImagePath, 'WriteMode', 'overwrite'); % Save average projection
        fprintf('FL image saved in: %s\n', resultImagePath);           % Display save path
    end

    % Save SMT tracking image stack
    if ~isempty(trackALL)
        saveStack(trackALL, fullfile(outputSubfolder, ['track' num2str(k-1) '.tif']));
    end

    % Print completion message for current folder
    fprintf('Subfolder %d: %s finished\n', k, subfolderName);
end

% Print final completion message
fprintf('All files finished!\n');

%% ------------------ Function: stackImages -------------------------------
% This function reads image files, crops to ROI, and stacks them into a 3D array
function stackedImages = stackImages(imagePaths, roi, nslice)
    if nargin < 3
        nslice = length(imagePaths);  % Use all if nslice not specified
    end

    % Get ROI coordinates
    x = roi(1);
    y = roi(2);
    width = roi(3);
    height = roi(4);

    images = cell(1, nslice);   % Preallocate image list

    for i = 1:nslice
        img = imread(imagePaths{i});                                  % Read image
        images{i} = img(y:y+height-1, x:x+width-1);                   % Crop ROI
    end

    stackedImages = cat(3, images{:});    % Stack into 3D array
end

%% ------------------ Function: saveStack -------------------------------
% This function saves a 3D image stack to a multi-page TIFF file
function saveStack(imageStack, outputPath)
    for i = 1:size(imageStack, 3)
        if i == 1
            imwrite(imageStack(:, :, i), outputPath, 'WriteMode', 'overwrite');  % Save first slice
        else
            imwrite(imageStack(:, :, i), outputPath, 'WriteMode', 'append');     % Append additional slices
        end
    end
end
