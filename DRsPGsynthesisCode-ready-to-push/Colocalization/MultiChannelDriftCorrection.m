%% Multichannel Image Preprocessing for Drift Correction
%
% This script performs preprocessing of multichannel fluorescence images,
% including BF (not used), and four fluorescence channels: Channel1 to Channel4.
% Each subfolder in the mainInputFolder represents a single field of view and contains multiple single-plane .tif images per channel.
%
% [Processing Steps]:
%   1. Initialize Fiji via the Miji interface;
%   2. Scan all subfolders under the specified main folder;
%   3. For each subfolder:
%       - Read the first 3 Z-plane images for each fluorescence channel;
%       - Create and save temporary Z-stacks for each channel;
%   4. Merge channels into a HyperStack;
%   5. Perform drift correction using the HyperStackReg plugin (Channel2 as reference);
%   6. Split and save the drift-corrected grayscale images;
%   7. Compute average intensity Z-projections for Channel2 to Channel4 and save them.
%
% [Requirements]:
%   - Fiji with Miji and HyperStackReg plugin installed;
%   - Correct Java path configuration for Miji in MATLAB;
%   - Image filenames must include channel keywords (e.g., "Channel1").
%
% [Output]:
%   - Channel Z-stacks: ChannelX_stack.tif;
%   - Drift-corrected images: ChannelX_000-0_2.tif;
%   - Z-projection images: AVG_ChannelX_000-0_2.tif;
%
% Cite:  
% MATLAB version: R2024a
% Recommended Fiji version: 1.54f or newer

clear; clc;

%% Section0: Initialization and Input Configuration
% This section sets up the Fiji-MATLAB connection using Miji, and specifies the input folder and image channels.

% Add Miji.jar to Java class path (provides MATLAB-Fiji interface)
javaaddpath 'C:\Fiji.app\scripts\mij.jar'; % Miji main interface JAR for controlling Fiji from MATLAB
% Add ImageJ core library (required for image operations)
javaaddpath 'C:\Fiji.app\jars\ij-1.54f.jar'; % Core ImageJ library, required for opening and saving images
% Add Miji scripts folder to MATLAB search path
addpath 'C:\Fiji.app\scripts'; % Ensures MATLAB can locate Miji-related functions

% (Optional) Add alternate Java paths for Miji/ImageJ (e.g., if installed elsewhere)
javaaddpath 'C:\MATLAB\SupportPackages\mij.jar';% Optional path to Miji JAR if installed with MATLAB support packages
javaaddpath 'C:\MATLAB\SupportPackages\ij-1.54f.jar';% Optional ImageJ JAR path

% Launch Fiji through Miji interface
Miji; 

% Set the main input directory containing subfolders (each subfolder = one field of view)
mainInputFolder = 'D:\ImageData\2025_ExampleProject';% Replace with the actual path to your dataset

% Search for all subfolders containing a specific keyword (e.g., "PBP1A")
subfolders = dir(mainInputFolder);                      
subfolders = subfolders([subfolders.isdir] & ...
    ~ismember({subfolders.name}, {'.', '..'}));        
subfolders = subfolders(contains({subfolders.name}, 'PBP1A')); 

% If no matching subfolders are found, stop execution
if isempty(subfolders)
    disp('No subfolders containing "PBP1A" were found. Script terminated.'); 
    return; % Exit the script
end

% Define the fluorescence channels to process (excluding BF channel)
channels = {'Channel1', 'Channel2', 'Channel3', 'Channel4'};  


%% Section1: Loop Through All Valid Subfolders
% This section processes each subfolder, performing Z-stack creation and drift correction for each fluorescence channel.

for sf = 1:length(subfolders)
    folderPath = fullfile(mainInputFolder, subfolders(sf).name); 
    fprintf('\nProcessing subfolder: %s\n', folderPath);          

    % Loop through all defined fluorescence channels (Channel1 to Channel4)
    for c = 1:length(channels)
        channelName = channels{c};% Current channel name (e.g., 'Channel2')

        % Search for all .tif images in this channel within the subfolder
        channelFiles = dir(fullfile(folderPath, ['*' channelName '*.tif']));  
        if isempty(channelFiles)
            fprintf('  No tif files found for %s, skipping...\n', channelName);  
            continue;
        end

        % Create Z-stack from the first 3 images of the current channel
        stackFile = fullfile(folderPath, [channelName '_stack.tif']); % Output path for the Z-stack
        createZStack(folderPath, stackFile, channelFiles); % Call function to generate and save Z-stack
    end

    % After all channels are processed, perform multichannel drift correction
    driftCorrectChannels(folderPath); % Call drift correction function

    % Close all image windows in Fiji to clean up memory
    ij.IJ.run('Close All'); 
end

% Display completion message once all subfolders have been processed
disp('All matching subfolders have been successfully processed!');


%% Function Definitions: createZStack and driftCorrectChannels
% These functions are used to preprocess multichannel Z-stack images:
%   - createZStack(): builds a Z-stack from individual image slices;
%   - driftCorrectChannels(): performs hyperstack merging, drift correction, and projection.

function createZStack(inputFolder, outputStackPath, channelFiles)
    % Build a Z-stack using the first up to 3 images in the given channel

    numFilesToUse = min(3, length(channelFiles));  % Use up to 3 images if available
    for i = 1:numFilesToUse
        filePath = fullfile(inputFolder, channelFiles(i).name);  % Full path to image
        ij.IJ.open(filePath);  % Open image in Fiji
    end

    ij.IJ.run('Images to Stack', 'name=Stack title=[] use'); % Convert opened images to Z-stack
    ij.IJ.saveAs('Tiff', outputStackPath);  % Save Z-stack to output path
    ij.IJ.run('Close All'); % Close all image windows in Fiji

    fprintf('Z-stack saved to %s\n', outputStackPath); 
end

function driftCorrectChannels(outputFolder)
    % Perform drift correction on multichannel Z-stacks using HyperStackReg

    % Define Z-stack file paths for each channel
    channel1Path = fullfile(outputFolder, 'Channel1_stack.tif');
    channel2Path = fullfile(outputFolder, 'Channel2_stack.tif');
    channel3Path = fullfile(outputFolder, 'Channel3_stack.tif');
    channel4Path = fullfile(outputFolder, 'Channel4_stack.tif');

    % Open all channel stacks in Fiji
    ij.IJ.open(channel1Path);
    ij.IJ.open(channel2Path);
    ij.IJ.open(channel3Path);
    ij.IJ.open(channel4Path);

    % Merge channels into a hyperstack (do not save this step)
    ij.IJ.run("Merge Channels...", ...
        "c1=Channel1_stack.tif c2=Channel2_stack.tif c3=Channel3_stack.tif c4=Channel4_stack.tif create");

    % Set hyperstack dimensions: 4 channels, 1 slice, 3 frames
    ij.IJ.run("Properties...", ...
        "channels=4 slices=1 frames=3 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");

    % Perform rigid-body drift correction using Channel2 as the reference
    ij.IJ.run("HyperStackReg ", ...
        "transformation=[Rigid Body] channel 2 show");

    % Split hyperstack into individual grayscale channel stacks
    ij.IJ.run("Split Channels");
    imageTitles = ij.WindowManager.getImageTitles();  % Get names of resulting image windows

    for i = 1:length(imageTitles)
        title = char(imageTitles(i));  % Convert Java string to MATLAB char

        % Determine save path for drift-corrected image based on title
        if contains(title, 'C1')
            savePath = fullfile(outputFolder, 'Channel1_000-0_2.tif');
        elseif contains(title, 'C2')
            savePath = fullfile(outputFolder, 'Channel2_000-0_2.tif');
        elseif contains(title, 'C3')
            savePath = fullfile(outputFolder, 'Channel3_000-0_2.tif');
        elseif contains(title, 'C4')
            savePath = fullfile(outputFolder, 'Channel4_000-0_2.tif');
        else
            continue;  % Skip if title doesn’t match any channel
        end

        ij.IJ.selectWindow(title);  % Bring the window to front
        ij.IJ.run("16-bit"); % Convert image to 16-bit
        ij.IJ.saveAs('Tiff', savePath);  % Save drift-corrected stack

        % Average Intensity Z Projection for Channels 2–4 ---
        if contains(title, 'C2') || contains(title, 'C3') || contains(title, 'C4')
            ij.IJ.run('Z Project...', 'projection=[Average Intensity]');  % Perform Z projection
            avgName = ['AVG_' strrep(title, '.tif', '') '.tif']; % Generate filename for projection
            avgPath = fullfile(outputFolder, avgName); % Output path
            ij.IJ.saveAs('Tiff', avgPath); % Save projection
            fprintf('  Saved Z-projection (Average Intensity): %s\n', avgPath);
        end
    end

    % Clean up: close all open image windows in Fiji
    ij.IJ.run("Close All");

    fprintf('Drift correction completed and saved to %s\n', outputFolder);
end
