%% Demograph Processing Script
% ---------------------------------------------------------------
% This script performs the following operations for multiple channels:
% 1. Loads the demograph matrix from a .mat file.
% 2. Smooths the matrix along the cell cycle by averaging columns.
% 3. Normalizes the resulting matrix by its maximum value.
% 4. Saves the normalized matrix as a .mat file and as a grayscale .tif image.
%
% The script processes three channels: Channel 2, Channel 3, and Channel 4.
% The input .mat files must contain variables: 
%   - dsDChannel2new
%   - dsDChannel3new
%   - dsDChannel4new
%
% Output:
%   - demo_sm_norm2new.mat / demo2.tif
%   - demo_sm_norm3new.mat / demo3.tif
%   - demo_sm_norm4new.mat / demo4.tif
%
% Dependencies:
%   - Requires the function `average_colN` to be available in the path.
% ---------------------------------------------------------------

%% Channel 2 Processing
N = 80;  % Number of bins (columns) after smoothing

% Load original demograph matrix
Demo_origin = load('dD_SDChannel2new.mat');
Demo_origin = Demo_origin.dsDChannel2new;

% Smooth the matrix along columns
Demo_smooth = average_colN(Demo_origin, N);

% Normalize by maximum value
Demo_normal2new = Demo_smooth ./ max(Demo_smooth(:));

% Save as image and .mat file
imwrite(Demo_normal2new, 'demo2.tif');
save('demo_sm_norm2new.mat', 'Demo_normal2new');

%% Channel 3 Processing
N = 80;

Demo_origin = load('dD_SDChannel3new.mat');
Demo_origin = Demo_origin.dsDChannel3new;

Demo_smooth = average_colN(Demo_origin, N);
Demo_normal3new = Demo_smooth ./ max(Demo_smooth(:));

imwrite(Demo_normal3new, 'demo3.tif');
save('demo_sm_norm3new.mat', 'Demo_normal3new');

%% Channel 4 Processing
N = 80;

Demo_origin = load('dD_SDChannel4new.mat');
Demo_origin = Demo_origin.dsDChannel4new;

Demo_smooth = average_colN(Demo_origin, N);
Demo_normal4new = Demo_smooth ./ max(Demo_smooth(:));

imwrite(Demo_normal4new, 'demo4.tif');
save('demo_sm_norm4new.mat', 'Demo_normal4new');

%% Convert and Save Normalized Matrix as 16-bit TIF
% ----------------------------------------------------------
% This section loads a normalized and smooth demograph matrix from a .mat file,
% rescales the data into the range of 16-bit unsigned integers (0–65535),
% and saves it as a 16-bit grayscale TIF image.
% ----------------------------------------------------------

% 1. Load the normalized matrix from .mat file
fileName = 'demo_sm_norm2new.mat';  % Replace with the actual file name if different
data = load(fileName);              % Load the contents of the .mat file

% 2. Extract the matrix
% Replace 'Demo_normal2new' with the actual variable name inside the .mat file if necessary
matrixData = data.Demo_normal2new;

% 3. Convert matrix to 16-bit range (0–65535)
% First normalize to [0,1], then scale to 16-bit and convert
matrixData_16bit = uint16(matrixData / max(matrixData(:)) * 65535);

% 4. Define output file name and save as 16-bit TIF
outputFileName = 'matrix_demo_sm_norm2new_16bit.tif';
imwrite(matrixData_16bit, outputFileName);

% 5. Display confirmation message
disp([' 16-bit TIF image has been saved as: ', outputFileName]);

