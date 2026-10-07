%% trajectoryAnalysisMain.m
%
% Main pipeline script to process original single-particle tracking data, associate
% trajectories with segmented cell ROIs, and generate summary visualization stacks.
%
% This script performs the following steps:
%   1. Load segmented cell outlines from a text file.
%   2. Load raw tracking results from tracking software (e.g., uTrack or TrackMate).
%   3. Load corresponding bright-field (BF) and fluorescence (FL) images.
%   4. Run the function `traceInROI` to associate trajectories with individual ROIs.
%   5. Save the refined trace data (`TrackRefine`) to a .mat file.
%   6. Run the `tracePlot.m` script to generate visualization figures and save them
%      as a multi-page TIFF stack.
%
% Requirements:
%   - ROIs must be exported as text outlines from cellpose.
%   - The images (BF/FL) must be spatially aligned (FL7_aligned.tif assumed to be registered).
%
% Inputs:
%   - Text file containing ROI outlines. such as 'BF7_cp_outlines.txt' 
%   - MATLAB file containing variable `tracksFinal`. such as'long-10-Coord-track7.mat'
%   - Bright-field image corresponding to the tracking. such as 'BF7.tif'
%   - Registered fluorescence image. such as 'FL7_aligned.tif'
%
% Outputs:
%   - 'TraceRefine-test.mat'     : Saved MATLAB structure for all ROIs and associated traces.
%   - 'test.tif'                 : Multi-page TIFF file with visual summaries (created by tracePlot.m).
%
% modified by by Xinxing Yang

clear
clc
% Load ROI outlines from text file
% These ROIs were segmented from BF image and saved as text
ROIs = roitextRead('BF7_cp_outlines.txt');
% Step 3: Load tracking results
load('long-10-Coord-track7.mat');
% Bright-field image (grayscale or RGB)
ImBF = imread('BF7.tif');
% Fluorescence image, aligned to BF
ImFL = imread('FL7_aligned.tif');

% Associate traces with ROIs using traceInROI
% traceInROI(tracksFinal, ROIs, timeInterval, maxTraceLength, minTraceLength, BFimage, nameTag, FLimage)
TrackRefine = traceInROI(tracksFinal, ROIs, ...
                         110, ...     % Frame interval (ms)
                         1000, ...    % Max allowed trace length
                         0, ...       % Minimum trace length
                         ImBF, ...    % BF image (for visualization)
                         'JJ', ...    % Name tag for trace IDs
                         ImFL);       % Aligned FL image

% Save the refined trace results
save('TraceRefine-test.mat', 'TrackRefine');

%Run the plotting script to generate TIFF stack
% This creates the visual summary plots and saves as test.tif
run('tracePlot.m');

