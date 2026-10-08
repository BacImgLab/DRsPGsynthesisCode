%% DemoS1Analysis.m
%
% Analyze a demograph from Sample 1 by fitting each vertical profile (column)
% with a two-component Gaussian model, in order to extract ridge positions
% for both the protein of interest and the membrane marker.
%
% For each column in the input demograph matrices:
%   - A two-Gaussian model is fitted to both protein and membrane signals
%   - The first Gaussian peak position is extracted
%   - Raw and fitted data are plotted and saved as a multi-frame TIFF
%   - All fitting results and parameters are saved in a .mat file
%
% INPUT:
%   - 'demo_S1_smooth_W.mat' : smoothed demograph for protein signal
%   - 'demo_S1_smooth_mem.mat' : smoothed demograph for membrane signal
%
% OUTPUT:
%   - 'Wmem.tif' : plot of each column's profile and its Gaussian fits
%   - 'WmemData.mat' : contains peak positions and fit parameters
%
% Author: Xinxing Yang
% Date: 2025-04-05

clear; clc; close all;

index_threshold = 30; % Threshold index (not used in this script)
filesave = 'Wmem.tif'; % Output TIFF file for fit plots
matsave = 'WmemData.mat'; % Output .mat file to store results

% Load smoothed protein demograph
data = load('demo_S1_smooth_W.mat');
fields = fieldnames(data);            
POImat = data.(fields{1}); % Extract the matrix from the loaded struct

% Load smoothed membrane demograph
data = load('demo_S1_smooth_mem.mat');
fields = fieldnames(data);            
MEMmat = data.(fields{1}); % Extract the matrix from the loaded struct

[row, col] = size(POImat); % Get size of the demograph (both matrices must match)

X = (1:row)'; % X-axis vector for profile fitting

mu1_mem = []; % Store membrane peak positions
mui_poi = []; % Store protein peak positions
ParameterAll = []; % Store all fit parameters and results

for ii = 1:col
    Imem = MEMmat(:, ii); % Membrane profile for this column
    Ipoi = POImat(:, ii); % Protein profile for this column

    % Fit both profiles using the double-Gaussian model
    [pk1mem, ~, p_fitmem, Fitresmem] = demofitGauss2(X, Imem);
    [pk1poi, ~, p_fitpoi, Fitrespoi] = demofitGauss2(X, Ipoi);

    % Store the first peak positions
    mui_poi(ii, 1) = pk1poi;
    mu1_mem(ii, 1) = pk1mem;

    % Store fitting parameters and result curves
    ParameterAll(ii, 1).p_fitmem = p_fitmem;
    ParameterAll(ii, 1).p_fitpoi = p_fitpoi;
    ParameterAll(ii, 1).Fitresmem = Fitresmem;
    ParameterAll(ii, 1).Fitrespoi = Fitrespoi;

    % Plot and save figure for this column
    h = figure('visible', 'off');
    plot(X, POImat(:, ii), 'ob'); hold on; % Raw protein profile
    plot(X, MEMmat(:, ii), 'dr'); % Raw membrane profile
    plot(Fitrespoi(:, 1), Fitrespoi(:, 2), '-b'); % Fitted protein curve
    plot(Fitresmem(:, 1), Fitresmem(:, 2), '-r'); % Fitted membrane curve
    title(num2str(ii));
    legend({'Protein', 'Membrane'});

    % Capture figure as image
    frame = getframe(h);
    imgplot = frame2im(frame); % Convert to RGB image
    [imind, cm] = rgb2ind(imgplot, 256); % Convert to indexed image

    % Save or append to multi-page TIFF
    if ii == 1
        imwrite(imgplot, filesave, 'tif', 'WriteMode', 'overwrite', 'Compression', 'none');
    else
        imwrite(imgplot, filesave, 'tif', 'WriteMode', 'append', 'Compression', 'none');
    end

    close(h); % Close figure to free memory
    ii % Display current column index
end

figure('visible', 'on'); % Open a figure (optional GUI reset)
close all; % Close all figures

% Save all fitting results and peak positions
save(matsave, 'mui_poi', 'mu1_mem', 'ParameterAll');
