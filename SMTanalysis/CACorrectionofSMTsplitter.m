%% CACorrectionofSMTsplitter
%
% This script performs:
% 1. Localization point matching from ThunderSTORM output
% 2. Polynomial transformation fitting (moving → fixed)
% 3. Applies the transformation to raw images
% 4. Saves aligned images
%
clear;
clc;
%% Tform Calculation
% parametera setting
file1 = '647.csv';      % Moving channel (e.g., 647nm)
file2 = '561.csv';      % Fixed channel (e.g., 561nm)
Pixelsize = 110;        % Pixel size in nm (ThunderSTORM must use same)
Disthresh = 500;        % Max allowed distance for matched points in nm

% read localization data from ThunderSTORM in Fiji
[Data1, ~] = xlsread(file1);  % Moving (647)
[Data2, ~] = xlsread(file2);  % Fixed (561)
Coords1 = Data1(:, 2:4);      % Columns: Frame, X, Y
Coords2 = Data2(:, 2:4);

% match points frame by fram
VectorM = [];  % Initialize matrix to store matched localizations

for ii = 1:max(Coords1(:,1))  % Loop through each frame
    IndexF1 = find(Coords1(:,1) == ii);  % Indices in Coords1 for frame ii
    IndexF2 = find(Coords2(:,1) == ii);  % Indices in Coords2 for frame ii
    
    Coords1F = Coords1(IndexF1,:);  % Frame ii, 647 channel
    Coords2F = Coords2(IndexF2,:);  % Frame ii, 561 channel
    
    for jj = 1:size(Coords1F,1)  % Loop through all localizations in 647
        CoordTemp = Coords1F(jj, 2:3);  % [X, Y] position
        
        % Compute Euclidean distances to all 561 points in same frame
        Dx = Coords2F(:,2) - CoordTemp(1);
        Dy = Coords2F(:,3) - CoordTemp(2);
        Difxy = sqrt(Dx.^2 + Dy.^2);
        
        [DifMin, DminIndex] = min(Difxy);  % Find closest 561 localization
        
        if DifMin < Disthresh  % Accept only close matches
            % Store: frame, X647, Y647, X561, Y561, dx, dy, distance
            VectorTemp = [Coords1F(jj,:), Coords2F(DminIndex,2:3), Dx(DminIndex), Dy(DminIndex), Difxy(DminIndex)];
            VectorM = [VectorM; VectorTemp];
        end
    end
end

% compite geoemtric transformation
% Convert coordinates from nanometers to pixels
fixedPoints = VectorM(:, 2:3) / Pixelsize;   % 561 channel (fixed)
movingPoints = VectorM(:, 4:5) / Pixelsize;  % 647 channel (moving)

% Fit 2nd-order polynomial geometric transformation
tform = fitgeotrans(movingPoints, fixedPoints, 'polynomial', 2);

%% apply transformation to images
% Set image folder path
imageFolder = 'D:\ImageData\SMT\Data'; % Replace with your data path

% Load the fixed image (561 channel)
image1 = imread(fullfile(imageFolder, 'caclutest647.tif'));  % fixed reference

% Find all matching 'FL*.tif' files (moving images)
channel2Files = dir(fullfile(imageFolder, 'FL7.tif'));  % only one file in this example
numFiles = length(channel2Files);

for i = 1:numFiles
    % Read moving image (e.g., 647 raw)
    image2 = imread(fullfile(imageFolder, channel2Files(i).name));
    
    % Set reference view to match fixed image size
    outputView = imref2d(size(image1));
    
    % Apply geometric transformation
    alignedImage2 = imwarp(image2, tform, 'OutputView', outputView);
    
    % Create output filename with "_aligned" suffix
    [~, name, ext] = fileparts(channel2Files(i).name);
    outputFileName = fullfile(imageFolder, [name '_aligned' ext]);
    
    % Save aligned image to disk
    imwrite(alignedImage2, outputFileName);
    
    % Print progress
    fprintf('Processing image %d of %d: %s\n', i, numFiles, channel2Files(i).name);
end
