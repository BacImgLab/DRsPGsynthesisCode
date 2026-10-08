%% Construct septal intensity profiles for multi-channel images of D. radiodurans stages 3-5
% This script processes fluorescence microscopy images of D. radiodurans cells at division stages 3 to 5.
% It performs the following steps:
% 1. User selects a folder containing segmented multi-channel TIFF images.
% 2. For each image:
%    - Reads the brightfield channel to compute the rotation angle and bounding box of the cell.
%    - Rotates all fluorescence channels to horizontally align the cell.
%    - Displays the rotated multi-channel image stack for user quality control and acceptance.
%    - If accepted, uses the membrane channel to interactively select septum positions:
%         * S0 (primary septum)
%         * S1 left and S1 right (secondary septa)
%    - Extracts the intensity profiles along each septum (S0, S1L, S1R) for all fluorescence channels,
%      using a specified linewidth for averaging.
%    - Calculates septal geometric parameters including:
%        * D (cell width or septal length based on bounding box or profile length),
%        * d (septal profile length),
%        * d/D ratio indicating septal constriction.
%    - Saves rotated images, septal positions, profiles, and geometric data for further analysis.
%
% Parameters:
%   Wline        - Width in pixels of the line used for averaging intensity profiles.
%   ChannelNum   - Number of image channels in the TIFF files.
%   Channelbf    - Index of the brightfield channel for rotation angle calculation.
%   ChannelMem   - Index of the membrane channel for septum selection and profile extraction.
%
% Outputs:
%   - Rotated multi-channel TIFF images saved with suffix '_rot.tif'.
%   - Data files (.mat) containing septal positions, rotation info, geometric parameters, and profiles.
%   - CSV files for septal intensity profiles of S0, S1L, and S1R.
%
% Dependencies:
%   - bf2rotateRec: calculates rotation angle and bounding box from brightfield image.
%   - DR4cShow: interactive GUI for user to select usable cells.
%   - S0S1select: interactive tool to select S0 septum position,S1 left and right septa positions..
%   - lineProfile: extracts intensity profile along a specified line with averaging width.
%
% This script supports detailed quantitative analysis of septal morphology and fluorescence intensity 
% distributions during late division stages (3-5) of D. radiodurans cells.

clear; clc;

% 1. setup the parameters
Wline = 7; % the linewidth for construct the intensity profile
ChannelNum = 4; % number of channels

Channelbf = 1; % the channel of brightfield
ChannelMem = 2; % the channel of membrane

% 2. select folder(or folders) for processing
selpath = uigetdir('','Select the folder contains segmented cells');
resultFolder = fullfile(selpath, 'Processed');
if ~exist(resultFolder, 'dir')
    mkdir(resultFolder);  % Create the folder if it doesn't exist
end

% Find all .tif files in the source folder
tifFiles = dir(fullfile(selpath, '*.tif'));

% 3. walk through all the images in the folder
for idf = 1 : numel(tifFiles)
    imageFilename = tifFiles(idf).name;
    [~, name, ext] = fileparts(imageFilename);        % get the name without extension
    imageFilepath = fullfile(selpath,imageFilename);  % full path for reading the original images
    newImgName = strcat(name, '_rot', '.tif');  
    newDataName = strcat(name, '_data', '.mat');
    newprofName1 = strcat(name, '_Sd', '.csv');
    newprofName2 = strcat(name, '_Sl', '.csv');
    
    imageSavename = fullfile(resultFolder,newImgName);  % path+filename for saving rotated images
    dataSavename = fullfile(resultFolder,newDataName);  % path+filename for saving septum positions and metadata
    profileSavename1 = fullfile(resultFolder,newprofName1);  % path+filename for saving intensity profile along Sd
    profileSavename2 = fullfile(resultFolder,newprofName2);  % path+filename for saving intensity profile along Sl
  
    % read brightfield channel and find the rotation angle and bounding box
    imgbf = imread(imageFilepath,Channelbf);
    [angle,rotatedBox] = bf2rotateRec(imgbf);

    % 4. read all the channels and rotate to align horizontally
    for ii = 1 : ChannelNum
        imgTemp = imread(imageFilepath,ii);
        imgRot = imrotate(imgTemp, -angle, 'bilinear', 'crop');
        ImageAllChannel(:,:,ii) = imgRot; % store rotated images for all channels
    end

    % 5. determine whether the current cell image is usable or not
    imgFlag = DR4cShow(ImageAllChannel); % interactive checking of the cell

    % 6. if cell is good, get the positions of S0 and S1 septa
    if imgFlag == 1
        % use membrane channel as reference for septa positions (lines)
        ImShowing = ImageAllChannel(:,:,ChannelMem);
        [xS0, yS0, xS1l, yS1l, xS1r, yS1r] = S0S1select(ImShowing,rotatedBox);

        % initialize profiles for septa
        Sdprofile = [];
        Slprofile = [];
        S1rprofile = [];

        % get the intensity profile for each septum in all channels
        for ii = 1 : ChannelNum
            imgC = ImageAllChannel(:,:,ii); % current channel image
            
            % profile along septum 0 (S0)
            [S0p,lineX0,lineY0,Lend] = lineProfile(imgC,xS0,yS0,Wline);
            % profile along septum 1 left (S1-left)
            [S1lp,lineX1,lineY1,Lenl] = lineProfile(imgC,xS1l,yS1l,Wline);
            % profile along septum 1 right (S1-right)
            [S1rp,lineX2,lineY2,Lenr] = lineProfile(imgC,xS1r,yS1r,Wline);
            
            Sdprofile = [Sdprofile, S0p'];
            Slprofile = [Slprofile, S1lp'];
            S1rprofile = [S1rprofile, S1rp'];
            
            % use membrane channel to calculate septal geometry D, d, and ratio
            if ii == ChannelMem
                % For S0, D is the cell width (rotatedBox width), d is the length of S0 profile
                ds0 = Lend;
                Ds0 = rotatedBox(3);
                DratioS0 = ds0/Ds0;
                
                % For S1 left, D is the length of the septum line
                % Calculate d based on intensity thresholding on profile
                MinS1lp = min(S1lp(round(length(S1lp)*0.25):end-round(length(S1lp)*0.25)));
                MaxS1lp = max(S1lp(1:round(length(S1lp)*0.5)));
                ThreshS1lp = 0.5*(MaxS1lp - MinS1lp) + MinS1lp; % 50% intensity threshold
                MaxS1lp2 = max(S1lp(round(length(S1lp)*0.5):end));
                % extract middle region below threshold
                S1lpTemp = S1lp(find(S1lp==MaxS1lp):find(S1lp==MaxS1lp2));
                S1lpMiddle = S1lpTemp(S1lpTemp < ThreshS1lp);
                ds1l = length(S1lpMiddle);
                DratioS1l = ds1l/Lenl;
                
                % For S1 right, calculate d similarly
                MinS1rp = min(S1rp(round(length(S1rp)*0.25):end-round(length(S1rp)*0.25)));
                MaxS1rp = max(S1rp(1:round(length(S1rp)*0.5)));
                ThreshS1rp = 0.5*(MaxS1rp - MinS1rp) + MinS1rp;
                MaxS1rp2 = max(S1rp(round(length(S1rp)*0.5):end));
                S1rpTemp = S1rp(find(S1rp==MaxS1rp):find(S1rp==MaxS1rp2));
                S1rpMiddle = S1rpTemp(S1rpTemp < ThreshS1rp);
                ds1r = length(S1rpMiddle);
                DratioS1r = ds1r/Lenr;
            end
        end
        
        % 6. save all results to a .mat file
        data.LACell = rotatedBox(3); % long axis length of the cell
        data.SACell = rotatedBox(4); % short axis length of the cell
        
        % lengths of septum lines
        data.LoS0 = sqrt(diff(xS0)^2 + diff(yS0)^2);
        data.LoS1L = sqrt(diff(xS1l)^2 + diff(yS1l)^2);
        data.LoS1R = sqrt(diff(xS1r)^2 + diff(yS1r)^2);
        
        % save septum endpoint coordinates
        data.xS0 = xS0; data.yS0 = yS0;
        data.xS1L = xS1l; data.yS1L = yS1l;
        data.xS1R = xS1r; data.yS1R = yS1r;
        
        % save full line coordinates for profiles
        data.XYS0 = [lineX0', lineY0'];
        data.XYS1L = [lineX1', lineY1'];
        data.XYS1R = [lineX2', lineY2'];
        
        % save rotated bounding box and angle
        data.rotatedBox = rotatedBox;
        data.angle = angle;
        
        % save width d, length D and ratio d/D for all septa
        data.ds0 = ds0;
        data.Ds0 = Ds0;
        data.dDratioS0 = DratioS0;
        
        data.ds1l = ds1l;
        data.Ds1l = Lenl;
        data.dDratioS1l = DratioS1l;
        
        data.ds1r = ds1r;
        data.Ds1r = Lenr;
        data.dDratioS1r = DratioS1r;
        
        % save intensity profiles for septa and all channels
        data.S0profile = Sdprofile;
        data.S1lprofile = Slprofile;
        data.S1rprofile = S1rprofile;

        save(dataSavename,'data');
        
        % save rotated images for all channels as multipage TIFF
        imwrite(ImageAllChannel(:,:,1), imageSavename, "tif", "WriteMode", "overwrite", "Compression", "none");
        for jj = 2 : ChannelNum
            imwrite(ImageAllChannel(:,:,jj), imageSavename, "tif", "WriteMode", "append", "Compression", "none");
        end
    end
end
