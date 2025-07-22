function TrackRefine = traceInROI(tracksFinal, ROIs, Pixelsize, ExpT, InterT, ImageBF, TraceName, ImageFL)
% traceInROI  Filter and process trajectories by assigning them to specified ROIs and calculating trajectory parameters.
%
% This function filters all trajectories in tracksFinal to find those belonging to any ROI (detected by Cellpose and formatted by roitextRead).
% Trajectories outside any ROI or near the edges of the field of view are removed.
% It calculates the mean squared displacement (MSD), diffusion coefficients, rotates trajectories to align with the ROI principal axis, 
% and crops corresponding bright-field (BF) and fluorescence (FL) images.
%
% Inputs:
%   - tracksFinal : Structure array containing all trajectories, typically obtained after spotLinking.m
%   - ROIs        : Cell array of ROI structures, formatted by roitextRead function
%   - Pixelsize   : Pixel size in nanometers (nm), used for coordinate conversion
%   - ExpT        : Exposure time in milliseconds (ms), used in MSD calculation
%   - InterT      : Interval time between frames in milliseconds (ms), used in MSD calculation
%   - ImageBF     : Bright-field image matrix for cropping ROI regions
%   - TraceName   : String used to name each trajectory
%   - ImageFL     : (Optional) Fluorescence image matrix for cropping ROI regions
%
% Outputs:
%   - TrackRefine : Cell array containing filtered trajectories with detailed info, 
%                   including rotated trajectories, MSD, diffusion coefficients, and cropped images.
%
% Notes:
%   - Trajectories must have at least some points inside an ROI.
%   - ROIs near the image edge are excluded to ensure complete cropping.
%   - Trajectories are rotated based on the PCA rotation matrix of their assigned ROI.
%   - MSD is calculated using the first 3 points of each trajectory and fitted linearly to estimate diffusion coefficient.
%   - Single-step displacements and diffusion coefficients are also calculated.
%   - Each trajectory is assigned a unique TraceId combining TraceName and current timestamp.
%
% Usage:
%   This function is typically used immediately after spotLinking.m to refine and analyze trajectories.
%
% Dependencies:
%   - MSDsingle2D (function to compute single trajectory 2D MSD)
%
% Author: Xinxing Yang
% Date: 2024-09-11 @ USTC

%% Input argument checks and default values
if nargin < 8
    ImageFL = [];  % Default: no fluorescence image
end
if nargin < 7
    TraceName = 'tj';   % Default trajectory name
end
if nargin < 6
    error('At least 6 inputs required: tracksFinal, ROIs, Pixelsize, ExpT, InterT, ImageBF');
end

%% Get all ROI center coordinates
ROIcenterCoord = zeros(length(ROIs), 2);
for jj = 1:length(ROIs)
    ROIcenterCoord(jj,:) = ROIs{jj,1}.CenterXY;
end

%% Loop through each trajectory, assign ROI and process
trindex = 1; % Index for valid trajectories
for ii = 1:length(tracksFinal)
    trackTemp = tracksFinal(ii).Coord; % Original trajectory coordinates
    
    % Convert coordinates to pixel units (add 0.5 for pixel center correction)
    trackTemp(:,2) = trackTemp(:,2)/Pixelsize + 0.5;
    trackTemp(:,3) = trackTemp(:,3)/Pixelsize + 0.5;
    
    % Calculate trajectory center coordinate
    CenterCoord = [mean(trackTemp(:,2)), mean(trackTemp(:,3))];
    
    % Count how many trajectory points fall inside each ROI
    Innumber = zeros(length(ROIs),1);
    for kk = 1:length(ROIs)
        Outline = ROIs{kk,1}.ROIsubpixel;
        Inindex = inpolygon(trackTemp(:,2), trackTemp(:,3), Outline(:,1), Outline(:,2));
        Innumber(kk) = sum(Inindex);
    end
    
    % Find ROI with maximum number of points inside
    [Maxnumber, ROIindex] = max(Innumber);
    
    % Only process if trajectory has points inside any ROI
    if Maxnumber > 0
        ROItemp = ROIs{ROIindex,1};
        Xroi = round(ROItemp.CenterXY(1));
        Yroi = round(ROItemp.CenterXY(2));
        CropSize = round(max(ROItemp.Cellsize));
        L = ROItemp.Cellsize(1);
        W = ROItemp.Cellsize(2);
        Mpca = ROItemp.PCAcoeff;
        
        % Exclude ROIs near image edges to ensure valid cropping
        if (Xroi - CropSize > 0) && (Xroi + CropSize < size(ImageBF,2)) && ...
           (Yroi - CropSize > 0) && (Yroi + CropSize < size(ImageBF,1))
       
            % Center and rotate trajectory relative to ROI
            trackRecenter = [trackTemp(:,2)-Xroi, trackTemp(:,3)-Yroi];
            trackRotate = trackRecenter * Mpca;
            trackRotateNm = trackRotate * Pixelsize;
            trackRotatePer = [trackRotate(:,1)/L*2, trackRotate(:,2)/W*2];
            
            % Calculate MSD for the trajectory
            s_MSD = MSDsingle2D(trackTemp, Pixelsize/1000, ExpT/1000, InterT/1000);
            if size(s_MSD,1) >= 3
                p = polyfit(s_MSD(1:3,1), s_MSD(1:3,2), 1);
                MSDFit = polyval(p, s_MSD(:,1));
                D_traj = p(1)/4;          % Diffusion coefficient (μm^2/s)
                C_traj = sqrt(p(2)) * 1000; % Intercept in nm
                s_MSD(:,5) = MSDFit;
            else
                D_traj = NaN;
                C_traj = NaN;
            end
            
            % Calculate single-step displacement diffusion coefficients
            TimeDiff = diff(trackTemp(:,1)) * (ExpT + InterT) / 1000; % in seconds
            Disp_single = sqrt(diff(trackRotateNm(:,1)).^2 + diff(trackRotateNm(:,2)).^2);
            D_single = Disp_single.^2 / 4 / 1e6 ./ TimeDiff;
            
            % Concatenate updated trajectory data
            trackupdate = [trackTemp, trackRotate, trackRotateNm, trackRotatePer, [TimeDiff;0], [Disp_single;0], [D_single;0]];
            
            % Save processed trajectory info
            TrackRefine{trindex,1}.TrackOriginal = trackTemp;
            TrackRefine{trindex,1}.Coord = trackupdate;
            TrackRefine{trindex,1}.TrackCenter = [CenterCoord; mean(trackRotate); mean(trackRotateNm); mean(trackRotatePer)];
            TrackRefine{trindex,1}.s_MSD = s_MSD;
            TrackRefine{trindex,1}.D_msd3 = D_traj;
            TrackRefine{trindex,1}.C_msd = C_traj;
            
            % Center and rotate ROI coordinates
            ROIcenter = [ROItemp.ROIsubpixel(:,1)-Xroi, ROItemp.ROIsubpixel(:,2)-Yroi];
            ROIRecRecenter = [ROItemp.ROIrectan(:,1)-Xroi, ROItemp.ROIrectan(:,2)-Yroi];
            ROIRotate = ROIcenter * Mpca;
            ROIRecRotate = ROIRecRecenter * Mpca;
            
            TrackRefine{trindex,1}.ROIOrigin = ROItemp;
            TrackRefine{trindex,1}.ROIcenter = ROIcenter;
            TrackRefine{trindex,1}.ROIRecRecenter = ROIRecRecenter;
            TrackRefine{trindex,1}.ROIRotate = ROIRotate;
            TrackRefine{trindex,1}.ROIRecRotate = ROIRecRotate;
            
            % Crop and rotate BF and FL images using PCA transform
            tform = affine2d([Mpca [0;0]; 0 0 1]);
            
            if ~isempty(ImageFL)
                ImFLtemp = ImageFL(Yroi-CropSize:Yroi+CropSize, Xroi-CropSize:Xroi+CropSize);
                TrackRefine{trindex,1}.ROIFL = ImFLtemp;
                G_FL = mode(ImFLtemp(:));
                ImFLtempRot = imwarp(ImFLtemp, tform, 'cubic', 'FillValues', G_FL);
                TrackRefine{trindex,1}.ROIFLrot = ImFLtempRot;
            end
            
            ImBFtemp = ImageBF(Yroi-CropSize:Yroi+CropSize, Xroi-CropSize:Xroi+CropSize);
            TrackRefine{trindex,1}.ROIBF = ImBFtemp;
            G_BF = mode(ImBFtemp(:));
            ImBFtempRot = imwarp(ImBFtemp, tform, 'cubic', 'FillValues', G_BF);
            TrackRefine{trindex,1}.ROIBFrot = ImBFtempRot;
            
            % Generate unique TraceId using TraceName and current timestamp
            currentTime = datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS');
            y = mod(year(currentTime), 100);
            m = month(currentTime);
            d = day(currentTime);
            h = hour(currentTime);
            mi = minute(currentTime);
            sec = second(currentTime);
            millisecond = mod(sec*1000, 1000);
            Strtime = sprintf('%02d%02d%02d%02d%02d%02d%02d', y, m, d, h, mi, round(sec), millisecond);
            TrackRefine{trindex,1}.TraceId = [TraceName Strtime];
            
            % Save constants for reference
            TrackRefine{trindex,1}.Pixelsize = Pixelsize;
            TrackRefine{trindex,1}.ExpT = ExpT;
            TrackRefine{trindex,1}.InterT = InterT;
            
            trindex = trindex + 1;
        end
    end
end

disp('All the trajectories found their home!');
end
