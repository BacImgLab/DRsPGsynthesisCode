%% statesClassifyDr
%
% This script is designed to classify the motion states of segmented trajectory data,
% with the goal to categorize segments into Directional motion, Diffusive motion, and Stationary states.
%
% Main functionalities include:
%  1. Simulating trajectory data within given boundaries and diffusion coefficients (optional)
%  2. Preprocessing and merging segmented trajectory information
%  3. Calculating directional motion probabilities (P_directional) and related statistics
%  4. Visualizing trajectories and their motion parameters
%  5. Classifying trajectory segments based on threshold values into motion states
%  6. Statistical analysis and plotting of velocity and dwell time histograms
%  7. Optional calculation and fitting of mean square displacement (MSD) for stationary trajectories
%
% Dependencies and required functions:
%  - rcdfCal: trajectory simulation
%  - mergeSegDataJ: merging segmented trajectory data
%  - segUpdate: computation of motion parameters and probabilities for segments
%  - segPRplot: visualization of trajectories and parameters
%  - MSDcalculate_2d, kusumi_xy: MSD calculation and fitting functions
%  - histLog_xy: plotting velocity histograms
%
% Input data:
%  - Pre-segmented trajectory data (segmentation performed using RefineTraceSegDr.APP)
%  - Simulated trajectory data (if applicable, generated in Section 0)
%  - mergedAllSeg.mat: merged structure of segmented trajectory data
%
% Output:
%  - SegUpdated.mat: structure containing updated motion probabilities and parameters
%  - Statistical data of motion states (velocity, dwell time, etc.)
%  - Various visualization plots (trajectory maps, velocity distributions, dwell time histograms, MSD fitting curves)
%
% Usage instructions:
%  - Ensure trajectories are segmented prior to running this script
%  - Adjust simulation parameters and thresholds according to experimental setup
%  - Run each section sequentially to complete the analysis
%  - Manual verification of threshold values and classification accuracy is recommended
%
% Author: Adapted from code by Xinxing Yang
% Version: 1 (2024-11-20)


%% Section 0: simulate a set of trajectories with a certain boundary and diffusion coeff 
%         1.This is not a neccessary step for every experiment. As long as the boundary does not change too much, you can use the same file generated in this step for other experiments.
%         2. A way to verify your simulation is to calculate the MSD of your stationary molecules in the end (section X). If the boundary and D are not so different from your setting here, it should be okay.
%         3. In this simulation,we consider the velocity = 0.
% Set the hyperparameters for your simulation, simulate, and save the file
clear; clc;
Frame_L = [3:200]; % the range of the possible trajectory length
ExpT = 1; % the exposure time (if there is dark interval, this should be total time interval)
D = 0.000050; % diffusion coefficient: in um^2/s
B = 70; % boundary size in nm (in DR，the B should bigger)
L_err = 15; % localization error in nm
N_traj = 2000; % number of trajectories for simulation in one condition
filenameSimu = 'SimuDR_2.mat'; % filename to save the simulation result
% simulation
[R_struc,Traj_struc,TimeMatrix,frameMatrix,SpeedMatrix] = rcdfCal(Frame_L,0,ExpT,D,B,L_err,N_traj);
save(filenameSimu,'R_struc','Traj_struc','TimeMatrix','frameMatrix','SpeedMatrix','Frame_L');

%% Section 1: Data preprocessing
%The purpose of this section is to organize the information of the seg single track into a new file, mergedAllSeg
%The purpose of this section is to calculate the orientation probability and the rest probability of the segmented trajectory

mergeSegDataJ(); % The information of the segmented single trace is integrated

%% Section 2: Calculate the Pdirect of the trajectory
%  set up bootstrapping parameters
Nboot = 100;    % number of the bootstrapping to get linear fitting
Pdrop = 0.1;    %dropout probability in bootstrapping

load('mergedAllSeg.mat'); 
load('SimuDR_2.mat');
mergedAllSeg = segUpdate(mergedAllSeg,Nboot,Pdrop,Traj_struc, R_struc, Frame_L, TimeMatrix); % calculate the R and P value of each segmets

% updata mergedAllSeg 
save('SegUpdated.mat', 'mergedAllSeg');
disp('The probobalitly of all trajectories being moving are calculated by Bootstrapping. Saved in SegUpdated.mat');

%% section 3 (optional): plot the trajectories one-by-one and save the images for checking. 
clear
clc
load('SegUpdated.mat');
segPRplot(mergedAllSeg);

%% Section 4: Classify the segments by their R and P
% note: the section might need to be done multiple times with human
% interference to make sure the threshold is reasonble and can be applied
% to different conditions
clear;
clc;

% Set the classification thresholds % this threshold is based on the Dr
% PBP1B and FtsW data
Rmax1 = 0.5;  % Maximum R value for directional tracks
Rmax2 = 0.3;   % Threshold for long directional segments
StDmax = 100;   % Threshold for standard deviation (stationary phase)
Pmin = 0.80;   % Minimum probability

% Set the parameters for histogram plotting
Nbin = 32;
LowB = -1.5;
HighB = 3;

% Load the finalAllSeg file (replace with the actual filename)
load('SegUpdated.mat');  % Replace with actual file path

% record the V, StD, R, P value, dwelltime, and the class of each segment
SegVSRPTC = [];
for i = 1:length(mergedAllSeg)
   
    % Extract the struct from the cell and access TrackData
    TrackData = mergedAllSeg{i};  % Accessing the struct inside the cell
    % Iterate over each track in TrackData
    
        % Extract relevant parameters for each track
        R = TrackData.Rxy;    % Ratio (R value)
        P = TrackData.P_directional;    % P-progressive
        StD = TrackData.StDXY; % Standard deviation
        V = TrackData.OverallSpeed;  % Velocity
        DT = TrackData.Length; % Dwell time
        Loc = TrackData.Location; % Location

        % Classification based on thresholds (using & and | for array compatibility)
        if (R <= Rmax2) | ((R > Rmax2) & (R < Rmax1) & (P >= Pmin))
            % Directional track
            StateC = 1; % Set state to 1 (Directional)
        elseif StD > StDmax
            % Diffusive track
            StateC = 2; % Set state to 2 (Diffusive)
        else
            % Stationary track
            StateC = 3; % Set state to 3 (Stationary)
        end
    
    SegVSRPTC(i,:) = [V,StD,R,P,DT,Loc,StateC];
end

% generate the histogram, V, and dwell time
Vd = SegVSRPTC(find(SegVSRPTC(:,7)==1),1); % find the speed of all the directional moving segments
Vs = SegVSRPTC(find(SegVSRPTC(:,7)==3),1); % find the speed of all the stationary segments
Td = SegVSRPTC(find(SegVSRPTC(:,7)==1),5); % find the dwell time of all the directional moving segments
Ts = SegVSRPTC(find(SegVSRPTC(:,7)==3),5); % find the dwell time of all the stationary segments

% calculate the histogram of dwell time of the directional and diffusive
[HisDx,HisDxx] = histcounts(Td,[0:2:100]);
[HisDf,HisDfx] = histcounts(Ts,[0:2:100]);

% calculate the histogram of speed

% StrainName = 'MG1655-FtsW^C-PBP1B^C-MTSES';
[HisRes,HisX,HisObj] = histLog_xy(Vd,LowB,HighB,Nbin);

% plot the R and P scatter 
h = figure('position',[100,100,1600,400]);
subplot(1,3,1)
semilogx(SegVSRPTC(:,3),SegVSRPTC(:,4),'.','markersize',15)
ylim([0,1.1])
xline(Rmax2,'r',['R2 =' num2str(Rmax2)],'LabelVerticalAlignment','bottom')
xline(Rmax1,'k',['R1 =' num2str(Rmax1)],'LabelVerticalAlignment','bottom')
yline(Pmin,'r',['P =' num2str(Pmin)])
yline(0.5,':k','50% line')
xlabel('Ratio')
ylabel('P-progressive')
set(gca,'fontsize',18)
% threshold the data
hold on
semilogx(SegVSRPTC(find(SegVSRPTC(:,7)==1),3),SegVSRPTC(find(SegVSRPTC(:,7)==1),4),'+','markersize',15)
semilogx(SegVSRPTC(find(SegVSRPTC(:,7)==3),3),SegVSRPTC(find(SegVSRPTC(:,7)==3),4),'d','markersize',10)

subplot(1,3,2)
% xlim([0.6,300])
% set(gca,'fontsize',14);
% fit the two log-normal distribution
BinEdges = HisObj.BinEdges;
BinCenter = BinEdges(1:end-1) + diff(BinEdges);
LogCenter = log10(BinCenter);
% h1 = figure('PaperUnits','inches','PaperPosition',[0 0 12 4],'PaperSize',[12 4])

hold on
bar(LogCenter,HisRes,1,'w','LineWidth',2)
% [fitresult, gof] = twoGaussFit(LogCenter, HisRes);
xlim([-1,2.5])
ylim([0,0.35])
xticks([0,1,2,3])
set(gca,'xticklabel',[1, 10, 100, 1000]);
set(gca,'fontsize',18);
xlabel('Speed (nm/s)','fontsize',24);
ylabel('Probability','fontsize',24);
grid on
title(['Vmean = ' num2str(mean(Vd),2) '; Vmedian = ' num2str(median(Vd),2)]);


subplot(1,3,3)
plot(HisDxx(1:end-1) + diff(HisDxx)/2,HisDx,'g','linewidth',3);
hold on
plot(HisDfx(1:end-1) + diff(HisDxx)/2,HisDf,'r','linewidth',3);
xline(mean(Td),'g',num2str(mean(Td),2));
xline(mean(Ts),'r',num2str(mean(Ts),2));
title(['Percentage of moving = ' num2str(length(Vd)/(length(Vd)+length(Vs))*100,3) '%']);
set(gca,'fontsize',18);
legend('Progressive','Stationary')
xlabel('Dwell Time/s')
% save all the results
[filenameS pathnameS] = uiputfile('.mat','Save the statistics to...');
save([pathnameS filenameS],'Rmax1','Rmax2','Pmin','SegVSRPTC','Vd','Vs','Td','Ts');
saveas(h, 'combined-plot.fig')
waitforbuttonpress
close all
output = [mean(Vd),mean(Vs),mean(Td),mean(Ts)]

%% Section 5 (optional): combine all the stationary ones and calculate the MSD, fit with kusumi
% This step is just for classification. Not neccessary for every time.
% find all the traces with R>0.25
clear;clc;
R0 = 0.3; % set the threshold
PixelS = 1; % nm
ExpT = 1;   % s
T_cut = 20; % the upper bound for fitting msd in sec

% Load the finalAllSeg file (replace with the actual filename)
load('SegUpdated.mat');  % Replace with actual file path

% Initialize the sMSD structure and the counter
CountIndex = 1;

% Iterate over each cell in mergedAllSeg
kk = 1; % Counter for the new structure
for i = 1:length(mergedAllSeg)
    % Extract the struct from the cell and access TrackData
    TrackData = mergedAllSeg{i};  % Accessing the struct inside the cell
    % Iterate over each track in TrackData
    for j = 1:length(TrackData)
        % Extract relevant parameters for each track
        R = TrackData(j).Rxy;    % Ratio (R value)
        T = TrackData(j).TimeT_TXY;    % time
        X = TrackData(j).TraceTx; % x-coordinates
        Y = TrackData(j).TraceTy; % y-coordinates

        % find the stationary ones
        if R > R0
            Dtrace = [];
            Dtrace(:,1) = T-T(1);
            Dtrace(:,2) = X;
            Dtrace(:,3) = Y;
            strucD.TracksROI(CountIndex).Coordinates = Dtrace;
            CountIndex = CountIndex + 1;
        end
    end
end

% calculate MSD
s_MSD = MSDcalculate_2d(strucD,PixelS,ExpT,0,'all');
% kusumi and plot
T_cut = 20
index = find(s_MSD(:,1) <= T_cut);
msd_2d = s_MSD(index,:);
KusumiFit =  kusumi_xy(msd_2d);

figure
errorbar(msd_2d(:,1),sqrt(msd_2d(:,2)),sqrt(msd_2d(:,4)),'ob','linewidth',1.5,'Color',[0 0 0.5430]);
hold on
plot(KusumiFit.D_fit(:,1),sqrt(KusumiFit.D_fit(:,2)),'-b','linewidth',1.5,'Color',[0 0 0.5430]);
legend({'MSD','Fit'},'Box','Off','Location','northwest');
ylabel('$\sqrt{MSD}$ (nm)','Interpreter','Latex')
xlabel('LagTime /s')
text(0.5, 0.2, {['D = ' num2str(KusumiFit.D) ' nm^2/s'],['L = ' num2str(KusumiFit.L) ' nm'],...
    ['sigma = ' num2str(KusumiFit.sigma) ' nm']}, 'Units', 'normalized', 'FontSize', 12, 'Color', 'black');
% 保存图像
saveas(gcf, 'MSD_Fit_Plot.png');


%% functions used in the code
function mergedAllSeg = segUpdate(mergedAllSeg,Nboot,Pdrop,Traj_struc, R_struc, Frame_L, TimeMatrix)
for i = 1:numel(mergedAllSeg)
    currentTrackData = mergedAllSeg{i}; 
    if isfield(currentTrackData, 'TimeT_TXY') && ...
       isfield(currentTrackData, 'TraceTx') && ...
       isfield(currentTrackData, 'TraceTy') && ...
       isfield(currentTrackData, 'Length')

        TimeT_TXY = currentTrackData.TimeT_TXY;  % time
        TraceTx = currentTrackData.TraceTx;  % X Coords
        TraceTy = currentTrackData.TraceTy;  % Y Coords
        Length = currentTrackData.Length;  % trace length
        Vx = currentTrackData.Vx;
        Vy = currentTrackData.Vy;

        % Bootstraping 
        [V, Displb, StDb, Ratiob, Nb] = tracedropoutXY(TimeT_TXY,TraceTx,TraceTy,Nboot,Pdrop);  % bootstrap the data 
        Vboot(1,1) = mean(V);
        Vboot(1,2) = std(V);
        Displboot(1,1) = mean(Displb);
        Displboot(1,2) = std(Displb);
        StDboot(1,1) = mean(StDb);
        StDboot(1,2) = std(StDb);
        Ratioboot(1,1) = mean(Ratiob);
        Ratioboot(1,2) = std(Ratiob);

        % start the probabiliry calculation
        Ratio_region = [Ratioboot(1,1) - Ratioboot(1,2), Ratioboot(1,1) + Ratioboot(1,2)];  % Ratio range
        frameL = Length;  
        [R_V, R_0, Traj_V, Traj_0] = addVoneFLtrajs(Traj_struc, R_struc, Frame_L, frameL, TimeMatrix, Vboot(1,1));  % R disturbation
        % calculate the probability of each mode
        P_0 = getProbR(R_0, Ratio_region);  
        P_V = getProbR(R_V, Ratio_region);  
        P_directional = P_V / (P_0 + P_V);  

        % save data
        mergedAllSeg{i}.V = sqrt(Vx^2+Vy^2);
        mergedAllSeg{i}.Vboot = Vboot;
        mergedAllSeg{i}.Displboot = Displboot;
        mergedAllSeg{i}.StDboot = StDboot;
        mergedAllSeg{i}.Ratioboot = Ratioboot;
        mergedAllSeg{i}.P_0 = P_0;
        mergedAllSeg{i}.P_V = P_V;
        mergedAllSeg{i}.P_directional = P_directional;
    end
end
end

%%
function segPRplot(mergedAllSeg)

filename = uiputfile('.tif','Save the segments and their R and P in a tif file');
for ii = 1 : size(mergedAllSeg,1) % loop all the trajectories
    %ii=1;
    TrackData = mergedAllSeg{ii,1};
    TimeT_TXY = TrackData.TimeT_TXY;  % time
    TraceTx = TrackData.TraceTx;  % X Coords
    TraceTy = TrackData.TraceTy;  % Y Coords
    Length = TrackData.Length;  % trace length
    Vx = TrackData.Vx;
    Vy = TrackData.Vy;
    V = TrackData.V;
    MSD = TrackData.MSD;
    P_directional = TrackData.P_directional;
    Displboot = TrackData.Displboot;
    StDboot = TrackData.StDboot;
    Ratioboot = TrackData.Ratioboot;   

    % initiate the figure for plot
    h = figure('visible', 'off');
    set(h, 'position',[100,100,1000,600]);
    
    % subplot 1-topleft: The trajectory along X
    hs1 = subplot('position',[0.08,0.58,0.41,0.38]);
    % fit the X curve and plot
    plot(TimeT_TXY, TraceTx, '-ob','LineWidth',2);  % Scatter plot of x and y
    hold on;
    %    Fit a linear model (y = mx + b)
    p = polyfit(TimeT_TXY, TraceTx, 1);  % p(1) is the slope, p(2) is the intercept
    % Create the fitted line
    x_fit = polyval(p, TimeT_TXY);
    plot(TimeT_TXY, x_fit, 'r-', 'LineWidth',3);  % Plot the fitted line
    % Label axes and add a title
    xlabel('Time/s');  % Change this to your desired X-axis label
    ylabel('X-position/nm');  % Change this to your desired Y-axis label
    title('X vesus time');
    % Display fit parameters on the plot
    fit_eq = sprintf('vx = %.2f nm/s', Vx);  % Equation of the fit
    text(0.05, 1.05, fit_eq, 'Units', 'normalized', 'FontSize', 12, 'Color', 'black');

    % subplot 2-buttomleft: The trajectory along Y
    hs2 = subplot('position',[0.08,0.08,0.41,0.38]);
    % fit the Y curve and plot
    plot(TimeT_TXY, TraceTy, '-ob', 'LineWidth', 2);  % Scatter plot of x and y
    hold on;
    %    Fit a linear model (y = mx + b)
    p = polyfit(TimeT_TXY, TraceTy, 1);  % p(1) is the slope, p(2) is the intercept
    % Create the fitted line
    y_fit = polyval(p, TimeT_TXY);
    plot(TimeT_TXY, y_fit, 'r-', 'LineWidth', 3);  % Plot the fitted line
    % Label axes and add a title
    xlabel('Time/s');  % Change this to your desired X-axis label
    ylabel('Y-position/nm');  % Change this to your desired Y-axis label
    title('Y vesus time');
    % Display fit parameters on the plot
    fit_eq = sprintf('vy = %.2f nm/s', Vy);  % Equation of the fit
    text(0.05, 1.05, fit_eq, 'Units', 'normalized', 'FontSize', 12, 'Color', 'black');
    
    % subplot 3-topright: The MSD curve
    hs3 = subplot('position',[0.55,0.58,0.41,0.38]);
    hold on
    errorbar(MSD(:,1),MSD(:,2),MSD(:,4),'linewidth',2);
    set(gca, 'FontSize', 8);
    xlabel('Lag Time (s)', 'FontSize', 9);
    ylabel('MSD \mum^2/s', 'FontSize', 9);
    % Display fit parameters on the plot
    Traceid = sprintf('Traceid = %d; v = %.2f nm/s', ii,V);  % Equation of the fit
    text(0.15, 1.05, Traceid, 'Units', 'normalized', 'FontSize', 12, 'Color', 'black');
    % subplot 4-bottomright: The ROI with trajectory color coded with time
    hs4 = subplot('position',[0.55,0.08,0.41,0.38]);
    hold on
    % plot the median trajectory in color copied from Josh-time
    colorCodeTracePlot(TimeT_TXY, [TraceTx,TraceTy],0.8)
    plot(TraceTx,TraceTy,'ok','markersize',1.5)
    plot(x_fit,y_fit,'r-', 'LineWidth',3)
    set(gca, 'FontSize', 8);
    % title(TraceID, 'FontSize', 10);
    axis equal
    xlabel('X', 'FontSize',9);
    ylabel('Y', 'FontSize', 9);
    % Display fit parameters on the plot
    R_P = sprintf('R = %.2f; P = %.2f; Sd = %.2f; D = %.2f', Ratioboot(1), P_directional, StDboot(1), Displboot(1));  % Equation of the fit
    text(0.15, 1.05, R_P, 'Units', 'normalized', 'FontSize', 12, 'Color', 'black');
    % save the figure to tiff file
    % Capture the current figure as an image
    frame = getframe(h);
    imgplot = frame2im(frame); % Convert frame to image
    
    % Convert to indexed image with colormap for TIFF compatibility
    [imind, cm] = rgb2ind(imgplot, 256);
    saveas(h,'test1.tif','tif');
    % Save the image as a multi-page TIFF
    if ii == 1
        % For the first image, create the TIFF file
        imwrite(imgplot, filename, 'tif', 'WriteMode', 'overwrite', 'Compression', 'none');
    else
        % Append to the existing TIFF file for subsequent images
        imwrite(imgplot, filename, 'tif', 'WriteMode', 'append', 'Compression', 'none');
    end
    
    % Close the figure to avoid memory issues
    close(h);
    ii
end
display('All traces saved!');
figure('visible', 'on');
close all
end




