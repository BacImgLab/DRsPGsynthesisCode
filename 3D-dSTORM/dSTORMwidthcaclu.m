%% dSTORMwidthcaclu.m
% This script batch processes CSV curve data in a specified folder,
% calculates the Full Width at Half Maximum (FWHM) for each curve,
% saves plots for each curve with FWHM annotation, 
% and outputs a summary CSV with all calculated widths.
%
clear;
clc;

% ----------- Set folder paths -----------
sourceFolder = 'C:\Users\WT_S0_width_compared';  % <-- Replace with your actual path
outputFolder = fullfile(sourceFolder, 'Plots');

% Create output folder if it does not exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Get list of all CSV files in the source folder
fileList = dir(fullfile(sourceFolder, '*.csv'));

% Initialize cell array to store all FWHM results (with header)
widthResults = {'FileName', 'HalfPeakWidth'};

% ----------- Process each CSV file -----------
for i = 1:length(fileList)
    fileName = fileList(i).name;
    filePath = fullfile(sourceFolder, fileName);

    % ----- Read data -----
    rawData = readmatrix(filePath);
    
    % Remove first row if it contains NaNs (e.g., header or invalid row)
    if isnan(rawData(1,1)) || isnan(rawData(1,2))
        rawData = rawData(2:end, :);
    end

    x = rawData(:,1);  % X-axis data (e.g., distance)
    y = rawData(:,2);  % Y-axis data (e.g., gray value)

    % ----- Calculate Full Width at Half Maximum (FWHM) -----
    y_max = max(y);
    y_min = min(y);
    halfHeight = (y_max + y_min) / 2;  % Half peak height

    % Find all points where curve crosses the half-height using linear interpolation
    X_cross = [];
    for j = 2:length(y)
        % Check if y crosses halfHeight between points j-1 and j
        if (y(j-1) - halfHeight) * (y(j) - halfHeight) < 0
            % Linear interpolation to find exact crossing position
            ratio = (halfHeight - y(j-1)) / (y(j) - y(j-1));
            x_interp = x(j-1) + ratio * (x(j) - x(j-1));
            X_cross(end+1) = x_interp;
        elseif y(j) == halfHeight
            % If y exactly equals halfHeight at point j, record that point
            X_cross(end+1) = x(j);
        end
    end

    % Calculate width as distance between first and last crossing points if >=2 crossings found
    if length(X_cross) >= 2
        X1 = min(X_cross);
        X2 = max(X_cross);
        width = abs(X2 - X1);
    else
        % Insufficient crossings to determine FWHM
        X1 = NaN;
        X2 = NaN;
        width = NaN;
    end

    % ----- Save results -----
    widthResults{end+1, 1} = fileName;
    widthResults{end, 2} = width;

    % ----- Visualization -----
    figure('Visible','off');
    plot(x, y, 'b-', 'LineWidth', 1.5); hold on;

    % Plot half-height line
    yline(halfHeight, 'k--', 'LineWidth', 1.2, 'DisplayName','Half Height');

    % Mark half-height crossing points and annotate width
    if ~isnan(X1) && ~isnan(X2)
        plot([X1 X2], [halfHeight halfHeight], 'r--', 'LineWidth', 1.5, 'DisplayName','FWHM line');
        plot(X1, halfHeight, 'ro', 'MarkerFaceColor','r');
        plot(X2, halfHeight, 'ro', 'MarkerFaceColor','r');
        text((X1+X2)/2, halfHeight, sprintf('Width = %.2f', width), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'r');
    end

    xlabel('Distance');
    ylabel('Gray Value');
    title(sprintf('%s - FWHM Analysis', fileName), 'Interpreter', 'none');
    legend('show');
    grid on;

    % Save figure as PNG
    [~, name, ~] = fileparts(fileName);
    saveas(gcf, fullfile(outputFolder, [name, '_HalfPeakWidth.png']));
    close(gcf);
end

% ----------- Save all FWHM results to CSV -----------
outputCSV = fullfile(outputFolder, 'AllWidths.csv');
writecell(widthResults, outputCSV);

disp('✅ All curves processed. Plots and FWHM data saved.');
