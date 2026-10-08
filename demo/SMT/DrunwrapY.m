function DrunwrapY()
% DrunwrapY
% -------------------------------------------------------------------------
% This interactive MATLAB GUI tool allows users to:
%   1. Load trajectories from `Driftout-TrackRefine-with-Coordsnew.mat`.
%   2. Visualize Bright Field (BF) and Fluorescence (FI) images side by side.
%   3. Navigate through traces and manually select two reference points.
%   4. Save all selected points to `selected_points.mat`.
%   5. Compute unwrapped coordinates for each trajectory based on the
%      manually chosen center and radius, and save results to
%      `TrackRefine_unwrapped.mat`.
%
% The unwrap procedure is based on a circular-to-linear coordinate
% transformation. The selected two points define the center and radius (R)
% used to transform each trajectory segment into an "unwrapped" coordinate
% system.
%
% Requirements:
%   - Input file: Driftout-TrackRefine-with-Coordsnew.mat
%   - Data structure inside the file must include a variable named
%     `TrackRefine`, which contains trajectories with fields:
%         - ROIBFrot, ROIFLrot: cropped BF/FI images
%         - Coord: raw trajectory coordinates
%         - ROIRotate, ROIRecRotate: ROI shapes
%         - Seg: segmented trajectories (with Coordsnew field)
% -------------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % Load data and perform initial checks
    % ---------------------------------------------------------------------
    if ~isfile('Driftout-TrackRefine-with-Coordsnew.mat')
        error('File Driftout-TrackRefine-with-Coordsnew.mat not found.');
    end
    s = load('Driftout-TrackRefine-with-Coordsnew.mat');
    if ~isfield(s, 'TrackRefine')
        error('TrackRefine field is missing in the loaded file.');
    end
    data = s.TrackRefine;
    numTraces = numel(data);
    traceIndex = 1;                   % Current trace index
    selectedPoints = [];              % Points selected on the current trace
    allSelectedPoints = cell(numTraces, 1); % Store selections for all traces

    % ---------------------------------------------------------------------
    % Initialize GUI
    % ---------------------------------------------------------------------
    fig = figure('Name', 'Trace Viewer', 'Position', [100, 100, 1800, 1000]);
    axBF = subplot(1,2,1);
    axFI = subplot(1,2,2);
    title(axBF, 'BF Image');
    title(axFI, 'FI Image');

    % Button layout parameters
    buttonY = 40;
    buttonW = 100;
    buttonH = 40;
    spacing = 20;
    centerX = 900;

    % Navigation and control buttons
    uicontrol('Style', 'pushbutton', 'String', 'Previous', ...
        'Position', [centerX - 4*buttonW - 3.5*spacing, buttonY, buttonW, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) changeTrace(-1));

    uicontrol('Style', 'pushbutton', 'String', 'Next', ...
        'Position', [centerX - 2*buttonW - spacing, buttonY, buttonW, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) changeTrace(1));

    uicontrol('Style', 'pushbutton', 'String', 'Cancel Selection', ...
        'Position', [centerX - buttonW/2, buttonY, buttonW+10, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) cancelSelection());

    uicontrol('Style', 'pushbutton', 'String', 'Confirm Selection', ...
        'Position', [centerX + buttonW + 1.5*spacing, buttonY, buttonW+10, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) confirmSelection());

    uicontrol('Style', 'pushbutton', 'String', 'Save Points', ...
        'Position', [centerX + 2*buttonW + 3*spacing, buttonY, buttonW+40, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) saveSelections());

    uicontrol('Style', 'pushbutton', 'String', 'Compute Unwrap', ...
        'Position', [centerX + 3*buttonW + 5*spacing + 40, buttonY, buttonW+40, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) computeUnwrap());

    % Register mouse-click callback for point selection
    set(fig, 'WindowButtonDownFcn', @onMouseClick);

    % Display first trace
    showTrace();

    % ---------------------- Internal Functions ---------------------------
    % Each of the following functions handles specific tasks in the GUI
    % such as trace visualization, navigation, point selection, and
    % coordinate transformation.
    % ---------------------------------------------------------------------

    function showTrace()
        % Clear both axes
        cla(axBF); cla(axFI);
        selectedPoints = [];

        % Skip traces that do not meet display criteria
        while ~shouldDisplayTrace(data{traceIndex})
            traceIndex = traceIndex + 1;
            if traceIndex > numTraces
                traceIndex = 1;
            end
        end

        % Extract data for current trace
        T = data{traceIndex};
        BF = T.ROIBFrot;
        FI = T.ROIFLrot;
        Trace = T.Coord;
        ROIRotate = T.ROIRotate;
        ROIRecRotate = T.ROIRecRotate;

        % Define cropping region based on ROI
        BFWidth = range(ROIRecRotate(:,1)) * 1.4;
        BFHeight = range(ROIRecRotate(:,2)) * 2;
        [H, W, ~] = size(BF);
        y1 = max(1, round((H - BFHeight)/2));
        y2 = min(H, y1 + round(BFHeight));
        x1 = max(1, round((W - BFWidth)/2));
        x2 = min(W, x1 + round(BFWidth));

        BF_crop = BF(y1:y2, x1:x2);
        FI_crop = FI(y1:y2, x1:x2);
        xlim_val = [-BFWidth/2, BFWidth/2];
        ylim_val = [-BFHeight/2, BFHeight/2];

        % Display BF image with ROI and trajectory overlay
        axes(axBF);
        imshow(BF_crop, [], 'XData', xlim_val, 'YData', ylim_val);
        hold on;
        plot(ROIRotate(:,1), ROIRotate(:,2), ':r', 'LineWidth', 2);
        plot(ROIRecRotate(:,1), ROIRecRotate(:,2), '-b', 'LineWidth', 1.5);
        plot(Trace(:,6), Trace(:,7), '-y', 'LineWidth', 1.5);
        set(axBF, 'YDir', 'normal');
        hold off;

        % Display FI image with overlays
        axes(axFI);
        imshow(FI_crop, [], 'XData', xlim_val, 'YData', ylim_val);
        hold on;
        plot(ROIRotate(:,1), ROIRotate(:,2), ':r', 'LineWidth', 2);
        plot(ROIRecRotate(:,1), ROIRecRotate(:,2), '-y', 'LineWidth', 1.5);
        plot(Trace(:,6), Trace(:,7), '-c', 'LineWidth', 1.5);
        set(axFI, 'YDir', 'normal');
        title(axFI, sprintf('FI Image - Trace %d / %d', traceIndex, numTraces));
        hold off;

        % Restore previously selected points if available
        if ~isempty(allSelectedPoints{traceIndex})
            selectedPoints = allSelectedPoints{traceIndex};
        else
            selectedPoints = [];
        end
        plotSelectedLine();
    end

    function valid = shouldDisplayTrace(traceData)
        % Only display traces that contain at least one segment
        % with Location == 3
        valid = false;
        if ~isfield(traceData, 'Seg') || isempty(traceData.Seg)
            return;
        end
        for k = 1:numel(traceData.Seg)
            if isfield(traceData.Seg(k), 'Location') && traceData.Seg(k).Location == 3
                valid = true;
                return;
            end
        end
    end

    function changeTrace(step)
        % Navigate to next/previous trace
        traceIndex = traceIndex + step;
        if traceIndex > numTraces
            traceIndex = 1;
        elseif traceIndex < 1
            traceIndex = numTraces;
        end
        showTrace();
    end

    function onMouseClick(~,~)
        % Select up to two points on the BF image
        ax = gca;
        if ax ~= axBF
            return;
        end
        pt = get(axBF, 'CurrentPoint');
        pt = pt(1,1:2);
        if size(selectedPoints,1) < 2
            selectedPoints(end+1,:) = pt;
            plotSelectedLine();
        else
            msgbox('Only two points are allowed. Use "Cancel Selection" to reselect.');
        end
    end

    function plotSelectedLine()
        % Display selected points and connecting line
        axes(axBF);
        hold on;
        delete(findobj(axBF, 'Type', 'scatter'));
        delete(findobj(axBF, 'Type', 'line', '-and', 'LineStyle', '--'));
        if size(selectedPoints,1) >= 1
            scatter(selectedPoints(:,1), selectedPoints(:,2), 100, 'r', 'filled');
        end
        if size(selectedPoints,1) == 2
            plot(selectedPoints(:,1), selectedPoints(:,2), '--w', 'LineWidth', 2);
        end
        hold off;
    end

    function cancelSelection()
        % Clear current selection
        selectedPoints = [];
        allSelectedPoints{traceIndex} = [];
        plotSelectedLine();
    end

    function confirmSelection()
        % Confirm selection of two points for the current trace
        if size(selectedPoints,1) ~= 2
            msgbox('Please select exactly two points first.');
            return;
        end
        allSelectedPoints{traceIndex} = selectedPoints;
        msgbox(sprintf('Points confirmed for trace %d.', traceIndex));
    end

    function saveSelections()
        % Save all confirmed selections to file
        result = struct('TraceIndex', {}, 'Point1', {}, 'Point2', {});
        for i = 1:numTraces
            if ~isempty(allSelectedPoints{i})
                result(end+1).TraceIndex = i;
                result(end).Point1 = allSelectedPoints{i}(1,:);
                result(end).Point2 = allSelectedPoints{i}(2,:);
            end
        end
        save('selected_points.mat', 'result');
        msgbox('All confirmed points saved to selected_points.mat');
    end

    function computeUnwrap()
        % Compute unwrapped coordinates for each trace
        for i = 1:numTraces
            if isempty(allSelectedPoints{i})
                continue;
            end

            pt1 = allSelectedPoints{i}(1,:);
            pt2 = allSelectedPoints{i}(2,:);
            Center = pt1;
            R = abs(pt2(2) - pt1(2));  % Initial radius estimate
            maxDisY = 0;

            % Find maximum vertical distance among all segments
            if isfield(data{i}, 'Seg')
                for k = 1:numel(data{i}.Seg)
                    if ~isfield(data{i}.Seg(k), 'Coordsnew')
                        continue;
                    end
                    coords = data{i}.Seg(k).Coordsnew;
                    if isempty(coords) || size(coords,2) < 3
                        continue;
                    end
                    for j = 1:size(coords,1)
                        disY = abs(coords(j,3) - Center(2));
                        if disY > maxDisY
                            maxDisY = disY;
                        end
                    end
                end
            end

            % Update radius
            if maxDisY > R
                R = maxDisY;
            end

            % Save Center and R to data structure
            data{i}.UnwrapCenter = Center;
            data{i}.UnwrapR = R;

            % Perform unwrapping for each segment
            if isfield(data{i}, 'Seg')
                for k = 1:numel(data{i}.Seg)
                    if ~isfield(data{i}.Seg(k), 'Coordsnew')
                        continue;
                    end
                    coords = data{i}.Seg(k).Coordsnew;
                    if isempty(coords) || size(coords,2) < 3
                        continue;
                    end

                    coords_unwrap = zeros(size(coords));
                    coords_unwrap(:,1) = coords(:,1);  % Time column
                    for j = 1:size(coords,1)
                        x = coords(j,2);
                        y = coords(j,3);
                        unwrapped = unwrapTraj(Center, R, [x, y]);
                        coords_unwrap(j,2:3) = unwrapped;
                    end

                    data{i}.Seg(k).Coords_unwrap = coords_unwrap;
                    data{i}.Seg(k).UnwrapCenter = Center;
                    data{i}.Seg(k).UnwrapR = R;
                end
            end
        end

        % Save unwrapped data
        TrackRefine = data;
        save('TrackRefine_unwrapped.mat', 'TrackRefine');
        msgbox('Trajectory unwrapping completed. Saved as TrackRefine_unwrapped.mat');
    end

end

% -------------------------------------------------------------------------
% Helper function: unwrapTraj
% -------------------------------------------------------------------------
function [Coord_unwrap] = unwrapTraj(Center, R, Coord_2d)
% UNWRAPTRAJ Perform circular-to-linear coordinate transformation
%   Input:
%       Center   - [Xc, Yc], center point for unwrapping
%       R        - Radius (scaling factor based on Y-distance)
%       Coord_2d - [X, Y], original 2D coordinate
%   Output:
%       Coord_unwrap - [Xnew, Ynew], unwrapped coordinate
%
% Transformation rule:
%   - X coordinate is translated relative to Center.
%   - Y coordinate is transformed using arc length (R * asin(ΔY / R)).

    Xc = Center(1);
    Yc = Center(2);
    Coord = zeros(1,2);

    % Translate X coordinate
    Coord(1) = Coord_2d(1) - Xc;

    % Transform Y coordinate
    DisY = abs(Coord_2d(2) - Yc);
    Ynew = R * asin(DisY / R);
    if Coord_2d(2) - Yc >= 0
        Coord(2) = Ynew;
    else
        Coord(2) = -Ynew;
    end

    Coord_unwrap = Coord;
end
