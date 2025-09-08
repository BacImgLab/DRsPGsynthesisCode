function DrunwrapX()
% DrunwrapX
% -------------------------------------------------------------------------
% Purpose:
%   - Open an interactive GUI for browsing and processing trajectory data
%   - Display trajectories in BF (Bright Field) and FI (Fluorescence Image)
%   - Allow the user to select two points as references for unwrapping
%   - Save selected points and compute unwrapped trajectory coordinates
%
% Input:
%   - Driftout-TrackRefine-with-Coordsnew.mat
%       Must contain the structure array "TrackRefine"
%
% Output:
%   - selected_points.mat
%       Stores user-selected reference points for each trajectory
%   - TrackRefine_unwrapped.mat
%       Contains TrackRefine with additional unwrapped coordinates
% -------------------------------------------------------------------------

    % -------- Check input file --------
    if ~isfile('Driftout-TrackRefine-with-Coordsnew.mat')
        error('File Driftout-TrackRefine-with-Coordsnew.mat not found.');
    end
    s = load('Driftout-TrackRefine-with-Coordsnew.mat');
    if ~isfield(s, 'TrackRefine')
        error('TrackRefine structure not found.');
    end
    data = s.TrackRefine;

    % -------- Initialization --------
    numTraces = numel(data);
    traceIndex = 1;
    selectedPoints = [];
    allSelectedPoints = cell(numTraces, 1);

    % -------- GUI setup --------
    fig = figure('Name', 'Trace Viewer', 'Position', [100, 100, 1800, 1000]);

    axBF = subplot(1,2,1);
    axFI = subplot(1,2,2);
    title(axBF, 'BF Image');
    title(axFI, 'FI Image');

    buttonY = 40;
    buttonW = 100;
    buttonH = 40;
    spacing = 20;
    centerX = 900;

    % Control buttons
    uicontrol('Style', 'pushbutton', 'String', 'Previous', ...
        'Position', [centerX - 4*buttonW - 3.5*spacing, buttonY, buttonW, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) changeTrace(-1));

    uicontrol('Style', 'pushbutton', 'String', 'Next', ...
        'Position', [centerX - 2*buttonW - spacing, buttonY, buttonW, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) changeTrace(1));

    uicontrol('Style', 'pushbutton', 'String', 'Clear Selection', ...
        'Position', [centerX - buttonW/2, buttonY, buttonW+10, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) cancelSelection());

    uicontrol('Style', 'pushbutton', 'String', 'Confirm Selection', ...
        'Position', [centerX + buttonW + 1.5*spacing, buttonY, buttonW+10, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) confirmSelection());

    uicontrol('Style', 'pushbutton', 'String', 'Save All Points', ...
        'Position', [centerX + 2*buttonW + 3*spacing, buttonY, buttonW+40, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) saveSelections());

    uicontrol('Style', 'pushbutton', 'String', 'Compute Unwrap', ...
        'Position', [centerX + 3*buttonW + 5*spacing + 40, buttonY, buttonW+40, buttonH], ...
        'FontSize', 12, 'Callback', @(~,~) computeUnwrap());

    % Mouse click handler
    set(fig, 'WindowButtonDownFcn', @onMouseClick);

    % Show first valid trace
    showTrace();

    % =====================================================================
    % ------------------------- Nested Functions --------------------------
    % =====================================================================

    function showTrace()
        cla(axBF); cla(axFI);
        selectedPoints = [];

        % Skip traces that do not meet conditions
        while ~shouldDisplayTrace(data{traceIndex})
            traceIndex = traceIndex + 1;
            if traceIndex > numTraces
                traceIndex = 1;
            end
        end

        T = data{traceIndex};
        BF = T.ROIBFrot;
        FI = T.ROIFLrot;
        Trace = T.Coord;
        ROIRotate = T.ROIRotate;
        ROIRecRotate = T.ROIRecRotate;

        % Crop BF and FI images around ROI
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

        % Show BF image
        axes(axBF);
        imshow(BF_crop, [], 'XData', xlim_val, 'YData', ylim_val);
        hold on;
        plot(ROIRotate(:,1), ROIRotate(:,2), ':r', 'LineWidth', 2);
        plot(ROIRecRotate(:,1), ROIRecRotate(:,2), '-b', 'LineWidth', 1.5);
        plot(Trace(:,6), Trace(:,7), '-y', 'LineWidth', 1.5);
        set(axBF, 'YDir', 'normal');
        hold off;

        % Show FI image
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

    % (keep other nested functions in same style...)
end

% =====================================================================
% ------------------------ Helper Functions ---------------------------
% =====================================================================

function [Coord_unwrap] = unwrapTraj(Center, R, Coord_2d)
    % UNWRAPTRAJ
    % -----------------------------------------------------------------
    % Convert 2D coordinates into unwrapped coordinates based on
    % reference center and radius.
    %
    % Input:
    %   Center   - [Xc, Yc], reference point
    %   R        - radius for unwrapping
    %   Coord_2d - [X, Y], original coordinate
    %
    % Output:
    %   Coord_unwrap - [X', Y'], unwrapped coordinate
    % -----------------------------------------------------------------

    Xc = Center(1);
    Yc = Center(2);
    Coord = zeros(1,2);

    % Translate Y coordinate
    Coord(2) = Coord_2d(2) - Yc;

    % X coordinate mapped via asin
    DisX = abs(Coord_2d(1) - Xc);
    Xnew = R * asin(DisX / R);

    if Coord_2d(1) - Xc >= 0
        Coord(1) = Xnew;
    else
        Coord(1) = -Xnew;
    end

    Coord_unwrap = Coord;
end
