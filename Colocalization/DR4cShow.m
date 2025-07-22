function imgFlag = DR4cShow(ImShowing)
% DR4cShow - Display and evaluate a 4-channel image (e.g., DR segmented cell)
%
% This function visualizes a 4-channel image by displaying each channel in a 2x2 subplot
% with adjusted contrast. After displaying, it prompts the user to evaluate 
% whether the displayed cell is valid or not.
%
% INPUT:
%   ImShowing - 3D image stack with 4 slices (channels)
% OUTPUT:
%   imgFlag   - Binary flag:
%                 1 = Cell is good (user clicked "Yes")
%                 0 = Cell is not good (user clicked "No")


% --- Input validation ---
if ndims(ImShowing) ~= 3
    error('Input must be a 3D image stack.');
end

if size(ImShowing, 3) < 4
    error('Input image must contain at least 4 slices (i.e., 4 channels).');
end

% --- Contrast limits for percentile stretching ---
LowP = 5;   % Lower percentile for contrast scaling
HighP = 99; % Upper percentile

% Default output flag (cell is assumed good)
imgFlag = 1;

% --- Create a figure for visualization ---
hFig = figure( ...
    'Name', 'DR cell checking...', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'Position', [900, 150, 1500, 1200], ...
    'WindowStyle', 'normal');

% --- Display each channel in a 2x2 subplot layout ---
for ch = 1:4
    img = ImShowing(:,:,ch);
    
    % Compute percentile-based contrast range, excluding zero background
    pixelVals = double(img(img > 0));
    if isempty(pixelVals)
        minVal = 0;
        maxVal = 1;
    else
        minVal = prctile(pixelVals, LowP);
        maxVal = prctile(pixelVals, HighP);
    end
    
    % Display the image with adjusted contrast
    subplot(2, 2, ch);
    imshow(img, [minVal, maxVal]);
    title(sprintf('Channel %d', ch));
end

% Bring figure to front (in case it's behind other windows)
figure(hFig);

% --- Prompt user to evaluate the cell ---
answer = questdlg('Is the current cell good?', ...
                  'Cell Evaluation', ...
                  'Yes', 'No', 'Yes');

% Set output flag based on user response
if strcmp(answer, 'No')
    imgFlag = 0;
else
    imgFlag = 1;
end

% --- Close the figure ---
close(hFig);

end
