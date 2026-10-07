function mergeSegData(parentDir)
% Merge segmented trajectory data from multiple subfolders.
%
% Inputs:
%   parentDir (optional) : string - Path to the root directory that contains
%                                     subfolders with 'TraceRefine-test.mat' files.
%
% Description:
%   This function recursively traverses all subdirectories under a given
%   parent directory and loads `TraceRefine-test.mat` files from each.
%   It extracts all valid segments (`Seg`) from the `TrackRefine` structure,
%   records them along with their original Trace IDs, and saves them to
%   a new file `mergedAllSeg.mat` for downstream analysis.
%
%   Useful when you have performed trajectory segmentation across
%   multiple conditions or batches and now want to consolidate the segments.
%
% Outputs (saved to file):
%   - mergedAllSeg : cell array - All individual trajectory segments
%   - SegmentID    : cell array - Corresponding Trace ID for each segment
%
% Notes:
%   - Will prompt the user to select a folder if `parentDir` is not provided.
%   - Issues warnings for folders missing expected files or fields.

clear;clc;

    % If no input argument, open folder selection UI
    if nargin == 0
        parentDir = uigetdir;  % Prompt user to select parent directory
    end

    % Validate folder selection
    if parentDir == 0
        error('No folder selected. Please select the root directory.');
    end

    % Initialize variables for collecting all segments and IDs
    mergedAllSeg = {};      % To store all valid segments
    SegmentID = {};         % To store corresponding Trace IDs
    Segcount = 1;           % Index for storing in mergedAllSeg

    % Get list of subfolders (excluding '.' and '..')
    subFolders = dir(parentDir);
    subFolders = subFolders([subFolders.isdir] & ~ismember({subFolders.name}, {'.', '..'}));

    % Loop through each subfolder
    for i = 1:numel(subFolders)
        subFolderPath = fullfile(parentDir, subFolders(i).name);  
        matFilePath = fullfile(subFolderPath, 'TraceRefine-test.mat');  % Expected file name

        if exist(matFilePath, 'file')
            data = load(matFilePath);  % Load .mat file

            if isfield(data, 'TrackRefine')  % Check expected structure
                trackRefine = data.TrackRefine;

                % Loop over each trace in TrackRefine
                for j = 1:numel(trackRefine)
                    currentStruct = trackRefine{j};              % One trajectory struct
                    currentTraceID = currentStruct.TraceId;      % Corresponding Trace ID

                    % If segment field exists and is not empty
                    if isfield(currentStruct, 'Seg') && ~isempty(currentStruct.Seg)
                        currentSeg = currentStruct.Seg;

                        % Loop through each segment and save
                        for k = 1:size(currentSeg, 2)
                            mergedAllSeg{Segcount, 1} = currentSeg(1, k);  % Store segment
                            SegmentID{Segcount, 1} = currentTraceID;       % Store Trace ID
                            Segcount = Segcount + 1;
                        end
                    end
                end
            else
                warning('The file %s does not contain the field ''TrackRefine''.', matFilePath);
            end
        else
            warning('The file %s does not exist.', matFilePath);
        end
    end

    % Save merged result to disk
    save('mergedAllSeg.mat', 'mergedAllSeg', 'SegmentID');
    disp('All segments and their trace IDs have been saved to mergedAllSeg.mat.');
end
