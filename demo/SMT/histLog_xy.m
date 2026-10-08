function [HisRes, HisX, HisObj] = histLog_xy(x, LowLog, HighLog, Nbin)
% Generate a normalized histogram of data `x` using logarithmic bins.
%
% This function creates a histogram for the input vector `x` with logarithmically spaced
% bins between 10^LowLog and 10^HighLog, using `Nbin` number of bins. The result is a
% normalized probability distribution (i.e., the sum of histogram bars equals 1).
% It returns the histogram values, bin edges, and histogram object.
%
% Inputs:
%   x        : [N x 1] or [1 x N] vector - The data to histogram (positive values only)
%   LowLog   : scalar                    - log10 of the lower boundary (e.g. -2 for 10^-2)
%   HighLog  : scalar                    - log10 of the upper boundary (e.g. 1 for 10^1)
%   Nbin     : integer                   - Number of logarithmic bins
%
% Outputs:
%   HisRes   : [1 x (Nbin-1)] vector     - Normalized histogram values (probabilities)
%   HisX     : [1 x Nbin] vector         - Bin edges in log space (10^LowLog to 10^HighLog)
%   HisObj   : histogram object          - Handle to the histogram object
%
% Notes:
% - Data in `x` must be positive to compute logarithmic histograms.
% - The y-axis is in probability scale (normalized to sum to 1).
% - The x-axis is automatically scaled to log10.
%
% Example:
%   x = rand(1,1000) * 100;
%   [HisRes, HisX, HisObj] = histLog_xy(x, 0, 2, 20);


    % Input validation
    if nargin < 4
        error('Not enough inputs. Required: x, LowLog, HighLog, Nbin.');
    end

    % Generate logarithmically spaced bin edges
    Xbin = logspace(LowLog, HighLog, Nbin);

    % Compute histogram bin edges (technically redundant, but used for clarity)
    [~, edges] = histcounts(x, Xbin);

    % Create histogram with defined edges and normalize to probability
    HisObj = histogram(x, edges, 'Normalization', 'probability');
    HisObj.FaceColor = 'w';  % Set bar color to white

    % Output the bin edges and corresponding histogram values
    HisX = HisObj.BinEdges;
    HisRes = HisObj.Values;

    % Set X-axis to logarithmic scale
    set(gca, 'XScale', 'log');
end

