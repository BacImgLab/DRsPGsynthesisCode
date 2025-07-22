function Prob = getProbR(R_sample, Bin, epsl)
% Estimate the probability that a random sample falls within a specific interval.
%
% This function calculates the probability that the values in `R_sample` fall within
% a defined interval specified by `Bin = [lowerBound, upperBound]`. A small optional 
% epsilon value (`epsl`) is added to the result to avoid returning exact zero probability 
% in sparse data cases (e.g., when used in simulations or log-scale computations).
%
% Inputs:
%   R_sample : [N x 1] or [1 x N] vector - Sample values (from simulation or experiment)
%   Bin      : [1 x 2] vector            - Range to compute the probability [min, max]
%   epsl     : scalar (optional)         - Small offset to avoid zero probability (default: 1e-5)
%
% Output:
%   Prob     : scalar                    - Estimated probability (∈ [0, 1])
%
% Example:
%   R_sample = randn(1000,1);            % Generate normal-distributed samples
%   Prob = getProbR(R_sample, [-1,1]);   % Probability of falling within [-1,1]
%

    % Assign default epsilon if not provided
    if nargin < 3
        epsl = 1e-5;
    end

    % Total number of samples
    N_tot = length(R_sample);

    % Count how many samples fall within the interval [min(Bin), max(Bin)]
    N_p = sum(R_sample >= min(Bin) & R_sample <= max(Bin));

    % Compute probability and add epsilon
    Prob = N_p / N_tot + epsl;

    % Ensure result does not exceed 1
    if Prob > 1
        Prob = 1;
    end
end

