function [V, Displ, StD, Ratio, N] = tracedropoutXY(Time, TraceX, TraceY, Nboot, pdrop)
% Estimate motion parameters by randomly dropping points in 2D trajectories.
%
% This function performs Nboot rounds of random subsampling on a 2D (X-Y) trajectory
% by dropping out a specified proportion (pdrop) of data points. Each remaining subset 
% is fit to a linear model separately in X and Y. The output includes velocity, displacement, 
% noise level, and signal-to-noise ratio for each bootstrapped subset.
%
% Inputs:
%   Time    : vector (Nx1)  - Time points
%   TraceX  : vector (Nx1)  - X position trace
%   TraceY  : vector (Nx1)  - Y position trace
%   Nboot   : scalar        - Number of bootstrap iterations
%   pdrop   : scalar [0–1]  - Proportion of points to randomly drop (e.g. 0.3 = drop 30%)
%
% Outputs:
%   V       : [Nboot x 1]   - Estimated velocity magnitude from each iteration
%   Displ   : [Nboot x 1]   - Total displacement from linear fit in 2D
%   StD     : [Nboot x 1]   - Residual standard deviation (combined X and Y)
%   Ratio   : [Nboot x 1]   - Noise-to-signal ratio (StD / Displ)
%   N       : [Nboot x 1]   - Number of data points retained in each iteration
%
% Notes:
%   - This function is useful to assess robustness of velocity estimation
%     under missing data scenarios (e.g., due to dropout, bleaching).
%   - The underlying 1D fitting function used is `linfitR`.

    % Ensure column vectors
    Time   = reshape(Time,  [length(Time), 1]);
    TraceX = reshape(TraceX,[length(TraceX), 1]);
    TraceY = reshape(TraceY,[length(TraceY), 1]);

    % Generate dropout masks for each bootstrap iteration
    RandN = rand(length(Time), Nboot);    % Uniform random matrix
    LeftN = RandN > pdrop;                % Logical matrix: points to retain

    % Preallocate output
    V      = zeros(Nboot, 1);
    Displ  = zeros(Nboot, 1);
    StD    = zeros(Nboot, 1);
    Ratio  = zeros(Nboot, 1);
    N      = zeros(Nboot, 1);

    % Bootstrap iterations
    for ii = 1:Nboot
        idxKeep = find(LeftN(:,ii));  % Indices of retained time points
        TimeTemp    = Time(idxKeep);
        TraceXTemp  = TraceX(idxKeep);
        TraceYTemp  = TraceY(idxKeep);
        N(ii,1)     = length(TimeTemp);  % Save number of points used

        % Linear fit in X and Y separately
        [~, pnX, DisplX, StDX, ~] = linfitR(TimeTemp, TraceXTemp);
        [~, pnY, DisplY, StDY, ~] = linfitR(TimeTemp, TraceYTemp);

        % Combine vector quantities into 2D magnitude
        V(ii,1)      = sqrt(pnX(1)^2 + pnY(1)^2);           % Total velocity
        Displ(ii,1)  = sqrt(DisplX^2 + DisplY^2);           % 2D displacement
        StD(ii,1)    = sqrt(StDX^2 + StDY^2);               % 2D noise magnitude
        Ratio(ii,1)  = StD(ii,1) / Displ(ii,1);             % Noise-to-signal ratio
    end
end
