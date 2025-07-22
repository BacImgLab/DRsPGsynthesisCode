function [R_struc, Traj_struc, TimeMatric, frameMatrix, SpeedMatrix] = rcdfCal(N_frame_list, V_x_list, T_exp, D_eff, Boundary, Local_err, N_traj)
% rcdfCal - Run multiple trajectory simulations with drift and diffusion,
% then compute the ratio (StD/Displacement) for each simulated trace.
%
% Inputs:
%   N_frame_list : array  - List of frame numbers to simulate (e.g. [10 20 30])
%   V_x_list     : array  - List of drift velocities (e.g. [0.01 0.02 0.05] µm/s)
%   T_exp        : scalar - Exposure time per frame (in seconds)
%   D_eff        : scalar - Diffusion coefficient (in µm²/s)
%   Boundary     : scalar - Confinement boundary size (in nm)
%   Local_err    : scalar - Localization error (in nm)
%   N_traj       : scalar - Number of simulated trajectories per condition
%
% Outputs:
%   R_struc      : cell   - Noise-to-signal ratio for each simulated trace (cell array indexed by [i,j])
%   Traj_struc   : cell   - Simulated trace data
%   TimeMatric   : cell   - Time vector corresponding to each frame length
%   frameMatrix  : matrix - Matrix of frame counts used for each condition
%   SpeedMatrix  : matrix - Matrix of velocities used for each condition
%
% Description:
%   For each combination of frame length and velocity, this function simulates N_traj
%   1D confined drift-diffusion trajectories using `confdirTrackGene`, fits them
%   linearly, and computes the ratio of standard deviation of residuals to total displacement.
%   This ratio is used as a measure of signal quality under varying simulation parameters.

    for ii = 1:length(N_frame_list)  % Iterate over frame lengths
        tic
        for jj = 1:length(V_x_list)  % Iterate over drift speeds
            N_frame = N_frame_list(ii);
            V_x = V_x_list(jj);

            R = zeros(1, N_traj);          % Preallocate ratio result
            Time = (1:N_frame)' * T_exp;   % Time vector
            Traj_temp = zeros(N_traj, N_frame);  % Store traces

            parfor kk = 1:N_traj  % Parallel loop over trajectories
                % Generate trajectory with drift + diffusion
                Traj = confdirTrackGene(N_frame, T_exp, D_eff, V_x, Boundary, Local_err);

                % Perform linear fit and calculate ratio
                [~, ~, ~, ~, Ratio] = linfitR_forcdf(Time, Traj);

                % Store absolute ratio
                R(kk) = abs(Ratio);
                Traj_temp(kk, :) = Traj;
            end

            % Store results
            R_struc{ii, jj} = R;
            Traj_struc{ii, jj} = Traj_temp;
            frameMatrix(ii, jj) = N_frame;
            SpeedMatrix(ii, jj) = V_x;

            disp(['Simulating: frame = ' num2str(ii) ' of ' num2str(length(N_frame_list)) ...
                ', V = ' num2str(jj) ' of ' num2str(length(V_x_list)) '...']);
        end

        % Store time axis for current frame length
        TimeMatric{ii, 1} = (1:N_frame) * T_exp;
        toc
    end
end

% -------------------------------------------------------------------------
function [FitTrace, p, Displ, StD, Ratio] = linfitR_forcdf(Time, Trace)
% linfitR_forcdf - Fit a trace linearly and compute noise-to-signal ratio.
%
% Inputs:
%   Time  : column vector - Time points
%   Trace : column vector - Simulated position at each time point
%
% Outputs:
%   FitTrace : fitted trace from linear regression
%   p        : fit parameters (p(1) = velocity, p(2) = intercept)
%   Displ    : total displacement based on the linear fit
%   StD      : standard deviation of the residuals
%   Ratio    : StD / Displacement (noise-to-signal)

    % Fit trace with linear polynomial (velocity + offset)
    p = polyfit(Time, Trace, 1);

    % Predict fitted trace from model
    FitTrace = polyval(p, Time);

    % Compute displacement (based on fitted velocity)
    Displ = p(1) * (Time(end) - Time(1));

    % Residual standard deviation (i.e., noise level)
    StD = std(Trace - FitTrace);

    % Noise-to-signal ratio
    Ratio = StD / Displ;
end
