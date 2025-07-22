function [R_V, R_0, Traj_V, Traj_0] = addVoneFLtrajs(Traj_struc, R_struc, Frame_L, frameL, TimeMatrix, V)
% Add directional velocity to a set of stationary trajectories and recalculate noise ratio.
%
% This function selects a group of simulated trajectories with 0 velocity from the output of `rcdfCal`
% that have a specified frame length (`frameL`). It then adds a constant velocity `V` to these traces 
% to simulate directional motion. The function recalculates the noise-to-displacement ratio (R) for 
% the new velocity-added traces and compares them to the original ones.
%
% Inputs:
%   Traj_struc : cell array   - Trajectories from `rcdfCal` (indexed by frame length and velocity)
%   R_struc    : cell array   - Corresponding R values from `rcdfCal`
%   Frame_L    : matrix       - Matrix of frame lengths
%   frameL     : scalar       - Desired frame length for selecting a subset
%   TimeMatrix : cell array   - Time arrays corresponding to each frame length
%   V          : scalar       - Velocity to be added (in same units as trajectories)
%
% Outputs:
%   R_V     : [1 x N]         - R values for new velocity-added trajectories
%   R_0     : [1 x N]         - Original R values from simulation with V = 0
%   Traj_V  : [N x T]         - Velocity-added trajectories
%   Traj_0  : [N x T]         - Original zero-velocity trajectories
%
% Author: Xinxing Yang

    % Step 1: Locate the index for requested frame length
    F_idx = find(Frame_L == frameL);  
    if isempty(F_idx)
        % Use the closest smaller value if exact match not found
        F_idX_temp = find(Frame_L <= frameL);
        F_idx = F_idX_temp(end);
        warning(['No exact match for frameL = ' num2str(frameL) ...
                 '. Using closest lower match: Frame_L = ' num2str(Frame_L(F_idx)) ' instead.']);
    end

    % Step 2: Extract original (V = 0) trajectories and R values
    Traj_0 = Traj_struc{F_idx, 1};  % Original traces (rows: traces, columns: time points)
    R_0    = R_struc{F_idx, 1};     % Corresponding R values
    TimeX  = TimeMatrix{F_idx};     % Time vector

    % Step 3: Add constant directional velocity to each trajectory
    V_dir  = TimeX * V;             % Create motion vector [T x 1]
    Traj_V = Traj_0 + V_dir;        % Add V to every trajectory (broadcasting)

    % Step 4: Recompute noise-to-signal ratio R for each velocity-added trace
    R_V = zeros(1, size(Traj_V, 1));
    for jj = 1:size(Traj_V, 1)
        Traj_temp = Traj_V(jj, :);  % One velocity-added trajectory
        [~, ~, ~, ~, Ratio] = linfitR_forcdf(TimeX, Traj_temp);
        R_V(1, jj) = abs(Ratio);   % Store absolute R value
    end
end

% Helper function
function [FitTrace, p, Displ, StD, Ratio] = linfitR_forcdf(Time, Trace)
% linfitR_forcdf - Fit a 1D trajectory to a linear model and calculate the noise ratio.
    p         = polyfit(Time, Trace, 1);                % Linear fit (velocity and offset)
    FitTrace  = polyval(p, Time);                       % Reconstruct fitted trace
    Displ     = p(1) * (Time(end) - Time(1));           % Displacement = v * Δt
    StD       = std(Trace - FitTrace);                  % Std deviation from fit
    Ratio     = StD / Displ;                            % Noise-to-signal ratio
end
