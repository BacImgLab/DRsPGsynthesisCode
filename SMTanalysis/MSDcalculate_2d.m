function s_MSD = MSDcalculate_2d(struc, pixelsize, frameT, darkT, CMode)
% Calculate mean squared displacement (MSD) for 2D particle trajectories.
%
% This function computes the MSD based on trajectory structures from tracking experiments.
% It works in 2D (X, Y), using either all possible frame intervals or only consecutive frames.
%
% Inputs:
%   struc      : struct        - Input tracking structure with 'TracksROI' field
%   pixelsize  : scalar [nm]   - Pixel size in nanometers
%   frameT     : scalar [s]    - Frame exposure time in seconds
%   darkT      : scalar [s]    - Dark time (delay between frames)
%   CMode      : string        - Calculation mode (default = 'all')
%                               'all'    : use all frame intervals
%                               'consec' : use only continuous sub-trajectories
%
% Output:
%   s_MSD : matrix [Nx4]       - MSD result table:
%       Column 1: Time delay (lag time)
%       Column 2: Mean squared displacement (MSD) in nm²
%       Column 3: Bootstrapped standard deviation of MSD
%       Column 4: Standard error of the bootstrapped mean
%
% Notes:
%   - Relies on an external function `TrajConsec` if using 'consec' mode.
%   - MSD is averaged over all pairs of points with the same frame lag.
%

    % Handle default calculation mode
    if nargin < 5
        CMode = 'all';
    end
    if nargin < 4
        error('Not enough input arguments. Require: struc, pixelsize, frameT, darkT.');
        return;
    end

    % Preprocess input structure depending on mode
    if strcmp(CMode, 'consec')
        struc_r = TrajConsec(struc);  % Get only continuous sub-trajectories
    else
        struc_r = struc;
    end
    S = struc_r.TracksROI;

    % Determine maximum length of MSD to compute (half of max trajectory length)
    for idxN = 1:length(S)
        Coordinates = S(idxN).Coordinates;
        N_msd(idxN) = Coordinates(end,1) - Coordinates(1,1);
    end
    N_msdmax = floor(max(N_msd) / 2);

    % Initialize MSD result storage
    MSD = nan(N_msdmax, 4);

    % Main loop to calculate MSD for each lag frame
    for idxM = 1:N_msdmax
        MSD_temp = [];  % temporary storage for square displacements

        for idxN = 1:length(S)
            Coordinates = S(idxN).Coordinates;
            frames = Coordinates(:,1) - Coordinates(1,1) + 1;

            for idxQ = 1:min((frames(end) - idxM), length(frames))
                Startf = frames(idxQ);
                Endf = find(frames == (Startf + idxM));

                if ~isempty(Endf)
                    Startp = Coordinates(idxQ, 2:3);  % [X, Y] start
                    Endp = Coordinates(Endf, 2:3);    % [X, Y] end

                    % Compute squared displacement in nanometers
                    SquareD = sum(((Startp - Endp) * pixelsize).^2);
                    MSD_temp = [MSD_temp; SquareD];
                end
            end
        end

        % Store statistics
        if isempty(MSD_temp)
            MSD(idxM, 2) = NaN;
        else
            MSD(idxM, 2) = mean(MSD_temp);
            stats = bootstrp(1000, @mean, MSD_temp);  % Bootstrap std estimation
            MSD(idxM, 3) = std(stats);
            MSD(idxM, 4) = std(stats) / sqrt(length(stats) - 1);  % standard error
        end
    end

    % Construct time vector for each lag
    Timebin = (1:N_msdmax) * (frameT + darkT);
    MSD(:, 1) = Timebin';

    % Return the MSD matrix
    s_MSD = MSD;

end
