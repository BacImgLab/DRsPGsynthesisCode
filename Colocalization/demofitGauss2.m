function [pk1, pk2, p_fit, Fitres] = demofitGauss2(X, intensity)
%
% This function fits a given 1D intensity profile (typically from a demograph)
% with a sum of two Gaussian functions. The two components are intended to
% capture the intensity distribution on the left and right halves of the profile.
%
% INPUTS:
%   X         - Vector of x-axis values (e.g., pixel positions)
%   intensity - Vector of intensity values corresponding to X; may include NaNs
%
% OUTPUTS:
%   pk1       - Center position (mu1) of the first Gaussian (typically the left peak)
%   pk2       - Center position (mu2) of the second Gaussian (typically the right peak)
%   p_fit     - Fitted parameters of the two Gaussians:
%               [A1, mu1, sigma1, A2, mu2, sigma2, offset]
%   Fitres    - Two-column matrix [X_real, Y_fit] containing the x-values used and
%               the corresponding fitted curve values
%
% Author: Xinxing Yang
% Date:   2025-04-05

    % Determine profile length and midpoint
    len = length(X);
    Index_mid = floor(len / 2);

    % Remove NaNs from the intensity data
    X_real = X(~isnan(intensity));
    I_real = intensity(~isnan(intensity));

    % Initialize parameters based on peak intensity
    [~, position] = max(I_real);
    locs = X_real(position);

    if locs > Index_mid
        mu2 = locs;
        mu1 = len - locs;
    else
        mu1 = locs;
        mu2 = len - locs;
    end

    % Estimate amplitudes at initial center guesses
    A1 = I_real(X_real == mu1);
    A2 = I_real(X_real == mu2);

    % If peak not found exactly at mu1/mu2, fall back to edge values
    if isempty(A1)
        A1 = X_real(1);
    end
    if isempty(A2)
        A2 = X_real(end);
    end

    % Initial guesses for standard deviation and offset
    sigma1 = 8;
    sigma2 = 8;
    offset = min(I_real);

    % Initial parameter vector
    % [A1, mu1, sigma1, A2, mu2, sigma2, offset]
    p0 = [A1(1), mu1, sigma1, A2(1), mu2, sigma2, offset];

    % Define the double Gaussian model
    doubleGauss = @(p, xdata) p(7) ...
        + p(1) * exp(-((xdata - p(2)).^2) / (2 * p(3)^2)) ...
        + p(4) * exp(-((xdata - p(5)).^2) / (2 * p(6)^2));

    % Optimization settings
    opts = optimoptions('lsqcurvefit', 'Display', 'off');

    % Lower and upper bounds for fitting parameters
    lb = [0, X_real(1) - 2, 0, 0, Index_mid, 0, -1];
    ub = [2, Index_mid, 15, 2, X_real(end) + 2, 15, 1];

    % Perform the curve fitting
    p_fit = lsqcurvefit(doubleGauss, p0, X_real, I_real, lb, ub, opts);

    % Extract fitted peak positions
    pk1 = p_fit(2);
    pk2 = p_fit(5);

    % Generate fitted curve
    Y = doubleGauss(p_fit, X_real);
    Fitres = [X_real, Y];
end
