function colorCodeTracePlot(Time, Trace, linewidth)
% Plot a 2D trajectory colored by time progression.
%
% This function visualizes a 2D trajectory (X, Y) where the line color
% represents the time evolution of the trace. The color changes smoothly
% from the start to the end of the trace, using the 'jet' colormap.
%
% Inputs:
%   Time      : vector [Nx1] - Time points corresponding to each position
%   Trace     : matrix [Nx2] - 2D positions; column 1 is X, column 2 is Y
%   linewidth : scalar       - Thickness of the plotted trajectory line
%
% Note:
%   - Trace should have two columns representing X and Y positions.
%   - Time must be the same length as Trace rows.
%

    % Prepare data for surface plot with color corresponding to time
    xS = [Trace(:,1), Trace(:,1)];          % duplicate X for surface strip
    yS = [Trace(:,2), Trace(:,2)];          % duplicate Y for surface strip
    zS = zeros(size(xS));                    % zero height for 2D plot
    cS = [Time - min(Time), Time - min(Time)]; % color data normalized from zero

    % Create colored line using surface plot (interp edge color)
    hs = surf(xS, yS, zS, cS, 'EdgeColor', 'interp');

    % Set the view to 2D
    view(2);

    % Set colormap and color axis limits
    colormap(gca, 'jet');
    caxis([min(cS(:,1)), max(cS(:,1))]);

    % Set the linewidth of the plot
    set(hs, 'LineWidth', linewidth);

    % Optional: Uncomment to overlay points on the line
    % hold on;
    % plot(Trace(:,1), Trace(:,2), 'ok');
    % hold off;
end

