function B = average_colN(A, N)
% average_colN - Averages the input matrix A along its columns to produce 
% a new matrix B with N columns.
%
% Inputs:
%   A - Input matrix of size (rows x cols)
%   N - Desired number of columns in the output matrix B
%
% Output:
%   B - Output matrix of size (rows x N), where each column is the average 
%       of a group of columns from A

    [rows, cols] = size(A);  % Get size of input matrix

    % If N is greater than number of columns, keep original number of columns
    if N > cols
        N = cols;
    end

    new_cols = N;               % Set number of output columns
    n = round(cols / new_cols); % Number of columns to average per group

    B = zeros(rows, new_cols);  % Preallocate output matrix

    for i = 1:new_cols
        % Calculate column index range for current group
        col_start = (i - 1) * n + 1;
        col_end = min(i * n, cols);  % Make sure not to exceed actual columns
        col_range = col_start:col_end;

        % Compute the mean of the selected columns, row-wise
        B(:, i) = mean(A(:, col_range), 2);
    end
end
