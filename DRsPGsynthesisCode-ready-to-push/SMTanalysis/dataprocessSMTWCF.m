%% dataprocessSMTWCF.m
%
% This script processes single-molecule tracking (SMT) data of FtsW
% under the condition 'FtsW-S1-Cef'. It fits the speed distribution
% using both single and double log-normal models via CDF fitting.
%
% The script includes:
%   1. Data loading and filtering by movement state
%   2. Empirical CDF calculation
%   3. Curve fitting to log-normal models
%   4. Plotting of CDF fits and residuals
%   5. Bootstrapping to estimate fitting uncertainties
%   6. PDF reconstruction for histogram comparison
%   7. Saving all fitting results into a structured MAT file
%

clear;
clc;

% ---------------------- Parameters ----------------------
Nbin   = 21;                    % Number of bins for log-spaced CDF
LowB   = 0;                     % Lower bound of log10(speed)
HighB  = 2;                     % Upper bound of log10(speed)
N      = 200;                   % Bootstrap iterations
xbin   = logspace(LowB, HighB, Nbin);
V_th   = 1;                     % Speed threshold (nm/s)
condition = 'FtsW-S1-Cef';      % Experiment label

% ---------------------- Load Data ----------------------
load('FtsW-all.mat');           % Load the SMT data

% Separate moving (flag = 1) and stationary (flag = 3) FtsW particles
idv = find(DataSMT(:,6)==3 & DataSMT(:,7)==1);  % Moving
dataV = DataSMT(idv,:);

ids = find(DataSMT(:,6)==3 & DataSMT(:,7)==3);  % Stationary
dataS = DataSMT(ids,:);

% Compute average for both groups
MeanVs = [mean(dataV); mean(dataS)];

% Filter velocities above threshold
Vx = dataV(:,1);
Vx = Vx(Vx > V_th);

% ---------------------- Empirical CDF ----------------------
CDF_real = CDF_logCalc(Vx, xbin);
Result.CDF_real = CDF_real;

% ---------------------- Single-population Log-normal Fit ----------------------
P_fit1 = lsqcurvefit(@logn1cdf, [0.5,2,0.5], xbin, CDF_real, [0.99,0,0], [1,5,5]);
CDF_real_fit1 = logn1cdf(P_fit1, xbin);
Residual1 = CDF_real_fit1 - CDF_real;

% Mean and std from log-normal parameters
Fit_Vdirx1 = [exp(P_fit1(2)+P_fit1(3)^2/2), ...
              sqrt((exp(P_fit1(3)^2)-1)*exp(2*P_fit1(2)+P_fit1(3)^2))];

% Store results
Result.CDF_real_fit1 = CDF_real_fit1;
Result.P_fit1 = P_fit1;
Result.Residual1 = Residual1;
Result.FIt_Vdirx1 = Fit_Vdirx1;

% ---------------------- Two-population Log-normal Fit ----------------------
P_fit = lsqcurvefit(@logn2cdf, [0.5,1.5,0.5,3,1], xbin, CDF_real, [0,0,0,0,0], [1,5,2,5,2]);

% Ensure pop1 is slower than pop2
if P_fit(2) > P_fit(4)
    P_fit2(1) = 1 - P_fit(1);
    P_fit2(2:3) = P_fit(4:5);
    P_fit2(4:5) = P_fit(2:3);
else
    P_fit2 = P_fit;
end

CDF_real_fit2 = logn2cdf(P_fit2, xbin);
Residual2 = CDF_real_fit2 - CDF_real;

Fit_Vdirx2 = [...
    P_fit2(1)*100, exp(P_fit2(2)+P_fit2(3)^2/2), ...
    sqrt((exp(P_fit2(3)^2)-1)*exp(2*P_fit2(2)+P_fit2(3)^2)), ...
    (1-P_fit2(1))*100, exp(P_fit2(4)+P_fit2(5)^2/2), ...
    sqrt((exp(P_fit2(5)^2)-1)*exp(2*P_fit2(4)+P_fit2(5)^2))
];

% Store results
Result.CDF_real_fit2 = CDF_real_fit2;
Result.Residual2 = Residual2;
Result.Fit_Vdirx2 = Fit_Vdirx2;
Result.P_fit2 = P_fit2;

% ---------------------- Plot CDF and Residuals ----------------------
h1 = figure('position', [100, 100, 1200, 800]);

subplot('position', [0.15, 0.4, 0.8, 0.5]);
semilogx(xbin, CDF_real_fit1, 'b', 'LineWidth', 2);
hold on;
scatter(xbin, CDF_real, 72, 'k', 'filled');
semilogx(xbin, CDF_real_fit2, 'm', 'LineWidth', 2);

text(1.1, 0.8, ['V1 = ' num2str(Fit_Vdirx2(2),4) ' nm/s; P1 = ' num2str(Fit_Vdirx2(1)/100,3)], 'FontSize', 18);
text(1.1, 0.7, ['V2 = ' num2str(Fit_Vdirx2(5),4) ' nm/s; P2 = ' num2str(Fit_Vdirx2(4)/100,3)], 'FontSize', 18);
text(1.1, 0.9, ['V = ' num2str(Fit_Vdirx1(1),4) ' nm/s'], 'FontSize', 18);

title(condition);
ylabel('CDF', 'FontSize', 24);
legend('Single fit', 'Data', 'Double fit', 'Location', 'best');
set(gca, 'FontSize', 20);

subplot('position', [0.15, 0.1, 0.8, 0.2]);
semilogx(xbin, Residual1, 'b', 'LineWidth', 2); hold on;
semilogx(xbin, Residual2, 'm', 'LineWidth', 2);
yline(0, ':', 'LineWidth', 2);
ylabel('Residual', 'FontSize', 24);
xlabel('Speed_{FtsW} (nm/s)', 'FontSize', 24);
set(gca, 'FontSize', 20);

% ---------------------- Bootstrapping ----------------------
[bootV, bootVdir] = bootstrp(N, @mean, Vx);

for idB = 1:N
    VdirxT = Vx(bootVdir(:,idB));
    CDF = CDF_logCalc(VdirxT, xbin);
    CDFboot(:, idB) = CDF;

    P_fit1B = lsqcurvefit(@logn1cdf, [1, 2, 1], xbin, CDF, [0.99, 0, 0], [1, 5, 5]);
    P_fit2Bt = lsqcurvefit(@logn2cdf, [0.5, 2, 0.5, 3, 1], xbin, CDF, [0, 0, 0, 0, 0], [1, 5, 2, 5, 2]);

    if P_fit2Bt(2) > P_fit2Bt(4)
        P_fit2B(1) = 1 - P_fit2Bt(1);
        P_fit2B(2:3) = P_fit2Bt(4:5);
        P_fit2B(4:5) = P_fit2Bt(2:3);
    else
        P_fit2B = P_fit2Bt;
    end

    P_fit1_boot(idB,:) = P_fit1B;
    P_fit2_boot(idB,:) = P_fit2B;

    Fit_Vdirx1_boot(idB,:) = [exp(P_fit1B(2)+P_fit1B(3)^2/2), ...
                              sqrt((exp(P_fit1B(3)^2)-1)*exp(2*P_fit1B(2)+P_fit1B(3)^2))];

    Fit_Vdirx2_boot(idB,:) = [...
        P_fit2B(1)*100, exp(P_fit2B(2)+P_fit2B(3)^2/2), ...
        sqrt((exp(P_fit2B(3)^2)-1)*exp(2*P_fit2B(2)+P_fit2B(3)^2)), ...
        (1 - P_fit2B(1))*100, exp(P_fit2B(4)+P_fit2B(5)^2/2), ...
        sqrt((exp(P_fit2B(5)^2)-1)*exp(2*P_fit2B(4)+P_fit2B(5)^2))];
end

Result.Vdirx1_boot = Fit_Vdirx1_boot;
Result.Vdirx2_boot = Fit_Vdirx2_boot;
Result.P_fit1_boot = P_fit1_boot;
Result.P_fit2_boot = P_fit2_boot;
Result.Vdirx1_sem = std(Fit_Vdirx1_boot, [], 1);
Result.Vdirx2_sem = std(Fit_Vdirx2_boot, [], 1);

% ---------------------- Histogram and PDF Plotting ----------------------
xbin_center = xbin(1:end-1) + diff(xbin)/exp(1);
Interval = diff(xbin);
HisRes = histcounts(Vx, xbin, 'Normalization', 'probability');

% --- Single-population PDF
P = Result.P_fit1(1);
mu = Result.P_fit1(2);
sigma = Result.P_fit1(3);
Pdf_pop1 = lognpdf(xbin_center, mu, sigma) * P;
His_pop1 = Pdf_pop1 .* Interval;
Result.His_1pop_final = [xbin_center', HisRes', His_pop1'];

% --- Double-population PDF
P1 = Result.P_fit2(1);
mu1 = Result.P_fit2(2); sigma1 = Result.P_fit2(3);
mu2 = Result.P_fit2(4); sigma2 = Result.P_fit2(5);
Y1_pdf = lognpdf(xbin_center, mu1, sigma1) * P1;
Y2_pdf = lognpdf(xbin_center, mu2, sigma2) * (1 - P1);
Y_P = (Y1_pdf + Y2_pdf) .* Interval;
Result.His_2pop_final = [xbin_center', HisRes', Y1_pdf.*Interval', Y2_pdf.*Interval', Y_P'];

% --- Plot PDF
figure;
semilogx(Result.His_2pop_final(:,1), Result.His_2pop_final(:,2)); hold on;
semilogx(Result.His_2pop_final(:,1), Result.His_2pop_final(:,3));
semilogx(Result.His_2pop_final(:,1), Result.His_2pop_final(:,4));
semilogx(Result.His_2pop_final(:,1), Result.His_2pop_final(:,5));

% ---------------------- Save Results ----------------------
save([condition '.mat'], 'Result');
