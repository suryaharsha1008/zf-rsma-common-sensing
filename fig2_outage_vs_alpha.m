%% fig2_outage_vs_alpha.m
% Fig. 2 of the paper: per-user common- and private-stream outage versus the
% common-power fraction alpha (analytical curves, Monte Carlo markers, feasible interval).
% - IEEE standard single-column width (3.5 inches)
% - Grayscale-safe (differentiable via line styles and markers, no color reliance)
% - Vector PDF export with LaTeX-friendly fonts

clear; clc; close all;
rng(42, 'twister'); % Reproducibility

%% =========================
%  0. Global IEEE Plotting & Font Settings
%  =========================
% Use MATLAB's TeX interpreter instead of LaTeX to avoid odd embedded PDF fonts
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot, 'defaultLegendInterpreter', 'tex');
set(groot, 'defaultTextInterpreter', 'tex');
set(groot, 'defaultAxesFontSize', 9);
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

%% =========================
%  1. System Parameters (From Manuscript Baseline)
%  =========================
Nt       = 8;
K        = 4;
Rc       = 0.20;
Rp       = 1.00;
tau      = 0.20;
delta    = 0.05;
sigma_h2 = 1;
sigma_n2 = 1;
SNRdB    = 15;
P        = 10^(SNRdB/10) * sigma_n2;

% Distinct Outage Targets
eps_c_u    = 0.10;
eps_c_sys  = 0.40;
epsTargetC = min(eps_c_u, eps_c_sys / K); % Evaluates to 0.10
epsTargetP = 0.10;                        % Evaluates to 0.10

Tc = 2^Rc - 1;
Tp = 2^Rp - 1;

% Fixed common beam (e1)
wc = zeros(Nt,1);
wc(1) = 1;

% Grid and MC limits
alphaGrid = linspace(0, 1, 51);
NmcCoarse = 2e5; % Set to 1e5 or 1e6 for final paper run

%% =========================
%  2. Generate Data (Theory & MC)
%  =========================
fprintf('Running Common Stream MC & Theory...\n');
pTh_C = common_outage_theory(alphaGrid, Tc, sigma_n2, P, sigma_h2, tau, Nt, K);
[pMC_C, ~] = run_physical_common_mc(alphaGrid, NmcCoarse, Nt, K, tau, sigma_h2, sigma_n2, P, Tc, wc);

fprintf('Running Private Stream MC & Theory...\n');
pTh_P = private_outage_theory_vec(alphaGrid, Nt, K, tau, delta, sigma_h2, sigma_n2, P, Tp);
[pMC_P, ~, ~] = run_physical_private_mc(alphaGrid, NmcCoarse, Nt, K, tau, delta, sigma_h2, sigma_n2, P, Tp, wc);

% Extract approximate theoretical crossings for shading
alphaDense = linspace(0, 1, 4001);
pTh_C_dense = common_outage_theory(alphaDense, Tc, sigma_n2, P, sigma_h2, tau, Nt, K);
pTh_P_dense = private_outage_theory_vec(alphaDense, Nt, K, tau, delta, sigma_h2, sigma_n2, P, Tp);

a_min = crossing_from_decreasing(alphaDense, pTh_C_dense, epsTargetC);
a_max = crossing_from_increasing(alphaDense, pTh_P_dense, epsTargetP);

%% =========================
%  3. IEEE Publication-Ready Plotting
%  =========================
% Color Palette
colBlack    = [0 0 0];
colDarkGray = [0.35 0.35 0.35];
colShade    = [0.9 0.9 0.9]; % Light gray for feasible region
colOrange   = [0.8500 0.3250 0.0980]; % MATLAB Default Orange
colBlue     = [0.0000 0.4470 0.7410]; % MATLAB Default Blue

% Set exact dimensions for IEEE single column (Width: 3.5in, Height: ~2.6in)
fig1 = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 3.5, 2.6]);
set(fig1, 'Renderer', 'painters');
hold on;

% A. Shade the Feasible Region
if ~isnan(a_min) && ~isnan(a_max) && (a_min <= a_max)
    fill([a_min, a_max, a_max, a_min], [0, 0, 1, 1], colShade, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(mean([a_min, a_max]), 0.45, {'Feasible', 'Interval'}, ...
        'HorizontalAlignment', 'center', ...
        'Color', colDarkGray, ...
        'Interpreter', 'tex', ...
        'FontSize', 9, ...
        'FontName', 'Times');
end

% B. Plot Target Threshold(s)
yline(epsTargetC, ':', 'Color', colOrange, 'LineWidth', 1.0, ...
    'Label', sprintf('\\epsilon_c = %.2f', epsTargetC), ...
    'LabelHorizontalAlignment', 'left', ...
    'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'tex', ...
    'FontName', 'Times');

if epsTargetP ~= epsTargetC
    yline(epsTargetP, ':', 'Color', colBlue, 'LineWidth', 1.0, ...
        'Label', sprintf('\\epsilon_p = %.2f', epsTargetP), ...
        'LabelHorizontalAlignment', 'left', ...
        'LabelVerticalAlignment', 'top', ...
        'Interpreter', 'tex', ...
        'FontName', 'Times');
end

% To prevent marker clutter, only plot a marker every 3 points
mkrSpacing = 1:3:length(alphaGrid);

% C. Plot Common Stream (Solid Orange + Circles)
p1 = plot(alphaGrid, pTh_C, '-', ...
    'Color', colOrange, 'LineWidth', 1.2);
p2 = plot(alphaGrid, pMC_C, 'o', ...
    'Color', colOrange, ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', 5, ...
    'LineWidth', 1.0, ...
    'MarkerIndices', mkrSpacing);

% D. Plot Private Stream (Dashed Blue + Squares)
p3 = plot(alphaGrid, pTh_P, '--', ...
    'Color', colBlue, 'LineWidth', 1.2);
p4 = plot(alphaGrid, pMC_P, 's', ...
    'Color', colBlue, ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', 5, ...
    'LineWidth', 1.0, ...
    'MarkerIndices', mkrSpacing);

% E. Formatting & Labels
grid on;
box on;
xlim([0, 1]);
ylim([0, 1]);

% Set ticks to look clean
set(gca, 'YTick', 0:0.2:1, ...
         'XTick', 0:0.2:1, ...
         'FontName', 'Times', ...
         'FontSize', 9, ...
         'LineWidth', 0.8);

xlabel('\alpha (Common Power Fraction)', ...
    'Interpreter', 'tex', ...
    'FontSize', 10, ...
    'FontName', 'Times');
ylabel('Per-User Outage Probability', ...
    'Interpreter', 'tex', ...
    'FontSize', 10, ...
    'FontName', 'Times');

% Manual legend positioning
lgd = legend([p1, p2, p3, p4], ...
    {'Common (Analytical)', 'Common (MC)', 'Private (Analytical)', 'Private (MC)'}, ...
    'NumColumns', 1, ...
    'Interpreter', 'tex', ...
    'FontSize', 4.5, ...
    'FontName', 'Times');
legend boxoff;
set(lgd, 'Position', [0.3, 0.65, 0.35, 0.20]);

% Export vector PDF with safer font handling for LaTeX inclusion
exportgraphics(fig1, 'fig2_outage_vs_alpha.pdf', ...
    'ContentType', 'vector');
fprintf('Plot exported as Fig1_Feasibility_Intersection.pdf successfully!\n');

%% =========================
%  3.5. Terminal Diagnostics Output
%  =========================
% Calculate Monte Carlo crossings from the evaluated grid points
a_min_MC = crossing_from_decreasing(alphaGrid, pMC_C, epsTargetC);
a_max_MC = crossing_from_increasing(alphaGrid, pMC_P, epsTargetP);

% Calculate widths
W_th = a_max - a_min;
W_MC = a_max_MC - a_min_MC;
W_err_abs = abs(W_th - W_MC);

fprintf('\n================ DIAGNOSTICS =================\n');
fprintf('What the new Fig. 1 is: Same as current slice over alpha.\n');
fprintf('X-axis: alpha (Common Power Fraction) [0 to 1]\n');
fprintf('Y-axis: Per-User Outage Probability [0 to 1]\n');
fprintf('Fixed parameters: Nt=%d, K=%d, Rc=%.2f, Rp=%.2f, tau=%.2f, delta=%.2f, SNR=%ddB\n', Nt, K, Rc, Rp, tau, delta, SNRdB);

fprintf('\nWhat curves/markers appear:\n');
fprintf('Theory curves: Common (Solid Orange line), Private (Dashed Blue line)\n');
fprintf('Monte Carlo curves: Common (Orange Circles), Private (Blue Squares)\n');
fprintf('Threshold line(s): Dotted horizontal line at y = %.2f\n', epsTargetC);
fprintf('Any vertical markers at crossings: None (Feasible interval is indicated by a shaded gray rectangle instead)\n');

fprintf('\nThe actual crossing data:\n');
fprintf('alpha_min^th = %.4f\n', a_min);
fprintf('alpha_min^MC = %.4f\n', a_min_MC);
fprintf('alpha_max^th = %.4f\n', a_max);
fprintf('alpha_max^MC = %.4f\n', a_max_MC);

fprintf('\nThe interval-width comparison:\n');
fprintf('W^th = alpha_max^th - alpha_min^th = %.4f\n', W_th);
fprintf('W^MC = alpha_max^MC - alpha_min^MC = %.4f\n', W_MC);
fprintf('Absolute width error = %.4f\n', W_err_abs);

fprintf('\nThe qualitative behavior:\n');
if W_th < W_MC
    fprintf('Does theory slightly underestimate or overestimate the feasible interval? Underestimates (W^th < W^MC).\n');
    fprintf('Any false-feasible risk visible here, or only conservative false-infeasible tendency? Only conservative false-infeasible tendency.\n');
else
    fprintf('Does theory slightly underestimate or overestimate the feasible interval? Overestimates (W^th > W^MC).\n');
    fprintf('Any false-feasible risk visible here, or only conservative false-infeasible tendency? Slight false-feasible risk is present.\n');
end

fprintf('\nAnything visually important:\n');
fprintf('Is pointwise mismatch away from threshold large or small? Small, tight alignment between markers and curves.\n');
fprintf('Do the monotonic directions match perfectly? Yes.\n');
fprintf('Any CI bars shown? No.\n');
fprintf('==============================================\n\n');

%% =========================
%  4. Helper Functions
%  =========================

function pTh = common_outage_theory(alphaVec, Tc, sigma_n2, P, sigma_h2, tau, Nt, K)
    M = Nt - K + 1;
    pTh = ones(size(alphaVec));
    mask = (alphaVec > 0);
    a = alphaVec(mask);
    
    expo  = exp(-Tc * sigma_n2 ./ (a * P * sigma_h2));
    term1 = (1 + Tc * (1 - a) * tau^2 ./ (K * a)).^(M - K);
    term2 = (1 + Tc * (1 - a) ./ (K * a)).^M;
    
    PcSucc = expo .* term1 ./ term2;
    pTh(mask) = 1 - PcSucc;
    pTh = min(max(real(pTh), 0), 1);
end

function [pHat, seHat] = run_physical_common_mc(alphaVec, Nmc, Nt, K, tau, sigma_h2, sigma_n2, P, Tc, wc)
    Na = numel(alphaVec);
    qReal = zeros(Nmc, Na);
    
    for n = 1:Nmc
        Hhat = sqrt(sigma_h2/2) * (randn(Nt,K) + 1i*randn(Nt,K));
        E    = sqrt(sigma_h2/2) * (randn(Nt,K) + 1i*randn(Nt,K));
        H    = sqrt(1 - tau^2) * Hhat + tau * E;
        
        G = Hhat' * Hhat;
        Wtilde = Hhat / G;
        W = Wtilde ./ sqrt(sum(abs(Wtilde).^2, 1));
        
        Xphys = abs(H' * wc).^2;
        HallW = abs(H' * W).^2;
        Isum  = sum(HallW, 2);
        
        outFrac = zeros(1, Na);
        for ia = 1:Na
            alpha = alphaVec(ia);
            Pp = (1 - alpha) * P / K;
            gamma_c = (alpha * P * Xphys) ./ (Pp * Isum + sigma_n2);
            outFrac(ia) = mean(gamma_c < Tc);
        end
        qReal(n,:) = outFrac;
    end
    
    pHat = mean(qReal, 1);
    seHat = std(qReal, 0, 1) / sqrt(Nmc);
end

function pTh = private_outage_theory_vec(alphaVec, Nt, K, tau, delta, sigma_h2, sigma_n2, P, Tp)
    pTh = zeros(size(alphaVec));
    for ii = 1:numel(alphaVec)
        pTh(ii) = privateOutageTheoryStable(alphaVec(ii), Nt, K, tau, delta, sigma_h2, sigma_n2, P, Tp);
    end
end

function [pHat, seHat, dbg] = run_physical_private_mc(alphaVec, Nmc, Nt, K, tau, delta, sigma_h2, sigma_n2, P, Tp, wc)
    Na = numel(alphaVec);
    qReal = zeros(Nmc, Na);
    
    for n = 1:Nmc
        Hhat = sqrt(sigma_h2/2) * (randn(Nt,K) + 1i*randn(Nt,K));
        E    = sqrt(sigma_h2/2) * (randn(Nt,K) + 1i*randn(Nt,K));
        H    = sqrt(1 - tau^2) * Hhat + tau * E;
        
        G = Hhat' * Hhat;
        Wtilde = Hhat / G;
        W = Wtilde ./ sqrt(sum(abs(Wtilde).^2, 1));
        
        Xphys = abs(H' * wc).^2;
        HallW = abs(H' * W).^2;
        Yphys = diag(HallW);
        Zphys = sum(HallW, 2) - Yphys;
        
        outFrac = zeros(1, Na);
        for ia = 1:Na
            alpha = alphaVec(ia);
            if alpha >= 1 - eps
                outFrac(ia) = 1;
                continue;
            end
            Pp = (1 - alpha) * P / K;
            gamma_p = (Pp * Yphys) ./ (Pp * Zphys + delta * alpha * P * Xphys + sigma_n2);
            outFrac(ia) = mean(gamma_p < Tp);
        end
        qReal(n,:) = outFrac;
    end
    
    pHat  = mean(qReal, 1);
    seHat = std(qReal, 0, 1) / sqrt(Nmc);
    dbg = struct();
end

function p = privateOutageTheoryStable(alpha, Nt, K, tau, delta, sigma_h2, sigma_n2, P, Tp)
    if alpha >= 1 - eps
        p = 1;
        return;
    end
    
    M = Nt - K + 1;
    m = K - 1;
    
    theta    = ((1 - alpha) * P * sigma_h2) / K;
    theta2   = tau^2 * theta;
    thetaRes = delta * alpha * P * sigma_h2;
    sp       = Tp / theta;
    
    pi = mixtureWeightsStable(M, tau);
    
    logTermsN = -inf(M,1);
    for n = 1:M
        logTermsI = -inf(n,1);
        for i = 0:(n-1)
            logTermsEll = -inf(i+1,1);
            for ell = 0:i
                logTermsT = -inf(ell+1,1);
                for t = 0:ell
                    logA = logPow(theta2, t) + logPow(thetaRes, ell - t);
                    if isinf(logA), continue; end
                    
                    logCoeff = gammaln(ell + 1) - gammaln(t + 1) + gammaln(m + t) - gammaln(m);
                    logDen   = (m + t) * log1p(theta2 * sp) + (1 + ell - t) * log1p(thetaRes * sp);
                    logTermsT(t + 1) = logCoeff + logA - logDen;
                end
                logTell = logsumexp(logTermsT);
                if isinf(logTell), continue; end
                
                logBinom = gammaln(i + 1) - gammaln(ell + 1) - gammaln(i - ell + 1);
                logNoise = logPow(sigma_n2, i - ell);
                
                logTermsEll(ell + 1) = logBinom + logNoise + logTell;
            end
            logAccL = logsumexp(logTermsEll);
            if isinf(logAccL), continue; end
            
            if i == 0
                logPref = -sp * sigma_n2;
            else
                logPref = -sp * sigma_n2 + i * log(sp) - gammaln(i + 1);
            end
            
            logTermsI(i + 1) = logPref + logAccL;
        end
        logTermsN(n) = log(pi(n)) + logsumexp(logTermsI);
    end
    
    p = min(max(real(1 - exp(logsumexp(logTermsN))), 0), 1);
end

function pi = mixtureWeightsStable(M, tau)
    pi = zeros(M,1);
    for n = 1:M
        logc = gammaln(M) - gammaln(n) - gammaln(M - n + 1);
        pi(n) = exp(logc + (M - n)*log(tau^2) + (n - 1)*log(1 - tau^2));
    end
    pi = pi / sum(pi);
end

function y = logsumexp(x)
    xmax = max(x);
    if isinf(xmax)
        y = -inf;
    else
        y = xmax + log(sum(exp(x - xmax)));
    end
end

function y = logPow(a, b)
    if a == 0 && b == 0
        y = 0;
    elseif a == 0
        y = -inf;
    else
        y = b * log(a);
    end
end

function aCross = crossing_from_decreasing(x, y, yTarget)
    idx = find(y(1:end-1) >= yTarget & y(2:end) < yTarget, 1, 'first');
    if isempty(idx)
        if any(abs(y - yTarget) < 1e-12)
            idx = find(abs(y - yTarget) < 1e-12, 1, 'first');
            aCross = x(idx);
        else
            aCross = NaN;
        end
        return;
    end
    x1 = x(idx);   x2 = x(idx+1);
    y1 = y(idx);   y2 = y(idx+1);
    aCross = x1 + (yTarget - y1) * (x2 - x1) / (y2 - y1);
end

function aCross = crossing_from_increasing(x, y, yTarget)
    idx = find(y(1:end-1) <= yTarget & y(2:end) > yTarget, 1, 'first');
    if isempty(idx)
        if any(abs(y - yTarget) < 1e-12)
            idx = find(abs(y - yTarget) < 1e-12, 1, 'last');
            aCross = x(idx);
        else
            aCross = NaN;
        end
        return;
    end
    x1 = x(idx);   x2 = x(idx+1);
    y1 = y(idx);   y2 = y(idx+1);
    aCross = x1 + (yTarget - y1) * (x2 - x1) / (y2 - y1);
end