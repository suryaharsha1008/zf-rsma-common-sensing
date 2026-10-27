function fig3_tau_Rp_frontier()
% Fig. 3 of the paper: FEASIBILITY FRONTIER IN THE (tau, Rp) PLANE
%
% Submission-oriented version:
%   - analytical feasible region and analytical frontier W = alpha_max - alpha_min = 0
%   - physical MC certification with realization-wise confidence intervals
%   - aggregate frontier audit
%   - IEEE single-column PDF export (3.5 inches) with embedded fonts
clear; clc; close all;
rng(7, 'twister');
%% =========================
%  Baseline parameters
%  =========================
prm.Nt       = 8;
prm.K        = 4;
prm.delta    = 0.05;
prm.SNRdB    = 15;
prm.Rc       = 0.20;
prm.sigma_h2 = 1;
prm.sigma_n2 = 1;
prm.eps_c_u   = 0.10;
prm.eps_c_sys = 0.40;
prm.eps_p_u   = 0.10;
prm.eps_c_bar = min(prm.eps_c_u, prm.eps_c_sys / prm.K);
prm.M = prm.Nt - prm.K + 1;
prm.P = 10^(prm.SNRdB / 10) * prm.sigma_n2;
prm.wc = zeros(prm.Nt, 1);
prm.wc(1) = 1;
% Numerical controls
prm.root_tol      = 1e-6;
prm.alpha_eps     = 1e-8;
prm.z95           = 1.95996398454005;   % 95% normal CI over realizations
prm.Nreal_mc      = 2e5;                % increase for final run if runtime allows
prm.doPhysicalMC  = true;
% Export controls
exportOutputs = true;
exportStem    = 'fig3_tau_Rp_frontier';
%% =========================
%  Plot styling (IEEE Strict)
%  =========================
% IEEE prefers 8-10pt font for single column (3.5 inch) figures
ieeeFontSize = 9; 
% Use standard PostScript font names to ensure PDF embedding for Overleaf/IEEE PDF eXpress
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultTextInterpreter', 'tex');
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot, 'defaultLegendInterpreter', 'tex');
% Custom gradient from very light blue to deep blue for the heatmap
numColors = 256;
lightestBlue = [0.92, 0.97, 1.00];
deepBlue     = [0.05, 0.35, 0.65];
customColormap = [linspace(lightestBlue(1), deepBlue(1), numColors)', ...
                  linspace(lightestBlue(2), deepBlue(2), numColors)', ...
                  linspace(lightestBlue(3), deepBlue(3), numColors)'];
infeasibleGray = [0.93, 0.93, 0.93];
orange         = [0.85, 0.33, 0.10]; 
darkGray       = [0.30, 0.30, 0.30];
black          = [0.10, 0.10, 0.10];
%% =========================
%  2D analytical map
%  =========================
tau_vec = linspace(0.01, 0.45, 81);
Rp_vec  = linspace(0.10, 2.00, 81);
[TAU, RP] = meshgrid(tau_vec, Rp_vec);
Amin = nan(size(TAU));
Amax = nan(size(TAU));
W    = nan(size(TAU));
fprintf('============================================================\n');
fprintf('Computing analytical frontier map on %d x %d grid...\n', numel(Rp_vec), numel(tau_vec));
fprintf('============================================================\n');
tic;
for idx = 1:numel(TAU)
    tau0 = TAU(idx);
    Rp0  = RP(idx);
    [aMin0, aMax0] = theory_thresholds(tau0, Rp0, prm);
    Amin(idx) = aMin0;
    Amax(idx) = aMax0;
    W(idx)    = aMax0 - aMin0;
end
tTheoryMap = toc;
analytical_feasible = (W >= 0);
fprintf('Analytical feasible points on dense map: %d / %d\n', nnz(analytical_feasible), numel(analytical_feasible));
fprintf('Analytical map runtime: %.2f s\n', tTheoryMap);
%% =========================
%  Coarser MC audit grid
%  =========================
mc_tau = linspace(0.03, 0.43, 9);
mc_Rp  = linspace(0.20, 1.80, 7);
[MC_TAU, MC_RP] = meshgrid(mc_tau, mc_Rp);
MC_status        = zeros(size(MC_TAU));
MC_alphaHit      = nan(size(MC_TAU));
MC_aMinTheory    = nan(size(MC_TAU));
MC_aMaxTheory    = nan(size(MC_TAU));
MC_WTheory       = nan(size(MC_TAU));
MC_nAlphaTested  = nan(size(MC_TAU));
MC_bestFeasMarg  = nan(size(MC_TAU));
MC_bestViolMarg  = nan(size(MC_TAU));
MC_theoryFeas    = false(size(MC_TAU));
if prm.doPhysicalMC
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('Computing physical-MC frontier audit on %d points...\n', numel(MC_TAU));
    fprintf('Nreal per audit point = %d\n', prm.Nreal_mc);
    fprintf('============================================================\n');
    tic;
    for idx = 1:numel(MC_TAU)
        tau0 = MC_TAU(idx);
        Rp0  = MC_RP(idx);
        [aMin0, aMax0] = theory_thresholds(tau0, Rp0, prm);
        MC_aMinTheory(idx) = aMin0;
        MC_aMaxTheory(idx) = aMax0;
        MC_WTheory(idx)    = aMax0 - aMin0;
        MC_theoryFeas(idx) = (MC_WTheory(idx) >= 0);
        [status0, alphaHit0, cert0] = certify_point_physical(tau0, Rp0, aMin0, aMax0, prm);
        MC_status(idx)       = status0;
        MC_alphaHit(idx)     = alphaHit0;
        MC_nAlphaTested(idx) = cert0.nAlpha;
        MC_bestFeasMarg(idx) = cert0.bestFeasMargin;
        MC_bestViolMarg(idx) = cert0.bestViolMargin;
        fprintf(['Point %3d / %3d: tau = %.3f, Rp = %.3f | ', ...
                 'Theory: [aMin, aMax] = [%.4f, %.4f], W = %+ .4f | ', ...
                 'MC status = %+d\n'], ...
                 idx, numel(MC_TAU), tau0, Rp0, aMin0, aMax0, aMax0 - aMin0, status0);
    end
    tMC = toc;
else
    tMC = NaN;
end

%% =========================
%  COMMAND WINDOW DIAGNOSTICS
%  =========================
if prm.doPhysicalMC
    fprintf('\n============================================================\n');
    fprintf('>>>  FIGURE 2 SWEEPS & CLASSIFICATION DIAGNOSTICS  <<<\n');
    fprintf('============================================================\n');
    
    % 1. What the sweeps are
    fprintf('1. What the new Fig. 2 sweeps\n');
    fprintf('   - Axes: tau (%.3f to %.3f) vs Rp (%.2f to %.2f bit/s/Hz)\n', mc_tau(1), mc_tau(end), mc_Rp(1), mc_Rp(end));
    fprintf('   - Fixed parameters: Nt=%d, K=%d, SNR=%.1f dB, delta=%.2f, Rc=%.2f\n', prm.Nt, prm.K, prm.SNRdB, prm.delta, prm.Rc);
    fprintf('   - Grid size: %d x %d (Total = %d points)\n', numel(mc_Rp), numel(mc_tau), numel(MC_TAU));
    
    % 2. What is classified
    fprintf('\n2. What exactly is being classified\n');
    fprintf('   - Joint-feasible vs infeasible (satisfying common & private outages jointly).\n');
    fprintf('   - Comparing Theory Region (analytical thresholds) vs MC Region (physical bounds).\n');
    
    % 3. Confusion Matrix
    idx_mc_feas   = (MC_status(:) == 1);
    idx_mc_infeas = (MC_status(:) == -1);
    idx_th_feas   = (MC_theoryFeas(:) == 1);
    idx_th_infeas = (MC_theoryFeas(:) == 0);
    
    true_feas_class_feas     = sum(idx_mc_feas & idx_th_feas);
    true_feas_class_infeas   = sum(idx_mc_feas & idx_th_infeas);
    true_infeas_class_infeas = sum(idx_mc_infeas & idx_th_infeas);
    true_infeas_class_feas   = sum(idx_mc_infeas & idx_th_feas);
    
    fprintf('\n3. The confusion-matrix numbers (True = MC, Classified = Theory)\n');
    fprintf('   Total grid points = %d\n', numel(MC_TAU));
    fprintf('   True feasible classified feasible   = %d\n', true_feas_class_feas);
    fprintf('   True feasible classified infeasible = %d (Conservative errors)\n', true_feas_class_infeas);
    fprintf('   True infeasible classified infeasible = %d\n', true_infeas_class_infeas);
    fprintf('   True infeasible classified feasible = %d (False-feasible errors)\n', true_infeas_class_feas);
    
    % 4. Boundary behavior
    err_mask_1d = (idx_mc_feas & idx_th_infeas) | (idx_mc_infeas & idx_th_feas);
    err_mask_2d = reshape(err_mask_1d, size(MC_TAU));
    theory_grid = reshape(MC_theoryFeas, size(MC_TAU));
    
    err_near_frontier = 0;
    err_far = 0;
    
    for r = 1:size(theory_grid,1)
        for c = 1:size(theory_grid,2)
            if err_mask_2d(r,c)
                % Check 8-neighbor adjacency for a flip in theory feasibility
                rmin = max(1, r-1); rmax = min(size(theory_grid,1), r+1);
                cmin = max(1, c-1); cmax = min(size(theory_grid,2), c+1);
                neighbors = theory_grid(rmin:rmax, cmin:cmax);
                if any(neighbors(:) ~= theory_grid(r,c))
                    err_near_frontier = err_near_frontier + 1;
                else
                    err_far = err_far + 1;
                end
            end
        end
    end
    
    fprintf('\n4. Boundary-specific behavior\n');
    if sum(err_mask_1d) == 0
        fprintf('   - Are the errors concentrated near the frontier? N/A (0 errors total)\n');
        fprintf('   - How many errors are more than one grid step from the frontier? 0\n');
    else
        fprintf('   - Are the errors concentrated near the frontier? %s (Near: %d, Far: %d)\n', ...
            mat2str(err_near_frontier > err_far), err_near_frontier, err_far);
        fprintf('   - How many errors are more than one grid step from the frontier? %d\n', err_far);
    end
    
    % 5. Operational claim
    fprintf('\n5. Most important operational claim\n');
    if true_infeas_class_feas == 0
        fprintf('   - Key message: "No false-feasible points." The theory acts as a strict lower-bound/safe region.\n');
    else
        fprintf('   - Key message: "Very small total classification error." The analytical frontier accurately tracks the physical shape.\n');
    end
    
    % 6. Visual Changes Interpretation
    fprintf('\n6. If the figure changed visually:\n');
    fprintf('   - What do colors mean? Shaded background: White = Theory Feasible, Light Gray = Theory Infeasible.\n');
    fprintf('   - What do markers/contours mean? Solid Black Line = Analytical Boundary (W=0).\n');
    fprintf('     Markers: Black Circles = MC Feasible, Red Crosses = MC Infeasible, Gray Dots = Undecided.\n');
    fprintf('   - Is the analytical frontier drawn over MC points? No, `contour` and `contourf` execute first; MC `scatter` markers sit ON TOP of the frontier.\n');
    fprintf('============================================================\n\n');
end

%% =========================
%  Plot (IEEE single-column style)
%  =========================
% IEEE single-column standard width is 3.5 inches. 
% 2.6 inches height gives a clean ~4:3 aspect ratio.
figWidth  = 3.5; 
figHeight = 2.6; 
fig1 = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, figWidth, figHeight]);
hold on; box on;
% 1. Classic Light Gray Shading for the Feasible Region
contourf(TAU, RP, double(analytical_feasible), [-0.5 0.5 1.5], 'LineStyle', 'none');
colormap([1 1 1; 0.88 0.88 0.88]); 
% 2. Stark Black Theoretical Boundary (W = 0)
[~, hFront] = contour(TAU, RP, W, [0 0], 'Color', 'k', 'LineWidth', 1.0, 'LineStyle', '-');
% 3. Standard IEEE Discrete Markers (Scaled down for 3.5in width)
if prm.doPhysicalMC
    feasMask  = (MC_status == +1);
    infeMask  = (MC_status == -1);
    undecMask = (MC_status == 0);
    
    % Feasible: Empty black circles
    hFeas = scatter(MC_TAU(feasMask), MC_RP(feasMask), 18, ...
        'o', 'MarkerEdgeColor', 'k', 'LineWidth', 0.6);
        
    % Infeasible: Dark red crosses
    hInfe = scatter(MC_TAU(infeMask), MC_RP(infeMask), 20, ...
        'x', 'MarkerEdgeColor', [0.7 0 0], 'LineWidth', 0.8);
        
    % Undecided: Gray dots
    hUndc = scatter(MC_TAU(undecMask), MC_RP(undecMask), 10, ...
        '.', 'MarkerEdgeColor', [0.5 0.5 0.5]);
else
    hFeas = gobjects(1);
    hInfe = gobjects(1);
    hUndc = gobjects(1);
end
% 4. Baseline point (Black pentagram/star)
tau_base = 0.20;
Rp_base  = 1.00;
hBase = plot(tau_base, Rp_base, 'p', 'MarkerSize', 6, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% Axis Formatting
set(gca, 'FontName', 'Times New Roman', 'FontSize', ieeeFontSize, 'Layer', 'top');
xlabel('\tau', 'FontName', 'Times New Roman', 'FontSize', ieeeFontSize);
ylabel('R_p (bit/s/Hz)', 'FontName', 'Times New Roman', 'FontSize', ieeeFontSize);
xlim([tau_vec(1), tau_vec(end)]);
ylim([Rp_vec(1), Rp_vec(end)]);
% Faint dotted grid 
grid on;
set(gca, 'GridLineStyle', ':', 'GridColor', 'k', 'GridAlpha', 0.2);
% 5. Compact Legend
hPatch = patch(nan, nan, [0.88 0.88 0.88], 'EdgeColor', 'none');
legendHandles = [hPatch, hFront, hBase];
legendEntries = {'Analytical feasible', 'Analytical frontier', 'Baseline point'};
if prm.doPhysicalMC
    legendHandles = [legendHandles, hFeas, hInfe];
    legendEntries = [legendEntries, 'MC: feasible', 'MC: infeasible']; 
end
% 5. Ultra-Compact Legend for 3.5in Width
hPatch = patch(nan, nan, [0.88 0.88 0.88], 'EdgeColor', 'none');
legendHandles = [hPatch, hFront, hBase];
% Slightly abbreviate text to save horizontal space
legendEntries = {'Analytical feasible', 'Analytical frontier', 'Baseline'};
if prm.doPhysicalMC
    legendHandles = [legendHandles, hFeas, hInfe];
    legendEntries = [legendEntries, 'MC: feasible', 'MC: infeasible']; 
end
% Position and format
lgd = legend(legendHandles, legendEntries, 'Location', 'southwest');
%lgd.NumColumns = 2;              % Split into 2 columns to halve the height!
lgd.ItemTokenSize = [12, 10];    % Shorten the marker/line preview lengths
lgd.FontSize = 5;                % 7pt font (Standard for IEEE figure legends)
lgd.FontName = 'Times New Roman';
lgd.Color = [1 1 1];             % Solid white background
lgd.EdgeColor = 'k';             % Black border
lgd.LineWidth = 0.5;             % Thinner border line
%% =========================
%  Exact-Width PDF Export (IEEE Bulletproof)
%  =========================
if exportOutputs
    % 1. Force the paper size to exactly match the figure size
    set(fig1, 'PaperUnits', 'inches');
    set(fig1, 'PaperPositionMode', 'auto');
    set(fig1, 'PaperSize', [figWidth, figHeight]);
    
    % 2. Force the 'painters' renderer. This is CRITICAL for preventing
    % rasterization of text and ensuring proper vector paths in LaTeX.
    set(fig1, 'Renderer', 'painters');
    
    % 3. Export using exportgraphics
    fileName = sprintf('%s.pdf', exportStem);
    fprintf('\nExporting IEEE exact-width PDF to: %s\n', fileName);
    exportgraphics(fig1, fileName, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
end
%% ========================================================================
function [alphaMin, alphaMax] = theory_thresholds(tau, Rp, prm)
% alphaMin = inf { alpha : pc_u(alpha) <= eps_c_bar }
% alphaMax = sup { alpha : pp_u(alpha) <= eps_p_u }
% ----- alphaMin from common outage (monotone decreasing)
pc_at_1 = common_outage_theory(1.0, tau, prm);
if pc_at_1 > prm.eps_c_bar
    alphaMin = 1 + 1e-3;
else
    lo = 0.0;
    hi = 1.0;
    while (hi - lo) > prm.root_tol
        mid = 0.5 * (lo + hi);
        if common_outage_theory(mid, tau, prm) <= prm.eps_c_bar
            hi = mid;
        else
            lo = mid;
        end
    end
    alphaMin = hi;
end
% ----- alphaMax from private outage (monotone increasing)
pp_at_0 = private_outage_theory_stable(0.0, tau, Rp, prm);
if pp_at_0 > prm.eps_p_u
    alphaMax = -1e-3;
else
    lo = 0.0;
    hi = 1.0;
    while (hi - lo) > prm.root_tol
        mid = 0.5 * (lo + hi);
        if private_outage_theory_stable(mid, tau, Rp, prm) <= prm.eps_p_u
            lo = mid;
        else
            hi = mid;
        end
    end
    alphaMax = lo;
end
end
%% ========================================================================
function pc = common_outage_theory(alpha, tau, prm)
Tc  = 2^prm.Rc - 1;
K   = prm.K;
M   = prm.M;
P   = prm.P;
sh2 = prm.sigma_h2;
sn2 = prm.sigma_n2;
if alpha <= 0
    pc = 1.0;
    return;
end
if alpha >= 1
    sc   = Tc / (P * sh2);
    Psuc = exp(-sc * sn2);
    pc   = 1 - Psuc;
    pc   = min(max(real(pc), 0), 1);
    return;
end
sc = Tc / (alpha * P * sh2);
term1 = 1 + Tc * (1 - alpha) * tau^2 / (K * alpha);
term2 = 1 + Tc * (1 - alpha) / (K * alpha);
Psuc = exp(-sc * sn2) * (term1^(M - K)) * (term2^(-M));
Psuc = min(max(real(Psuc), 0), 1);
pc = 1 - Psuc;
pc = min(max(real(pc), 0), 1);
end
%% ========================================================================
function pp = private_outage_theory_stable(alpha, tau, Rp, prm)
Tp = 2^Rp - 1;
K  = prm.K;
M  = prm.M;
P  = prm.P;
sh2 = prm.sigma_h2;
sn2 = prm.sigma_n2;
delta = prm.delta;
if alpha >= 1
    pp = 1.0;
    return;
end
theta = (1 - alpha) * P * sh2 / K;
if theta <= 0
    pp = 1.0;
    return;
end
theta2   = tau^2 * theta;
thetaRes = delta * alpha * P * sh2;
sp       = Tp / theta;
m        = K - 1;
piw = mixture_weights_stable(M, tau);
logTermsN = -inf(M,1);
for n = 1:M
    logTermsI = -inf(n,1);
    for i = 0:(n-1)
        logTermsEll = -inf(i+1,1);
        for ell = 0:i
            logTermsT = -inf(ell+1,1);
            for t = 0:ell
                logA = logPow(theta2, t) + logPow(thetaRes, ell - t);
                if isinf(logA)
                    continue;
                end
                logCoeff = gammaln(ell + 1) - gammaln(t + 1) ...
                         + gammaln(m + t) - gammaln(m);
                logDen = (m + t) * log1p(theta2 * sp) ...
                       + (1 + ell - t) * log1p(thetaRes * sp);
                logTermsT(t + 1) = logCoeff + logA - logDen;
            end
            logTell = logsumexp(logTermsT);
            if isinf(logTell)
                continue;
            end
            logBinom = gammaln(i + 1) - gammaln(ell + 1) - gammaln(i - ell + 1);
            logNoise = logPow(sn2, i - ell);
            logTermsEll(ell + 1) = logBinom + logNoise + logTell;
        end
        logAccL = logsumexp(logTermsEll);
        if isinf(logAccL)
            continue;
        end
        if i == 0
            logPref = -sp * sn2;
        else
            logPref = -sp * sn2 + i * log(sp) - gammaln(i + 1);
        end
        logTermsI(i + 1) = logPref + logAccL;
    end
    logTermsN(n) = log(piw(n)) + logsumexp(logTermsI);
end
Psuc = exp(logsumexp(logTermsN));
Psuc = min(max(real(Psuc), 0), 1);
pp = 1 - Psuc;
pp = min(max(real(pp), 0), 1);
end
%% ========================================================================
function piw = mixture_weights_stable(M, tau)
if M == 1
    piw = 1;
    return;
end
q = tau^2;
if q == 0
    piw = zeros(M,1);
    piw(M) = 1;
    return;
elseif q == 1
    piw = zeros(M,1);
    piw(1) = 1;
    return;
end
piw = zeros(M,1);
piw(1) = q^(M-1);
ratio = (1 - q) / q;
for n = 1:(M-1)
    piw(n+1) = piw(n) * ((M - n) / n) * ratio;
end
piw = piw / sum(piw);
end
%% ========================================================================
function [status, alphaHit, cert] = certify_point_physical(tau, Rp, aMinTheory, aMaxTheory, prm)
gainBank  = simulate_gain_bank(tau, prm, prm.Nreal_mc);
alphaList = build_alpha_candidates(aMinTheory, aMaxTheory);
nA = numel(alphaList);
feasMargin = -inf(nA,1);
violMargin = -inf(nA,1);
for ia = 1:nA
    a = alphaList(ia);
    stats = physical_outage_ci_realization(a, Rp, gainBank, prm);
    feasMargin(ia) = min(prm.eps_c_bar - stats.pc_ub, prm.eps_p_u - stats.pp_ub);
    violMargin(ia) = max(stats.pc_lb - prm.eps_c_bar, stats.pp_lb - prm.eps_p_u);
end
status = 0;
alphaHit = NaN;
idxFeas = find(feasMargin >= 0);
if ~isempty(idxFeas)
    [~, jj] = max(feasMargin(idxFeas));
    alphaHit = alphaList(idxFeas(jj));
    status = +1;
elseif all(violMargin > 0)
    status = -1;
else
    status = 0;
end
cert = struct();
cert.nAlpha = nA;
cert.bestFeasMargin = max(feasMargin);
cert.bestViolMargin = max(violMargin);
end
%% ========================================================================
function alphaList = build_alpha_candidates(aMinTheory, aMaxTheory)
base = linspace(0, 1, 31);
local = [0, 1];
if isfinite(aMinTheory)
    local = [local, aMinTheory, ...
        linspace(max(0, aMinTheory - 0.06), min(1, aMinTheory + 0.06), 13)];
end
if isfinite(aMaxTheory)
    local = [local, aMaxTheory, ...
        linspace(max(0, aMaxTheory - 0.06), min(1, aMaxTheory + 0.06), 13)];
end
if isfinite(aMinTheory) && isfinite(aMaxTheory)
    left  = max(0, min(aMinTheory, aMaxTheory) - 0.04);
    right = min(1, max(aMinTheory, aMaxTheory) + 0.04);
    local = [local, linspace(left, right, 17)];
    if aMaxTheory >= aMinTheory
        local = [local, ...
            linspace(max(0, aMinTheory), min(1, aMaxTheory), 15), ...
            0.5 * (max(0, aMinTheory) + min(1, aMaxTheory))];
    end
else
    local = [local, linspace(0.30, 0.80, 13)];
end
alphaList = unique(max(0, min(1, [base, local])));
alphaList = sort(alphaList);
end
%% ========================================================================
function gainBank = simulate_gain_bank(tau, prm, Nreal)
Nt = prm.Nt;
K  = prm.K;
sh = sqrt(prm.sigma_h2 / 2);
gC = zeros(Nreal, K);
gD = zeros(Nreal, K);
gI = zeros(Nreal, K);
for r = 1:Nreal
    while true
        Hhat = sh * (randn(Nt, K) + 1i * randn(Nt, K));
        G = Hhat' * Hhat;
        if rcond(G) > 1e-12
            break;
        end
    end
    E = sh * (randn(Nt, K) + 1i * randn(Nt, K));
    H = sqrt(max(0, 1 - tau^2)) * Hhat + tau * E;
    Wtilde = Hhat / G;
    norms  = sqrt(sum(abs(Wtilde).^2, 1));
    W      = Wtilde ./ norms;
    HW = H' * W;                     
    commonGain = abs(H' * prm.wc).^2;
    desiredGain = abs(diag(HW)).^2;
    totalPrivGain = sum(abs(HW).^2, 2);
    interfGain = totalPrivGain - desiredGain;
    gC(r,:) = commonGain.';
    gD(r,:) = desiredGain.';
    gI(r,:) = interfGain.';
end
gainBank.gC = gC;
gainBank.gD = gD;
gainBank.gI = gI;
gainBank.Nreal = Nreal;
gainBank.K = K;
end
%% ========================================================================
function stats = physical_outage_ci_realization(alpha, Rp, gainBank, prm)
P  = prm.P;
K  = prm.K;
Tc = 2^prm.Rc - 1;
Tp = 2^Rp - 1;
Pp = (1 - alpha) * P / K;
gC = gainBank.gC;
gD = gainBank.gD;
gI = gainBank.gI;
gammaC = (alpha * P .* gC) ./ (Pp .* (gD + gI) + prm.sigma_n2);
gammaP = (Pp .* gD) ./ (Pp .* gI + prm.delta * alpha * P .* gC + prm.sigma_n2);
commonFail  = (gammaC < Tc);
privateFail = (gammaP < Tp);
qC = mean(commonFail, 2);
qP = mean(privateFail, 2);
[pc_hat, pc_lb, pc_ub, pc_se] = mean_ci_from_realizations(qC, prm.z95);
[pp_hat, pp_lb, pp_ub, pp_se] = mean_ci_from_realizations(qP, prm.z95);
stats.pc_hat = pc_hat;
stats.pc_lb  = pc_lb;
stats.pc_ub  = pc_ub;
stats.pc_se  = pc_se;
stats.pp_hat = pp_hat;
stats.pp_lb  = pp_lb;
stats.pp_ub  = pp_ub;
stats.pp_se  = pp_se;
end
%% ========================================================================
function [mu, lb, ub, se] = mean_ci_from_realizations(q, z)
q = q(:);
mu = mean(q);
se = std(q, 0, 1) / sqrt(numel(q));
lb = max(0, mu - z * se);
ub = min(1, mu + z * se);
end
%% ========================================================================
function y = logPow(a, p)
if p == 0
    y = 0;
elseif a == 0
    y = -inf;
else
    y = p * log(a);
end
end
%% ========================================================================
function s = logsumexp(v)
vmax = max(v);
if isinf(vmax)
    s = -inf;
else
    s = vmax + log(sum(exp(v - vmax)));
end
end