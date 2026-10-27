%% fig4_high_snr_frontier.m
% Fig. 4 of the paper.
% Generates a single-column (3.5 inch), stacked two-panel figure:
% Panel (a): Feasible-interval width vs SNR (asymptotic plateau).
% Panel (b): Asymptotic outage curves (perfect vs imperfect SIC).
%
% Includes strict IEEE formatting and LaTeX font embedding.
clear; clc; close all;
rng(17, 'twister');
%% ================== 1. IEEE Standard Font & Interpreter Rules =========
% Force MATLAB to use Times and its internal LaTeX engine globally
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 9); % 9pt is standard for IEEE figures
set(groot, 'defaultTextFontSize', 9);
% Force LaTeX interpreter so Overleaf recognizes the embedded fonts
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
%% ===================== 2. Baseline parameters ============================
par.Nt        = 8;
par.K         = 4;
par.tau       = 0.20;
par.delta     = 0.05; % Baseline residual SIC
par.Rc        = 0.20;
par.Rp        = 1.00;
par.sigma_h2  = 1.0;
par.sigma_n2  = 1.0;
par.eps_c_u   = 0.10;
par.eps_c_sys = 0.40;
par.eps_p_u   = 0.10;
% Derived quantities
par.eps_c_bar = min(par.eps_c_u, par.eps_c_sys / par.K);
par.M         = par.Nt - par.K + 1;
par.m         = par.K - 1;
par.Tc        = 2^par.Rc - 1;
par.Tp        = 2^par.Rp - 1;
par.root_tol  = 1e-8;
par.alpha_eps = 1e-12;
% Targets
qc = 1 - par.eps_c_bar;
qp = 1 - par.eps_p_u;
%% ======================== 3. Grids & Styling =============================
SNRdB = 0:1:35;
Pgrid = 10.^(SNRdB/10) * par.sigma_n2;
alpha_grid = linspace(1e-4, 1 - 1e-4, 4001);
delta_list = [0.00, par.delta];
skyBlue  = [0.35, 0.78, 0.98];
orange   = [0.90, 0.45, 0.12];
black    = [0.10, 0.10, 0.10];
gray     = [0.55, 0.55, 0.55];
%% ================== 4. Compute Data: Finite SNR Sweep ====================
fprintf('Computing Panel (a): Finite-SNR sweep...\n');
alpha_min = nan(size(Pgrid));
alpha_max = nan(size(Pgrid));
W         = nan(size(Pgrid));
is_feas   = false(size(Pgrid));
for ii = 1:numel(Pgrid)
    P = Pgrid(ii);
    alpha_min(ii) = find_alpha_min(P, par);
    alpha_max(ii) = find_alpha_max(P, par);
    if isfinite(alpha_min(ii)) && isfinite(alpha_max(ii)) && (alpha_min(ii) <= alpha_max(ii))
        is_feas(ii) = true;
        W(ii) = alpha_max(ii) - alpha_min(ii);
    end
end
% High-SNR limits for Panel A baseline
alpha_min_inf_base = find_alpha_min_inf(par);
alpha_max_inf_base = find_alpha_max_inf(par);
if isfinite(alpha_min_inf_base) && isfinite(alpha_max_inf_base) && (alpha_min_inf_base <= alpha_max_inf_base)
    W_inf_base = alpha_max_inf_base - alpha_min_inf_base;
else
    W_inf_base = NaN;
end
first_feas = find(is_feas, 1, 'first');
%% ================== 5. Compute Data: Asymptotic Bifurcation =================
fprintf('Computing Panel (b): Asymptotic Bifurcation...\n');
pi_n = mixture_weights_stable(par.M, par.tau);
Psi_c = psi_c_vec(alpha_grid, par);
Psi_p = zeros(numel(delta_list), numel(alpha_grid));
for id = 1:numel(delta_list)
    Psi_p(id, :) = psi_p_vec(alpha_grid, delta_list(id), par, pi_n);
end
pc_inf = 1 - Psi_c;
pp_inf = 1 - Psi_p;
% Find crossing points for annotations
alpha_min_inf_bifurc = threshold_from_increasing(alpha_grid, Psi_c, qc);
alpha_max_inf_bifurc = nan(size(delta_list));
for id = 1:numel(delta_list)
    if all(Psi_p(id, :) >= qp)
        alpha_max_inf_bifurc(id) = 1.0;
    elseif any(Psi_p(id, :) >= qp)
        alpha_max_inf_bifurc(id) = threshold_from_decreasing(alpha_grid, Psi_p(id, :), qp);
    end
end
%% ================== 6. Generate IEEE Formatted Figure =======================
% Dimensions tuned for a single IEEE column
figWidth  = 3.5; % 3.5 inches wide
figHeight = 4.5; % 4.5 inches tall (accommodates 2 panels comfortably)
fig = figure('Color', 'w', ...
             'Units', 'inches', ...
             'Position', [1, 1, figWidth, figHeight], ...
             'PaperUnits', 'inches', ...
             'PaperPositionMode', 'auto'); % Critical for exact size export
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% -----------------------------------------------------------------------
% Panel (a): Feasible width versus SNR
% -----------------------------------------------------------------------
ax1 = nexttile(1);
hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');
set(ax1, 'LineWidth', 0.8, 'GridAlpha', 0.15);
h1 = plot(ax1, SNRdB(is_feas), W(is_feas), '-', 'Color', skyBlue, 'LineWidth', 1.5);
plot(ax1, SNRdB(is_feas), W(is_feas), 'o', 'Color', skyBlue, 'MarkerFaceColor', skyBlue, 'MarkerEdgeColor', 'w', 'MarkerSize', 3.5);
h2 = plot(ax1, SNRdB(~is_feas), zeros(1, nnz(~is_feas)), 'x', 'Color', orange, 'LineWidth', 1.2, 'MarkerSize', 5);
if isfinite(W_inf_base)
    h3 = yline(ax1, W_inf_base, '--', 'Color', black, 'LineWidth', 1.2);
    text(ax1, max(SNRdB)-1, W_inf_base + 0.02, sprintf('$W_{\\infty}=%.3f$', W_inf_base), ...
        'HorizontalAlignment', 'right');
else
    h3 = [];
end
if ~isempty(first_feas)
    xline(ax1, SNRdB(first_feas), ':', 'Color', gray, 'LineWidth', 1.0);
end
xlabel(ax1, 'SNR (dB)');
ylabel(ax1, 'Feasible width $W(P)$');
title(ax1, '(a) Feasible-interval width vs SNR');
xlim(ax1, [SNRdB(1), SNRdB(end)]);
ylim(ax1, [0, min(1, max(W(is_feas)) + 0.05)]);
if isfinite(W_inf_base)
    legend(ax1, [h1 h2 h3], {'Analytical $W(P)$', 'Infeasible SNR', 'Asymptotic limit'}, 'Location', 'southeast', 'Box', 'off');
else
    legend(ax1, [h1 h2], {'Analytical $W(P)$', 'Infeasible SNR'}, 'Location', 'southeast', 'Box', 'off');
end
% -----------------------------------------------------------------------
% Panel (b): Asymptotic Outage / Bifurcation
% -----------------------------------------------------------------------
ax2 = nexttile(2);
hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');
set(ax2, 'LineWidth', 0.8, 'GridAlpha', 0.15);
hc  = plot(ax2, alpha_grid, pc_inf, '-', 'Color', black, 'LineWidth', 1.5);
hp0 = plot(ax2, alpha_grid, pp_inf(1,:), '--o', ...
    'Color', skyBlue, ...
    'LineWidth', 1.5, ...
    'MarkerSize', 4, ...
    'MarkerFaceColor', 'w', ...
    'MarkerIndices', 1:400:length(alpha_grid)); % <--- The magic trick
hp1 = plot(ax2, alpha_grid, pp_inf(2,:), '-.', 'Color', orange, 'LineWidth', 1.5);
% Target line
yline(ax2, par.eps_p_u, ':', 'Color', gray, 'LineWidth', 1.0);
text(ax2, 0.02, par.eps_p_u + 0.06, '${\epsilon}_c^u=\epsilon_p^u=0.10$');
% =========================================================================
% Vertical guides & Line-carried Labels
% =========================================================================
% 1. Alpha Min Vertical Line and Label
xline(ax2, alpha_min_inf_bifurc, ':', 'Color', gray, 'LineWidth', 1.5, 'HandleVisibility', 'off');

text(ax2, alpha_min_inf_bifurc - 0.04, 0.5, '$\alpha_{\min,\infty}^{\mathrm{ana}}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 10, ...
    'Rotation', 90, ...
    'HorizontalAlignment', 'center', ...
    'Color', black);

% 2. Alpha Max Vertical Line and Label (if it exists)
if ~isnan(alpha_max_inf_bifurc(2))
    xline(ax2, alpha_max_inf_bifurc(2), ':', 'Color', orange, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Removed sprintf() and added the missing closing brace '}'
    text(ax2, alpha_max_inf_bifurc(2) - 0.04, 0.5, '$\alpha_{\max,\infty}^{\mathrm{ana}}$', ...
        'Interpreter', 'latex', ...
        'FontSize', 10, ...
        'Rotation', 90, ...
        'HorizontalAlignment', 'center', ...
        'Color', orange);
end
xlabel(ax2, 'Common power fraction $\alpha$');
ylabel(ax2, 'Asymptotic outage prob.');
title(ax2, '(b) Perfect-SIC vs imperfect-SIC bifurcation');
% Adjusted legend location so it doesn't block the curves in the smaller frame
lgd = legend(ax2, [hc, hp0, hp1], ...
    {'Com. outage', 'Pri. outage ($\delta=0$)', 'Pri. outage ($\delta=0.05$)'}, ...
    'Box', 'off');
lgd.Location = 'northeast';
lgd.FontSize = 6;
drawnow; 
pos = lgd.Position;
nudge_amount = 0.28; 
pos(1) = pos(1) - nudge_amount; 
lgd.Position = pos;
xlim(ax2, [0 1]);
ylim(ax2, [0 1.0]);
% -----------------------------------------------------------------------
% Export
% -----------------------------------------------------------------------
fprintf('Exporting IEEE-ready Figure...\n');
pdfFileName = 'fig4_high_snr_frontier.pdf';
exportgraphics(fig, pdfFileName, 'ContentType', 'vector', 'BackgroundColor', 'w');
fprintf('Successfully exported %s at %.1f x %.1f inches.\n', pdfFileName, figWidth, figHeight);

%% ================== 8. Diagnostic Output to Command Window ===============
fprintf('\n================ DIAGNOSTIC DATA ================\n');

fprintf('\n--- Underlying Endpoint Values & Finite-SNR Width ---\n');
fprintf(' SNR (dB) | alpha_min(P) | alpha_max(P) | W(P)\n');
fprintf('---------------------------------------------------\n');
for i = 1:numel(SNRdB)
    fprintf(' %8d | %12.4f | %12.4f | %12.4f\n', SNRdB(i), alpha_min(i), alpha_max(i), W(i));
end

fprintf('\n--- Asymptotic Bifurcation Data ---\n');
fprintf('Asymptotic alpha_min_inf: %.4f\n', alpha_min_inf_bifurc);
fprintf('Perfect SIC (delta=%.2f): alpha_max_inf = %.4f\n', delta_list(1), alpha_max_inf_bifurc(1));
fprintf('Imperfect SIC (delta=%.2f): alpha_max_inf = %.4f\n', delta_list(2), alpha_max_inf_bifurc(2));

W_inf_perfect = alpha_max_inf_bifurc(1) - alpha_min_inf_bifurc;
W_inf_imperfect = alpha_max_inf_bifurc(2) - alpha_min_inf_bifurc;

fprintf('\nAsymptotic Width W_inf (Perfect SIC): %.4f\n', W_inf_perfect);
fprintf('Asymptotic Width W_inf (Imperfect SIC): %.4f\n', W_inf_imperfect);

if W_inf_imperfect < 1 - 1e-4
    fprintf('Does W_inf stay bounded away from 1? YES (W_inf = %.4f < 1)\n', W_inf_imperfect);
else
    fprintf('Does W_inf stay bounded away from 1? NO\n');
end

fprintf('\n--- The One-Sentence Message ---\n');
fprintf('"High-SNR does not remove the obstruction because residual interference scales linearly with the transmit power, introducing a hard asymptotic ceiling on alpha_max that permanently bounds the feasible interval W_inf away from 1."\n');
fprintf('=================================================\n\n');

%% ========================== 7. Local Functions ===========================
function pc = common_outage_approx(alpha, P, par)
    if alpha <= 0, pc = 1.0; return; end
    if alpha > 1, pc = NaN; return; end
    sc = par.Tc / (alpha * P * par.sigma_h2);
    succ = exp(-sc * par.sigma_n2) * (1 + par.Tc * (1 - alpha) * par.tau^2 / (par.K * alpha))^(par.M - par.K) * (1 + par.Tc * (1 - alpha) / (par.K * alpha))^(-par.M);
    pc = 1 - min(max(real(succ), 0), 1);
end
function pp = private_outage_approx_stable(alpha, P, par)
    if alpha >= 1, pp = 1.0; return; end
    if alpha < 0, pp = NaN; return; end
    theta = (1 - alpha) * P * par.sigma_h2 / par.K;
    if theta <= 0, pp = 1.0; return; end
    theta2 = par.tau^2 * theta;
    thetaRes = par.delta * alpha * P * par.sigma_h2;
    sp = par.Tp / theta;
    piw = mixture_weights_stable(par.M, par.tau);
    logTermsN = -inf(par.M,1);
    for n = 1:par.M
        logTermsI = -inf(n,1);
        for i = 0:(n-1)
            logTermsEll = -inf(i+1,1);
            for ell = 0:i
                logTermsT = -inf(ell+1,1);
                for t = 0:ell
                    logA = logPow(theta2, t) + logPow(thetaRes, ell - t);
                    if isinf(logA), continue; end
                    logCoeff = gammaln(ell + 1) - gammaln(t + 1) + gammaln(par.m + t) - gammaln(par.m);
                    logDen = (par.m + t) * log1p(theta2 * sp) + (1 + ell - t) * log1p(thetaRes * sp);
                    logTermsT(t + 1) = logCoeff + logA - logDen;
                end
                logTell = logsumexp(logTermsT);
                if isinf(logTell), continue; end
                logBinom = gammaln(i + 1) - gammaln(ell + 1) - gammaln(i - ell + 1);
                logNoise = logPow(par.sigma_n2, i - ell);
                logTermsEll(ell + 1) = logBinom + logNoise + logTell;
            end
            logAccL = logsumexp(logTermsEll);
            if isinf(logAccL), continue; end
            if i == 0, logPref = -sp * par.sigma_n2; else, logPref = -sp * par.sigma_n2 + i * log(sp) - gammaln(i + 1); end
            logTermsI(i + 1) = logPref + logAccL;
        end
        logTermsN(n) = log(piw(n)) + logsumexp(logTermsI);
    end
    succ = exp(logsumexp(logTermsN));
    pp = 1 - min(max(real(succ), 0), 1);
end
function pis = mixture_weights_stable(M, tau)
    if M == 1, pis = 1; return; end
    q = tau^2;
    if q == 0, pis = zeros(M,1); pis(M) = 1; return; end
    if q == 1, pis = zeros(M,1); pis(1) = 1; return; end
    pis = zeros(M,1); pis(1) = q^(M-1);
    ratio = (1 - q) / q;
    for n = 1:(M-1), pis(n+1) = pis(n) * ((M - n) / n) * ratio; end
    pis = pis / sum(pis);
end
function aMin = find_alpha_min(P, par)
    f = @(a) common_outage_approx(a, P, par) - par.eps_c_bar;
    if f(1.0) > 0, aMin = NaN; return; end
    if f(par.alpha_eps) <= 0, aMin = 0.0; return; end
    lo = par.alpha_eps; hi = 1.0;
    for it = 1:100, mid = 0.5 * (lo + hi); if f(mid) > 0, lo = mid; else, hi = mid; end; end
    aMin = hi;
end
function aMax = find_alpha_max(P, par)
    f = @(a) private_outage_approx_stable(a, P, par) - par.eps_p_u;
    if f(0.0) > 0, aMax = NaN; return; end
    if f(1 - par.alpha_eps) <= 0, aMax = 1.0; return; end
    lo = 0.0; hi = 1 - par.alpha_eps;
    for it = 1:100, mid = 0.5 * (lo + hi); if f(mid) <= 0, lo = mid; else, hi = mid; end; end
    aMax = lo;
end
function psi = psi_c_vec(alpha, par)
    psi = (1 + par.Tc * (1 - alpha) * par.tau^2 ./ (par.K * alpha)).^(par.M - par.K) ./ (1 + par.Tc * (1 - alpha) ./ (par.K * alpha)).^par.M;
end
function Psi = psi_p_vec(alpha, delta, par, pi_n)
    a = par.Tp * par.tau^2;
    Psi = zeros(size(alpha));
    for ia = 1:numel(alpha)
        b = par.Tp * delta * par.K * alpha(ia) / (1 - alpha(ia));
        acc = 0;
        for n = 1:par.M
            sum_n = 0;
            for ell = 0:(n - 1)
                for t = 0:ell
                    rising = exp(gammaln(par.m + t) - gammaln(par.m));
                    term = (rising / factorial(t)) * pow0(a, t) * pow0(b, ell - t) / ((1 + a)^(par.m + t) * (1 + b)^(1 + ell - t));
                    sum_n = sum_n + term;
                end
            end
            acc = acc + pi_n(n) * sum_n;
        end
        Psi(ia) = acc;
    end
    Psi = min(max(real(Psi), 0), 1);
end
function val = psi_c(alpha, par)
    if alpha <= 0 || alpha >= 1, val = NaN; return; end
    val = (1 + par.Tc * (1 - alpha) * par.tau^2 / (par.K * alpha))^(par.M - par.K) * (1 + par.Tc * (1 - alpha) / (par.K * alpha))^(-par.M);
end
function val = psi_p(alpha, par)
    if alpha <= 0 || alpha >= 1, val = NaN; return; end
    pis = mixture_weights_stable(par.M, par.tau);
    a = par.Tp * par.tau^2;
    b = par.Tp * par.delta * par.K * alpha / (1 - alpha);
    val = 0.0;
    for n = 1:par.M
        inner = 0.0;
        for ell = 0:(n - 1)
            for t = 0:ell
                rising = exp(gammaln(par.m + t) - gammaln(par.m));
                inner = inner + (rising / factorial(t)) * a^t * b^(ell - t) / ((1 + a)^(par.m + t) * (1 + b)^(1 + ell - t));
            end
        end
        val = val + pis(n) * inner;
    end
end
function aMinInf = find_alpha_min_inf(par)
    qc = 1 - par.eps_c_bar;
    f = @(a) psi_c(a, par) - qc;
    if f(1 - par.alpha_eps) < 0, aMinInf = NaN; return; end
    if f(par.alpha_eps) >= 0, aMinInf = 0.0; return; end
    lo = par.alpha_eps; hi = 1 - par.alpha_eps;
    for it = 1:120, mid = 0.5 * (lo + hi); if f(mid) < 0, lo = mid; else, hi = mid; end; end
    aMinInf = hi;
end
function aMaxInf = find_alpha_max_inf(par)
    qp = 1 - par.eps_p_u;
    f = @(a) psi_p(a, par) - qp;
    if f(par.alpha_eps) < 0, aMaxInf = NaN; return; end
    if f(1 - par.alpha_eps) >= 0, aMaxInf = 1.0; return; end
    lo = par.alpha_eps; hi = 1 - par.alpha_eps;
    for it = 1:120, mid = 0.5 * (lo + hi); if f(mid) >= 0, lo = mid; else, hi = mid; end; end
    aMaxInf = lo;
end
function a_star = threshold_from_increasing(x, y, target)
    idx = find(y >= target, 1, 'first');
    if isempty(idx), a_star = NaN; elseif idx == 1, a_star = x(1); else
        x1 = x(idx - 1); x2 = x(idx); y1 = y(idx - 1); y2 = y(idx);
        if y2 == y1, a_star = x2; else, a_star = x1 + (target - y1) * (x2 - x1) / (y2 - y1); end
    end
end
function a_star = threshold_from_decreasing(x, y, target)
    idx = find(y >= target, 1, 'last');
    if isempty(idx), a_star = NaN; elseif idx == numel(x), a_star = x(end); else
        x1 = x(idx); x2 = x(idx + 1); y1 = y(idx); y2 = y(idx + 1);
        if y2 == y1, a_star = x1; else, a_star = x1 + (target - y1) * (x2 - x1) / (y2 - y1); end
    end
end
function y = pow0(a, p)
    if p == 0, y = 1; else, y = a^p; end
end
function y = logPow(a, p)
    if p == 0, y = 0; elseif a == 0, y = -inf; else, y = p * log(a); end
end
function s = logsumexp(v)
    vmax = max(v);
    if isinf(vmax), s = -inf; else, s = vmax + log(sum(exp(v - vmax))); end
end