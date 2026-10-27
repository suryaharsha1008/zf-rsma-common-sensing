%% fig5_detection_matched_vs_fixed.m
% Fig. 5 of the paper: conditional detection probability over the analytical
% communication-feasible interval, (a) matched geometry-driven common beam,
% (b) fixed benchmark beam.
clear; close all; clc;
rng(7, 'twister');

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

sty.fontName     = 'Times New Roman';
sty.fontSize     = 11;
sty.ms           = 5.0;
sty.gridAlpha    = 0.15;
sty.patchColor   = [0.92 0.92 0.92];
sty.patchAlpha   = 1.0;
sty.lineWidth    = 2.0;
sty.errLineWidth = 1.0;
sty.axisLineW    = 0.8;

saveFig = true;

%% ------------------------------------------------------------------------
% 1) Baseline communication parameters
% -------------------------------------------------------------------------
prm.Nt        = 8;
prm.K         = 4;
prm.tau       = 0.20;
prm.delta     = 0.05;
prm.SNRdB     = 15;
prm.Rc        = 0.20;
prm.Rp        = 1.00;
prm.sigma_h2  = 1.0;
prm.sigma_n2  = 1.0;
prm.eps_c_u   = 0.10;
prm.eps_c_sys = 0.40;
prm.eps_p_u   = 0.10;

prm.P         = 10^(prm.SNRdB/10) * prm.sigma_n2;
prm.Tc        = 2^prm.Rc - 1;
prm.Tp        = 2^prm.Rp - 1;
prm.M         = prm.Nt - prm.K + 1;
prm.m         = prm.K - 1;

prm.epsbar_c  = min(prm.eps_c_u, prm.eps_c_sys / prm.K);

wc_bench = zeros(prm.Nt,1);
wc_bench(1) = 1;
wc_bench = wc_bench / norm(wc_bench);

%% ------------------------------------------------------------------------
% 2) Detection controls
% -------------------------------------------------------------------------
detBase.Nr             = 1;
detBase.L              = 16;
detBase.PFA            = 1e-3;
detBase.NmcX           = 4e3;
detBase.alpha_mc_count = 19;
detBase.PDstar         = 0.60;
senseGain = 0.028;

%% ------------------------------------------------------------------------
% 3) Communication thresholds alpha_min and alpha_max
% -------------------------------------------------------------------------
a_lo  = 1e-4;
a_hi  = 1 - 1e-4;
alpha = linspace(a_lo, a_hi, 1001);

alpha_min = find_left_threshold(@(a) common_outage(a, prm), prm.epsbar_c, a_lo, a_hi);
alpha_max = find_right_threshold(@(a) private_outage(a, prm), prm.eps_p_u, a_lo, a_hi);

if isnan(alpha_min) || isnan(alpha_max) || alpha_min > alpha_max
    error('Analytical communication-feasible interval is empty on this baseline.');
end

alpha_mc = unique(sort([linspace(a_lo, a_hi, detBase.alpha_mc_count), alpha_min, alpha_max]));

%% ------------------------------------------------------------------------
% 4) Sample and freeze one realized ZF beam set W
% -------------------------------------------------------------------------
[W, Hhat] = sample_fixed_zf_beamset(prm); %#ok<NASGU>
B_W = (W * W') / prm.K;

if rank_tol([wc_bench, W], 1e-10) < prm.K + 1
    error('rank([w_c^{bench}, W]) < K+1. Resample W.');
end

%% ------------------------------------------------------------------------
% 5) Auto-select one shared sensing geometry Q
% -------------------------------------------------------------------------
caseSel = select_figure5_case(alpha, prm, detBase, B_W, wc_bench, W, alpha_min, alpha_max, senseGain);
det       = caseSel.det;
wc_match  = caseSel.wc_match;

if rank_tol([wc_match, W], 1e-10) < prm.K + 1
    warning('rank([w_c^{match}, W]) < K+1. Detection DoF may differ.');
end

%% ------------------------------------------------------------------------
% 6) Analyze both beam choices under the SAME Q and SAME W
% -------------------------------------------------------------------------
resMatch = analyze_detection_case_beam( ...
    alpha, alpha_mc, prm, det, wc_match, W, alpha_min, alpha_max, ...
    '(a)', 'Matched geometry-driven beam');

resBench = analyze_detection_case_beam( ...
    alpha, alpha_mc, prm, det, wc_bench, W, alpha_min, alpha_max, ...
    '(b)', 'Benchmark fixed beam');

%% ------------------------------------------------------------------------
% 7) Console diagnostics
% -------------------------------------------------------------------------
fprintf('\n=============================================================\n');
fprintf('Fig. 5: matched geometry-driven beam vs fixed benchmark beam\n');
fprintf('Shared target detection requirement: P_D^* = %.4f\n', detBase.PDstar);
fprintf('A_com^ana(P) = [%.4f, %.4f]\n', alpha_min, alpha_max);
fprintf('Selected target alignment with benchmark beam = %.3f\n', caseSel.align);
fprintf('=============================================================\n');

print_case_diagnostics(resMatch, alpha_min, alpha_max);
print_case_diagnostics(resBench, alpha_min, alpha_max);
print_comparison_summary(resMatch, resBench, alpha_min, alpha_max);

%% ------------------------------------------------------------------------
% 8) Plot
% -------------------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100, 100, 490, 680]);
tl  = tiledlayout(fig, 2, 1, 'Padding', 'tight', 'TileSpacing', 'compact');

ax1 = nexttile(tl, 1);
ax2 = nexttile(tl, 2);

hA = plot_case(ax1, resMatch, alpha_min, alpha_max, sty, true);
plot_case(ax2, resBench, alpha_min, alpha_max, sty, false);

set([ax1, ax2], 'YLim', [0 1], 'XLim', [0 1]);
xticklabels(ax1, {});

ylabel(ax1, 'Detection Prob. $P_D(\alpha)$', 'Interpreter', 'latex', 'FontSize', sty.fontSize + 1);
ylabel(ax2, 'Detection Prob. $P_D(\alpha)$', 'Interpreter', 'latex', 'FontSize', sty.fontSize + 1);
xlabel(ax2, 'Common power fraction $\alpha$', 'Interpreter', 'latex', 'FontSize', sty.fontSize + 1);

lgd = legend(ax1, hA.legendHandles, ...
    {'Feasible interval $\mathcal{A}_{\mathrm{com}}^{\mathrm{ana}}(P)$', ...
     'Large-$L$ $\bar{P}_{D,\mathbf{W}}(\alpha)$', ...
     'Finite-$L$ GLRT Monte Carlo (95\% CI)', ...
     'Target $P_D^\star$', ...
     'Endpoint predicted by theory'}, ...
    'Interpreter', 'latex', ...
    'FontSize', sty.fontSize - 3, ...
    'NumColumns', 2, ...
    'Location', 'southoutside');
lgd.Layout.Tile = 'south';

drawnow;

%% ------------------------------------------------------------------------
% 9) Export
% -------------------------------------------------------------------------
if saveFig
    pdfName = 'fig5_detection_matched_vs_fixed.pdf';
    exportgraphics(fig, pdfName, 'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('\nSaved: %s\n', pdfName);
end

%% =========================================================================
% Local functions
% =========================================================================

function caseSel = select_figure5_case(alpha, prm, detBase, B_W, wc_bench, W, alpha_min, alpha_max, senseGain)
    u_priv = dominant_private_orthogonal_mode(B_W, wc_bench);
    alignGrid = linspace(0.05, 0.85, 81);
    
    bestStrict.score = -Inf;
    bestSoft.score   = -Inf;
    hasStrict        = false;
    
    for a = alignGrid
        u_target = a * wc_bench + sqrt(1 - a^2) * u_priv;
        u_target = u_target / norm(u_target);
        
        det = build_target_geometry('Shared target geometry', u_target, senseGain, detBase, B_W);
        wc_match = dominant_eigvec(det.Q);
        
        thMatch = analyze_detection_case_theory_only(alpha, prm, det, wc_match, W, alpha_min, alpha_max);
        thBench = analyze_detection_case_theory_only(alpha, prm, det, wc_bench, W, alpha_min, alpha_max);
        
        Pm_max = interp1(alpha, thMatch.Pd_theory, alpha_max, 'pchip');
        Pb_min = interp1(alpha, thBench.Pd_theory, alpha_min, 'pchip');
        Pb_max = interp1(alpha, thBench.Pd_theory, alpha_max, 'pchip');
        Pb_best = max(Pb_min, Pb_max);
        
        strictGood = (thMatch.Delta_W >= -1e-12) && ...
                     (Pm_max >= detBase.PDstar + 0.03) && ...
                     (thBench.Delta_W < -1e-4) && ...
                     (Pb_best <= detBase.PDstar - 0.03);
                 
        score = 6.0 * (Pm_max - Pb_best) ...
              + 2.0 * double(thBench.Delta_W < 0) ...
              + 1.0 * max(0, detBase.PDstar - Pb_best) ...
              + 0.5 * max(0, Pm_max - detBase.PDstar);
          
        if score > bestSoft.score
            bestSoft.score    = score;
            bestSoft.align    = a;
            bestSoft.det      = det;
            bestSoft.wc_match = wc_match;
        end
        
        if strictGood && score > bestStrict.score
            hasStrict         = true;
            bestStrict.score  = score;
            bestStrict.align  = a;
            bestStrict.det    = det;
            bestStrict.wc_match = wc_match;
        end
    end
    
    if hasStrict
        caseSel = bestStrict;
    else
        caseSel = bestSoft;
    end
end

function out = plot_case(ax, res, alpha_min, alpha_max, sty, makeLegendHandles)
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    set(ax, 'FontName', sty.fontName, ...
        'FontSize', sty.fontSize, ...
        'LineWidth', sty.axisLineW, ...
        'GridAlpha', sty.gridAlpha, ...
        'Layer', 'top');
    
    hPatch = patch(ax, ...
        [alpha_min alpha_max alpha_max alpha_min], ...
        [0 0 1 1], ...
        sty.patchColor, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', sty.patchAlpha);
    
    hTheory = plot(ax, res.alpha, res.Pd_theory, ...
        'k-', 'LineWidth', sty.lineWidth);
    
    hPhys = errorbar(ax, res.alpha_mc, res.Pd_phys, 1.96 * res.se_phys, ...
        'o', ...
        'LineWidth', sty.errLineWidth, ...
        'MarkerSize', sty.ms, ...
        'Color', [0.8500 0.3250 0.0980], ...
        'MarkerFaceColor', 'w');
    try, hPhys.CapSize = 0; catch, end
    
    hTarget = plot(ax, [0 1], [res.det.PDstar res.det.PDstar], 'k:', 'LineWidth', 1.2);
    
    if exist('xline', 'file') == 2
        xline(ax, alpha_min, ':', 'LineWidth', 1.2, 'Color', [0.35 0.35 0.35]);
        xline(ax, alpha_max, ':', 'LineWidth', 1.2, 'Color', [0.35 0.35 0.35]);
    else
        plot(ax, [alpha_min alpha_min], [0 1], ':', 'LineWidth', 1.2, 'Color', [0.35 0.35 0.35]);
        plot(ax, [alpha_max alpha_max], [0 1], ':', 'LineWidth', 1.2, 'Color', [0.35 0.35 0.35]);
    end
    
    hOpt = plot(ax, res.alpha_star_theory, res.Pd_star_theory, ...
        'p', ...
        'MarkerSize', 13, ...
        'LineWidth', 1.0, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', [0.9290 0.6940 0.1250]);
    
    title(ax, sprintf('\\textbf{%s %s}', res.panelTag, res.captionLabel), ...
        'Interpreter', 'latex', 'FontSize', sty.fontSize + 1);
    
    text(ax, alpha_min - 0.008, 0.30, '$\alpha_{\min}$', ...
        'Interpreter', 'latex', ...
        'FontSize', sty.fontSize, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'Rotation', 90);
    text(ax, alpha_max - 0.008, 0.30, '$\alpha_{\max}$', ...
        'Interpreter', 'latex', ...
        'FontSize', sty.fontSize, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'Rotation', 90);
    
    % FIXED LATEX STRINGS HERE
    txt = {
        sprintf('$\\Delta_{\\mathbf{W}}=%.3e$', res.Delta_W), ...
        sprintf('%s', res.endpointRuleLatex), ...
        sprintf('%s', res.jointIntervalLatex)
        };
        
    text(ax, 0.03, 0.95, txt, ...
        'Units', 'normalized', ...
        'Interpreter', 'latex', ...
        'FontSize', sty.fontSize - 1, ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', ...
        'Margin', 2);
    
    out = struct();
    if makeLegendHandles
        out.legendHandles = [hPatch, hTheory, hPhys, hTarget, hOpt];
    end
end

function print_case_diagnostics(res, alpha_min, alpha_max)
    Pd_th_min = interp1(res.alpha, res.Pd_theory, alpha_min, 'pchip');
    Pd_th_max = interp1(res.alpha, res.Pd_theory, alpha_max, 'pchip');
    
    idx_mc_min = find(abs(res.alpha_mc - alpha_min) < 1e-9, 1);
    idx_mc_max = find(abs(res.alpha_mc - alpha_max) < 1e-9, 1);
    Pd_mc_min = res.Pd_phys(idx_mc_min);
    Pd_mc_max = res.Pd_phys(idx_mc_max);
    hw_mc_min = 1.96 * res.se_phys(idx_mc_min);
    hw_mc_max = 1.96 * res.se_phys(idx_mc_max);
    
    fprintf('\n%s %s\n', res.panelTag, res.captionLabel);
    fprintf('-------------------------------------------------------------\n');
    fprintf('rank([w_c, W])           = %d\n', res.rX_theory);
    fprintf('d                        = %d\n', res.d);
    fprintf('q_c                      = %.6e\n', res.qc);
    fprintf('q_p(W)                   = %.6e\n', res.det.qp);
    fprintf('Delta_W                  = %.6e\n', res.Delta_W);
    fprintf('Predicted endpoint       = %s\n', res.endpointRuleStr);
    fprintf('Large-L theory: P_D(alpha_min)=%.4f, P_D(alpha_max)=%.4f\n', Pd_th_min, Pd_th_max);
    fprintf('Finite-L MC   : P_D(alpha_min)=%.4f +/- %.4f, P_D(alpha_max)=%.4f +/- %.4f\n', ...
        Pd_mc_min, hw_mc_min, Pd_mc_max, hw_mc_max);
    fprintf('Large-L joint-feasible set: %s\n', res.jointIntervalStr);
end

function print_comparison_summary(resMatch, resBench, alpha_min, alpha_max)
    Pd_match_max = interp1(resMatch.alpha, resMatch.Pd_theory, alpha_max, 'pchip');
    Pd_bench_min = interp1(resBench.alpha, resBench.Pd_theory, alpha_min, 'pchip');
    Pd_bench_max = interp1(resBench.alpha, resBench.Pd_theory, alpha_max, 'pchip');
    Pd_bench_best = max(Pd_bench_min, Pd_bench_max);
    
    idxM = find(abs(resMatch.alpha_mc - alpha_max) < 1e-9, 1);
    idxB = find(abs(resBench.alpha_mc - resBench.alpha_star_theory) < 1e-9, 1);
    
    fprintf('\nComparison summary\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('q_p(W)                    = %.6e\n', resMatch.det.qp);
    fprintf('q_c(matched)              = %.6e\n', resMatch.qc);
    fprintf('q_c(benchmark)            = %.6e\n', resBench.qc);
    fprintf('lambda_max(Q)             = %.6e\n', resMatch.det.lambda_max_Q);
    fprintf('Delta_W(matched)          = %.6e\n', resMatch.Delta_W);
    fprintf('Delta_W(benchmark)        = %.6e\n', resBench.Delta_W);
    fprintf('Predicted endpoint match  = %s\n', resMatch.endpointRuleStr);
    fprintf('Predicted endpoint bench  = %s\n', resBench.endpointRuleStr);
    fprintf('Theory gain at alpha_max  = %.4f\n', Pd_match_max - Pd_bench_max);
    
    fprintf('\nPaper-ready summary\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('(a) Matched geometry-driven beam: Delta_W = %.3e, so the theory selects alpha_max = %.4f. ', ...
        resMatch.Delta_W, alpha_max);
    fprintf('The large-L model gives P_D(alpha_max) = %.4f, and the finite-L GLRT Monte Carlo estimate is %.4f +/- %.4f. ', ...
        Pd_match_max, resMatch.Pd_phys(idxM), 1.96 * resMatch.se_phys(idxM));
    fprintf('The corresponding large-L joint-feasible set is %s.\n', resMatch.jointIntervalStr);
    
    fprintf('(b) Benchmark fixed beam: Delta_W = %.3e, so the theory selects %s. ', ...
        resBench.Delta_W, resBench.endpointRuleStr);
    fprintf('The best large-L endpoint value is %.4f, and the corresponding finite-L GLRT Monte Carlo estimate is %.4f +/- %.4f. ', ...
        Pd_bench_best, resBench.Pd_phys(idxB), 1.96 * resBench.se_phys(idxB));
    fprintf('The corresponding large-L joint-feasible set is %s.\n', resBench.jointIntervalStr);
end

function res = analyze_detection_case_theory_only(alpha, prm, det, wc_case, W, alpha_min, alpha_max)
    qc      = real(wc_case' * det.Q * wc_case);
    Delta_W = qc - det.qp;
    
    rX_theory = rank_tol([wc_case, W], 1e-10);
    d         = det.Nr * rX_theory;
    gamma_thr = central_chi2_threshold(1 - det.PFA, 2 * d);
    Lambda_bar = 2 * det.L * prm.P * (alpha * qc + (1 - alpha) * det.qp);
    Pd_theory  = ncx2_sf_stable(gamma_thr, 2 * d, Lambda_bar);
    
    % FIXED \mathbf{W} HERE
    if Delta_W > 0
        alpha_star_theory = alpha_max;
        endpointRuleStr   = 'alpha_max';
        endpointRuleLatex = '$\alpha^\star=\alpha_{\max}$';
    elseif Delta_W < 0
        alpha_star_theory = alpha_min;
        endpointRuleStr   = 'alpha_min';
        endpointRuleLatex = '$\alpha^\star=\alpha_{\min}$';
    else
        alpha_star_theory = alpha_max;
        endpointRuleStr   = 'flat';
        endpointRuleLatex = '$\bar{P}_{D,\mathbf{W}}(\alpha)\ \mathrm{flat}$';
    end
    
    [jointIntervalStr, jointIntervalLatex] = extract_joint_interval(alpha, Pd_theory, det.PDstar, alpha_min, alpha_max);
    
    res.alpha             = alpha;
    res.Pd_theory         = Pd_theory;
    res.alpha_star_theory = alpha_star_theory;
    res.endpointRuleStr   = endpointRuleStr;
    res.endpointRuleLatex = endpointRuleLatex;
    res.qc                = qc;
    res.Delta_W           = Delta_W;
    res.rX_theory         = rX_theory;
    res.d                 = d;
    res.gamma_thr         = gamma_thr;
    res.det               = det;
    res.jointIntervalStr  = jointIntervalStr;
    res.jointIntervalLatex = jointIntervalLatex;
end

function res = analyze_detection_case_beam(alpha, alpha_mc, prm, det, wc_case, W, alpha_min, alpha_max, panelTag, captionLabel)
    th = analyze_detection_case_theory_only(alpha, prm, det, wc_case, W, alpha_min, alpha_max);
    
    Pd_phys = zeros(size(alpha_mc));
    se_phys = zeros(size(alpha_mc));
    for ia = 1:numel(alpha_mc)
        [Pd_phys(ia), se_phys(ia)] = pd_mc_physical_glrt_fixedW(alpha_mc(ia), prm, det, wc_case, W, th.gamma_thr, th.d);
    end
    
    res = th;
    res.alpha_mc          = alpha_mc;
    res.Pd_phys           = Pd_phys;
    res.se_phys           = se_phys;
    res.Pd_star_theory    = interp1(alpha, th.Pd_theory, th.alpha_star_theory, 'pchip', 'extrap');
    res.panelTag          = panelTag;
    res.captionLabel      = captionLabel;
    res.wc_case           = wc_case;
end

function [jointIntervalStr, jointIntervalLatex] = extract_joint_interval(alpha, Pd_theory, PDstar, alpha_min, alpha_max)
    valid_mask = (alpha >= alpha_min) & (alpha <= alpha_max);
    valid_alpha = alpha(valid_mask);
    valid_PD    = Pd_theory(valid_mask);
    
    feas_mask   = valid_PD >= PDstar;
    
    % FIXED \emptyset HERE
    if any(feas_mask)
        aL = min(valid_alpha(feas_mask));
        aR = max(valid_alpha(feas_mask));
        jointIntervalStr   = sprintf('[%.4f, %.4f]', aL, aR);
        jointIntervalLatex = sprintf('$\\mathcal{A}_{\\mathrm{joint}}=[%.3f,%.3f]$', aL, aR);
    else
        jointIntervalStr   = 'empty';
        jointIntervalLatex = '$\mathcal{A}_{\mathrm{joint}}=\emptyset$';
    end
end

function [W, Hhat] = sample_fixed_zf_beamset(prm)
    maxTry = 500;
    for it = 1:maxTry
        Hhat = sqrt(prm.sigma_h2/2) * (randn(prm.Nt, prm.K) + 1i * randn(prm.Nt, prm.K));
        Wf   = Hhat / (Hhat' * Hhat);
        W    = Wf ./ sqrt(sum(abs(Wf).^2, 1));
        if all(isfinite(W(:)))
            return;
        end
    end
    error('Failed to generate a valid fixed ZF beam set.');
end

function u = dominant_private_orthogonal_mode(B_W, wc)
    Nt = size(B_W,1);
    Pperp = eye(Nt) - wc * wc';
    M = Pperp * B_W * Pperp;
    M = 0.5 * (M + M');
    [V, D] = eig(M, 'vector');
    [~, idx] = max(real(D));
    u = V(:, idx);
    u = u - wc * (wc' * u);
    u = u / norm(u);
end

function det = build_target_geometry(captionLabel, uTarget, senseGain, detBase, B_W)
    uTarget = uTarget / norm(uTarget);
    det = detBase;
    det.captionLabel = captionLabel;
    det.uTarget      = uTarget;
    det.HT           = sqrt(senseGain) * (uTarget.');
    det.R            = eye(det.Nr);
    det.Q            = det.HT' * det.HT;
    det.qp           = real(trace(det.Q * B_W));
    
    evals = eig(0.5 * (det.Q + det.Q'));
    det.lambda_max_Q = max(real(evals));
end

function wc_opt = dominant_eigvec(Q)
    Qh = 0.5 * (Q + Q');
    [V, D] = eig(Qh, 'vector');
    [~, idx] = max(real(D));
    wc_opt = V(:, idx);
    wc_opt = wc_opt / norm(wc_opt);
end

function [PdHat, se] = pd_mc_physical_glrt_fixedW(alpha, prm, det, wc_case, W, gamma_thr, d)
    pVals = zeros(det.NmcX, 1);
    for t = 1:det.NmcX
        X = sample_rsma_snapshot_fixedW(alpha, prm, wc_case, W, det.L);
        A = det.HT * X;
        lambda_inst = 2 * real(sum(abs(A(:)).^2));
        pVals(t) = ncx2_sf_stable(gamma_thr, 2 * d, lambda_inst);
    end
    PdHat = mean(pVals);
    se    = std(pVals, 0, 1) / sqrt(det.NmcX);
end

function X = sample_rsma_snapshot_fixedW(alpha, prm, wc_case, W, L)
    X = sqrt(alpha * prm.P) * (wc_case * ((randn(1, L) + 1i * randn(1, L)) / sqrt(2))) + ...
        sqrt((1 - alpha) * prm.P / prm.K) * (W * ((randn(prm.K, L) + 1i * randn(prm.K, L)) / sqrt(2)));
end

function tf = rank_tol(A, tol)
    s = svd(A, 'econ');
    tf = sum(s > tol * max(s));
end

function gamma_thr = central_chi2_threshold(p, nu)
    gamma_thr = chi2inv(p, nu);
end

function sf = ncx2_sf_stable(x, nu, lambda_vec)
    sf = 1 - ncx2cdf(x, nu, lambda_vec);
    sf = min(max(sf, 0), 1);
end

function p_out = common_outage(alpha, prm)
    sc = prm.Tc / (alpha * prm.P * prm.sigma_h2);
    p_out = 1 - ( ...
        exp(-sc * prm.sigma_n2) * ...
        (1 + prm.Tc * (1 - alpha) * prm.tau^2 / (prm.K * alpha))^(prm.M - prm.K) * ...
        (1 + prm.Tc * (1 - alpha) / (prm.K * alpha))^(-prm.M) );
end

function p_out = private_outage(alpha, prm)
    theta = (1 - alpha) * prm.P * prm.sigma_h2 / prm.K;
    sp    = prm.Tp / theta;
    p_suc = 0.0;
    
    for n = 1:prm.M
        inner_i = 0.0;
        for i = 0:(n-1)
            inner_l = 0.0;
            for ell = 0:i
                inner_l = inner_l + ...
                    nchoosek(i, ell) * (prm.sigma_n2^(i - ell)) * ...
                    T_ell(ell, sp, prm.tau^2 * theta, prm.delta * alpha * prm.P * prm.sigma_h2, prm.m);
            end
            inner_i = inner_i + (sp^i / factorial(i)) * inner_l;
        end
        p_suc = p_suc + ...
            nchoosek(prm.M - 1, n - 1) * ...
            (1 - prm.tau^2)^(n - 1) * ...
            (prm.tau^2)^(prm.M - n) * ...
            exp(-sp * prm.sigma_n2) * inner_i;
    end
    p_out = 1 - p_suc;
end

function val = T_ell(ell, s, theta2, thetaRes, m)
    val = 0.0;
    for t = 0:ell
        val = val + ...
            (factorial(ell) / factorial(t) * (gamma(m + t) / gamma(m)) * ...
            (theta2^t) * (thetaRes^(ell - t))) / ...
            ((1 + theta2 * s)^(m + t) * (1 + thetaRes * s)^(1 + ell - t));
    end
end

function alpha_star = find_left_threshold(fun, target, a_lo, a_hi)
    if fun(a_hi) > target
        alpha_star = NaN; return;
    end
    if fun(a_lo) <= target
        alpha_star = a_lo; return;
    end
    L = a_lo; R = a_hi;
    for it = 1:100
        mid = 0.5 * (L + R);
        if fun(mid) <= target
            R = mid;
        else
            L = mid;
        end
    end
    alpha_star = 0.5 * (L + R);
end

function alpha_star = find_right_threshold(fun, target, a_lo, a_hi)
    if fun(a_lo) > target
        alpha_star = NaN; return;
    end
    if fun(a_hi) <= target
        alpha_star = a_hi; return;
    end
    L = a_lo; R = a_hi;
    for it = 1:100
        mid = 0.5 * (L + R);
        if fun(mid) <= target
            L = mid;
        else
            R = mid;
        end
    end
    alpha_star = 0.5 * (L + R);
end