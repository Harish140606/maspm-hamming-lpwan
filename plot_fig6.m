%PLOT_FIG6 Reproduce Fig. 6: cross-correlations of different-length PSF chirps.
%
% PAPER REFERENCE
%   Section II-C, Fig. 6 ("For sufficiently different Li and Lj,
%   cross-correlations of g_hat_i[k] and g_hat_j[k] for constant-
%   envelope PSFs given by Eq. (11) will have large TBPs").
%
% EQUATIONS USED
%   Eq. (11): constant-envelope PSF chirp construction (via
%             pulse_shaping_filter.m).
%   Section II-C discussion: for sufficiently different Li, Lj, the
%   cross-correlation (g_hat_i * g_j)[k] has a large time-bandwidth
%   product (i.e. is spread out / low-peakedness), in contrast with
%   the compact same-PSF autocorrelation (small TBP).
%
% PARAMETERS USED
%   beta = cfg.beta_fig9 (1/4); two PSF lengths Li != Lj (illustrative,
%   Li=64, Lj=128); up-chirp and down-chirp directions.
%
% SIMULATION EXPLANATION
%   Top row: cross-correlation of an up-chirp PSF (length Li) with a
%   down-chirp matched filter (length Lj) -- reproduces the paper's
%   top panel of Fig. 6 (up-chirp * down-chirp / "own-chirp"), which
%   shows a broad, low-level, spread-out result of width ~|Lj-Li|.
%   Bottom row: cross-correlation of two up-chirp PSFs of different
%   lengths -- reproduces the paper's bottom panel, width ~Li+Lj.
%   Both are contrasted against the compact (small-TBP) SAME-length,
%   same-direction autocorrelation, shown for reference.

clear; clc;
cfg = config();
psf = pulse_shaping_filter();

beta = cfg.beta_fig9;
Ns = cfg.Ns_fig9;
Li = 64;
Lj = 128;

[zeta_i_up, ~, ~] = psf.design_chirp(Li, beta, Ns, +1);   % up-chirp, length Li
[zeta_j_down, ~, ~] = psf.design_chirp(Lj, beta, Ns, -1); % down-chirp, length Lj
[zeta_j_up, ~, ~] = psf.design_chirp(Lj, beta, Ns, +1);   % up-chirp, length Lj

% Reference: compact same-length, same-direction autocorrelation.
[zeta_ref, ~, ~] = psf.design_chirp(Li, beta, Ns, +1);
ac_ref = xcorr(zeta_ref, zeta_ref);
ac_ref = abs(ac_ref) / max(abs(ac_ref));

% Top panel: up-chirp (Li) vs down-chirp (Lj) cross-correlation.
cc_updown = xcorr(zeta_i_up, zeta_j_down);
cc_updown = abs(cc_updown) / max(abs(cc_updown));

% Bottom panel: up-chirp (Li) vs up-chirp (Lj) cross-correlation.
cc_upup = xcorr(zeta_i_up, zeta_j_up);
cc_upup = abs(cc_upup) / max(abs(cc_upup));

fig = figure('Position',[100 100 950 750], 'Color','w');

subplot(2,1,1); hold on; box on;
lags_ref = -(Li-1):(Li-1);
lags_ud  = -(Lj-1):(Lj-1);  % xcorr zero-pads to max(Li,Lj)=Lj
plot(lags_ref, ac_ref, 'k-', 'LineWidth', 1.2);
plot(lags_ud, cc_updown, 'r-', 'LineWidth', 1.4);
xlabel('lag k'); ylabel('|correlation| (normalized)');
title(sprintf('Fig. 6 (top) reproduction: up-chirp (L_i=%d) * down-chirp (L_j=%d)  --  large-TBP, spread width \\approx|L_j-L_i|=%d', Li, Lj, abs(Lj-Li)));
legend('same-PSF autocorrelation (compact, small TBP, reference)', ...
       'cross-correlation up*down (large TBP)', 'Location','best');
grid on;

subplot(2,1,2); hold on; box on;
lags_uu  = -(Lj-1):(Lj-1);  % xcorr zero-pads to max(Li,Lj)=Lj
plot(lags_ref, ac_ref, 'k-', 'LineWidth', 1.2);
plot(lags_uu, cc_upup, 'b-', 'LineWidth', 1.4);
xlabel('lag k'); ylabel('|correlation| (normalized)');
title(sprintf('Fig. 6 (bottom) reproduction: up-chirp (L_i=%d) * up-chirp (L_j=%d)  --  large-TBP, spread width \\approx L_i+L_j=%d', Li, Lj, Li+Lj));
legend('same-PSF autocorrelation (compact, small TBP, reference)', ...
       'cross-correlation up*up (large TBP)', 'Location','best');
grid on;

% Peakedness metric: fraction of samples (within each curve's own lag
% span) whose magnitude exceeds half the peak -- small for a compact
% (small-TBP) autocorrelation, large for a spread-out (large-TBP)
% cross-correlation. This quantifies the qualitative "spread vs
% compact" contrast the figure illustrates.
peakedness_ref    = mean(ac_ref    > 0.5*max(ac_ref));
peakedness_updown = mean(cc_updown > 0.5*max(cc_updown));
peakedness_upup   = mean(cc_upup   > 0.5*max(cc_upup));

if ~exist('results','dir'); mkdir('results'); end
saveas(fig, fullfile('results','fig6_cross_correlation_tbp.png'));
fprintf('Saved results/fig6_cross_correlation_tbp.png\n');
fprintf('Fraction of lags with |corr|>0.5*peak (small = compact/small-TBP, large = spread/large-TBP):\n');
fprintf('  same-PSF autocorrelation (reference): %.4f\n', peakedness_ref);
fprintf('  up*down cross-correlation:            %.4f\n', peakedness_updown);
fprintf('  up*up cross-correlation:               %.4f\n', peakedness_upup);
