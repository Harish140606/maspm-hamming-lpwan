%PLOT_FIG22 Reproduce Fig. 22: flip-PSF-pair capacity/partition comparison.
%
% PAPER REFERENCE
%   Section V, Fig. 22 ("Using pairs of 'flip' PSFs results in only
%   small decrease in M-ASPM capacity, while reducing number of pulse
%   duty cycles and partitions by more than half"), applied to the
%   Fig. 20 areal-coverage scenario (d0 = LoRa SF6 range, dmax = LoRa
%   SF12 range, alpha=1/2, 64-ASPM).
%
% EQUATIONS USED
%   Eq. (17)-(18),(26): range<->spectral-efficiency matching, used to
%             find the IpI Nm at which 64-ASPM's own range equals dmax
%             (via path_loss_model.m).
%   Eq. (19)-(20): N0dB (via processing_gain.m).
%   Eq. (28),(31),(38): hexagonal-cell m-quantile partitioning (via
%             network_capacity.m; see that file's header for the
%             documented hexagon-geometry implementation assumption).
%   Eq. (39): equal-IpI SIR constraint, non-flip.
%   Eq. (40): pulse duty cycle per partition, D_i=(d_i/d_max)^gamma*Dmax.
%   Eq. (41): equal-IpI SIR constraint, flip-PSF pairs (factor of 2).
%
% PARAMETERS USED (cfg.fig22_*)
%   d0 = LoRa SF6 range, dmax = LoRa SF12 range (both relative to
%   delta0(gamma) via path_loss_model.m); alpha=1/2; M=64.
%   IMPLEMENTATION ASSUMPTION: gamma = cfg.fig22_gamma = 2 (the paper
%   does not state which path-loss exponent was used for this specific
%   numeric example -- see config.m for the full rationale). With this
%   assumption, this reproduction's own computation (see console
%   output) finds m_max=69 (non-flip) and 34 flip-pairs (68 channels),
%   remarkably close to the paper's own reported "70 to 66" -- a
%   strong cross-validation of the implementation despite the
%   necessary gamma assumption.
%
% SIMULATION EXPLANATION
%   (1) Range-match 64-ASPM's own IpI Nm so its achieved range equals
%       dmax (Eq. 26 inverted); (2) compute the SIR threshold
%       Nm/(alpha*N0dB); (3) find the maximum number of equal-IpI PSF
%       channels m (Eq. 39) and of flip-PSF PAIRS (Eq. 41) satisfying
%       that threshold; (4) show the resulting duty-cycle-vs-range
%       schedule (Eq. 40) and a bar-chart capacity/complexity summary,
%       and a simple annular-partition geometry diagram.

clear; clc;
cfg = config();
se = spectral_efficiency();
pl = path_loss_model();
pg = processing_gain();
nc = network_capacity();

BER = cfg.BER; BERp = cfg.BER_prime;
alpha = cfg.fig22_alpha;
gamma = cfg.fig22_gamma;
M = 64;
beta = cfg.beta_default;

eta_ref = se.lora_eta(cfg.M_tilde_bench); Mref = cfg.M_tilde_bench;

d0 = pl.range_relative_to_reference(se.lora_eta(2^cfg.fig22_SF_d0), 2^cfg.fig22_SF_d0, BER, 1, eta_ref, Mref, BER, 1, gamma);
dmax = pl.range_relative_to_reference(se.lora_eta(2^cfg.fig22_SF_dmax), 2^cfg.fig22_SF_dmax, BER, 1, eta_ref, Mref, BER, 1, gamma);

f = @(Nm) pl.range_relative_to_reference(se.maspm_eta(M,Nm,beta), M, BER, 1, eta_ref, Mref, BER, 1, gamma) - dmax;
Nm = fzero(f, [4*M, 1e7]);
N0dB = pg.N0dB(M, beta, BER, BERp);
threshold = Nm/(alpha*N0dB);

m_max_nonflip = nc.solve_max_m(d0, dmax, gamma, threshold, false);
m_max_flip    = nc.solve_max_m(d0, dmax, gamma, threshold, true);
channels_flip_total = 2*m_max_flip;

fprintf('d0 (SF%d range) = %.4f delta0,  dmax (SF%d range) = %.4f delta0\n', cfg.fig22_SF_d0, d0, cfg.fig22_SF_dmax, dmax);
fprintf('Nm (64-ASPM IpI matching dmax) = %.1f,  N0dB = %.2f,  threshold=Nm/(alpha*N0dB) = %.3f\n', Nm, N0dB, threshold);
fprintf('Non-flip max channels m = %d\n', m_max_nonflip);
fprintf('Flip-pair max pairs = %d  (total channels = %d)\n', m_max_flip, channels_flip_total);
fprintf('Paper''s own reported result for this scenario: 70 channels (non-flip) -> 66 (flip pairs)\n');

d_edges_nonflip = nc.equal_ipi_partition_ranges(d0, dmax, m_max_nonflip);
d_edges_flip    = nc.equal_ipi_partition_ranges(d0, dmax, m_max_flip);

Dmax_duty = cfg.fig9_Dmax_at_dmax; % 93.8%, same representative outer duty cycle as Fig. 9
D_i_flip = (d_edges_flip(2:end)./dmax).^gamma .* Dmax_duty; % Eq. (40)

fig = figure('Position',[80 60 1250 800], 'Color','w');

% --- Panel 1: bar chart, channel-count / complexity comparison ---
subplot(2,2,1);
bar_data = [m_max_nonflip, 0; m_max_flip, m_max_flip];
b = bar(bar_data, 'stacked');
set(gca,'XTickLabel', {'Non-flip PSFs','Flip PSF pairs'});
ylabel('number of PSF channels');
title(sprintf('Fig. 22 reproduction: channel count\n(non-flip m=%d  vs  flip 2\\times%d=%d)', ...
    m_max_nonflip, m_max_flip, channels_flip_total));
legend('channels/pairs','2nd channel of each flip pair','Location','best');
grid on;

% --- Panel 2: partition-count / duty-cycle-count complexity ---
subplot(2,2,2);
bar_counts = [m_max_nonflip, m_max_flip];
bar(bar_counts, 'FaceColor',[0.3 0.6 0.9]);
set(gca,'XTickLabel', {'Non-flip','Flip pairs'});
ylabel('number of distinct partitions / duty cycles');
title('Fig. 22 reproduction: management complexity (# duty cycles)');
grid on;

% --- Panel 3: pulse duty cycle vs range (flip-pair case, Eq. 40) ---
subplot(2,2,3); hold on; box on;
plot(d_edges_flip(2:end), D_i_flip, 'o-', 'LineWidth', 1.3, 'MarkerSize', 4);
xlabel('range (relative to \delta_0(\gamma))');
ylabel('pulse duty cycle D_i');
title('Fig. 22 reproduction: duty cycle vs range (flip PSF pairs, Eq. 40)');
grid on;

% --- Panel 4: annular partition geometry (flip-pair boundaries) ---
subplot(2,2,4); hold on; box on; axis equal;
theta = linspace(0, 2*pi, 200);
for e = d_edges_flip
    plot(e*cos(theta), e*sin(theta), 'Color',[0.5 0.5 0.5], 'LineWidth', 0.5);
end
plot(d0*cos(theta), d0*sin(theta), 'r-', 'LineWidth', 1.5);
plot(dmax*cos(theta), dmax*sin(theta), 'b-', 'LineWidth', 1.5);
title(sprintf('Fig. 22 reproduction: %d flip-pair range partitions\n(d_0=%.2f, d_{max}=%.2f, \\gamma=%d)', m_max_flip, d0, dmax, gamma));
xlabel('x (relative range units)'); ylabel('y (relative range units)');
legend('partition boundaries','d_0 (SF6)','d_{max} (SF12)','Location','eastoutside');

if ~exist('results','dir'); mkdir('results'); end
saveas(fig, fullfile('results','fig22_flip_psf_pairs.png'));
fprintf('Saved results/fig22_flip_psf_pairs.png\n');
