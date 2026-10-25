%PLOT_FIG5 Reproduce Fig. 5: Spectral efficiency vs range, LoRa and M-ASPM.
%
% PAPER REFERENCE
%   Section II-B, Fig. 5 ("Spectral efficiency vs range for LoRa and
%   M-ASPM (noncoherent detection)").
%
% EQUATIONS USED
%   Eq. (14): LoRa spectral efficiency eta_tilde(M_tilde)=log2(M_tilde)/M_tilde.
%   Eq. (15)-(16): M-ASPM spectral efficiency and its maximum.
%   Eq. (17)-(18): BER-based range matching (via path_loss_model.m).
%   Eq. (21),(26): range ~ eta^(-1/gamma) for power-law path loss.
%
% PARAMETERS USED
%   BER = 1e-4 (cfg.BER), gamma = 2 (cfg.gamma), M = 16, 256 (a stand-in
%   for the paper's "256-ASPM" curve, added alongside 16- and 64-ASPM
%   as in Fig. 5), beta = 1/3 (cfg.beta_default).
%
% SIMULATION EXPLANATION
%   For a sweep of M-ASPM spectral efficiencies eta (equivalently, IpI
%   values Np via Eq. 15), the achieved range is computed via
%   path_loss_model's range_relative_to_reference, which implements
%   Eq. (26): d/delta0(gamma) = (D*eta0/eta)^(1/gamma) with D=1 (full
%   duty cycle, spread-spectrum IpI-controlled range, Sec. II-B). LoRa
%   points are placed on the same axis using the same range-matching
%   relation with eta = eta_tilde(2^SF), one point per SF.
%
% NOTE
%   For M-ASPM, the achievable spectral efficiency is capped above by
%   eta_max(M) = eta_tilde(M)/(1+beta) (Eq. 16, at Np=4M, the shortest
%   allowed IpI) -- the M-ASPM curves in this reproduction start at
%   that maximum efficiency (shortest range plotted) and decrease
%   toward longer range as Np increases, matching the paper's Fig. 5
%   shape (M-ASPM curves are flat/capped at short range, then slope
%   down at longer range, exactly as in the source figure).

clear; clc;
cfg = config();
se = spectral_efficiency();
pl = path_loss_model();

BER = cfg.BER;
gamma = cfg.gamma;
beta = cfg.beta_default;
eta_ref = se.lora_eta(cfg.M_tilde_bench);
Mref = cfg.M_tilde_bench;

fig = figure('Position',[100 100 850 600], 'Color','w'); hold on; box on;

% --- LoRa points (piecewise-constant per SF, as in Fig. 5) ---
range_lora = zeros(1,numel(cfg.SF_list));
eta_lora = zeros(1,numel(cfg.SF_list));
for i = 1:numel(cfg.SF_list)
    Mt = 2^cfg.SF_list(i);
    eta_lora(i) = se.lora_eta(Mt);
    range_lora(i) = pl.range_relative_to_reference(eta_lora(i), Mt, BER, 1, eta_ref, Mref, BER, 1, gamma);
end
% Draw LoRa as horizontal steps (spectral efficiency constant across
% the range interval associated with each SF), matching Fig. 5's style.
for i = 1:numel(cfg.SF_list)
    r_lo = range_lora(i);
    if i < numel(cfg.SF_list)
        r_hi = range_lora(i+1);
    else
        r_hi = range_lora(i)*3;
    end
    plot([r_lo r_hi], [eta_lora(i) eta_lora(i)], 'k--', 'LineWidth', 1.3, 'HandleVisibility','off');
    plot(r_lo, eta_lora(i), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 4, 'HandleVisibility','off');
    text(r_lo, eta_lora(i)*1.15, sprintf('SF%d', cfg.SF_list(i)), 'FontSize', 8, 'HorizontalAlignment','center');
end
h_lora = plot(nan, nan, 'k--', 'LineWidth', 1.3); % legend proxy

% --- M-ASPM curves ---
M_curve_list = [16, 64, 256];
colors = lines(numel(M_curve_list));
h_maspm = gobjects(1,numel(M_curve_list));
for mi = 1:numel(M_curve_list)
    M = M_curve_list(mi);
    eta_max = se.maspm_eta_max(M, beta);
    Np_min = cfg.Np_min_factor * M;
    Np_sweep = Np_min .* logspace(0, 2.5, 200); % IpI from minimum up to 300x minimum
    eta_sweep = se.maspm_eta(M, Np_sweep, beta);
    r_sweep = pl.range_relative_to_reference(eta_sweep, M, BER, 1, eta_ref, Mref, BER, 1, gamma);
    h_maspm(mi) = plot(r_sweep, eta_sweep, '-', 'Color', colors(mi,:), 'LineWidth', 1.6);
    fprintf('%3d-ASPM: eta_max=%.5f at Np_min=%d, range/delta0=%.3f\n', M, eta_max, Np_min, r_sweep(1));
end

set(gca,'XScale','log','YScale','log');
xlabel('range relative to LoRa w/ SF=7 (linear ratio; paper axis is in dB)');
ylabel('spectral efficiency \eta');
title('Fig. 5 reproduction: Spectral efficiency vs range (power-law path loss, \gamma=2)');
grid on;
legend([h_lora, h_maspm], [{'LoRa'}, arrayfun(@(m) sprintf('%d-ASPM', m), M_curve_list, 'UniformOutput', false)], ...
    'Location','southwest');

if ~exist('results','dir'); mkdir('results'); end
saveas(fig, fullfile('results','fig5_spectral_efficiency_vs_range.png'));
fprintf('Saved results/fig5_spectral_efficiency_vs_range.png\n');
