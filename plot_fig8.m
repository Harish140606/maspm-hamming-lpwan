%PLOT_FIG8 Reproduce Fig. 8: SIR margins vs spectral efficiency and range.
%
% PAPER REFERENCE
%   Section II-C, Fig. 8 ("SIR margins for inter-PSF interference in
%   M-ASPM scale as Delta_SIR ~ eta^-1 ~ Np, and for power-law path
%   loss as Delta_SIR ~ d^gamma").
%
% EQUATIONS USED
%   Eq. (20): Delta_SIR(eta) = const/eta.
%   Eq. (21): Delta_SIR ~ eta^-1 ~ Np ~ d^gamma.
%   Eq. (26): range <-> spectral efficiency mapping (via path_loss_model.m).
%
% PARAMETERS USED
%   M = 16, 64; BER'/BER = 1e-3/1e-4 (cfg.BER, cfg.BER_prime); gamma=2.
%
% SIMULATION EXPLANATION
%   Top panel: Delta_SIR (dB) vs eta^-1 = B/fb (dB), a direct plot of
%   Eq. (20)-(21) (linear relationship in dB-dB, i.e. a straight line
%   of unit slope, since Delta_SIR ~ eta^-1).
%   Bottom panel: the same Delta_SIR values, but with the x-axis
%   converted from eta to "range relative to LoRa SF=7" via Eq. (26)
%   (path_loss_model.range_relative_to_reference), with LoRa SF6-12
%   tick reference points annotated exactly as in the paper's Fig. 8.

clear; clc;
cfg = config();
se = spectral_efficiency();
pg = processing_gain();
pl = path_loss_model();

BER = cfg.BER; BERp = cfg.BER_prime; gamma = cfg.gamma; beta = cfg.beta_default;
Ms = [16, 64];
colors = lines(numel(Ms));
eta_ref = se.lora_eta(cfg.M_tilde_bench); Mref = cfg.M_tilde_bench;

fig = figure('Position',[100 100 800 850], 'Color','w');

% ---- Top panel: Delta_SIR vs eta^-1 (B/fb) ----
subplot(2,1,1); hold on; box on;
h1 = gobjects(1,numel(Ms));
for mi = 1:numel(Ms)
    M = Ms(mi);
    Np_min = cfg.Np_min_factor*M;
    Np_sweep = Np_min .* logspace(0, 2, 200);
    eta_sweep = se.maspm_eta(M, Np_sweep, beta);
    delta_sir = pg.delta_sir_of_eta(eta_sweep, M, BER, BERp);
    h1(mi) = plot(10*log10(1./eta_sweep), 10*log10(delta_sir), '-', 'Color', colors(mi,:), 'LineWidth',1.6);
end
xlabel('\eta^{-1} = B/f_b  (dB)');
ylabel('\Delta_{SIR}  (dB)');
title(sprintf('Fig. 8 (top) reproduction: \\Delta_{SIR} vs spectral efficiency, BER''/BER=%.0e/%.0e', BERp, BER));
grid on;
legend(h1, arrayfun(@(m) sprintf('%d-ASPM',m), Ms, 'UniformOutput', false), 'Location','northwest');

% ---- Bottom panel: Delta_SIR vs range relative to LoRa SF=7 ----
subplot(2,1,2); hold on; box on;
h2 = gobjects(1,numel(Ms));
for mi = 1:numel(Ms)
    M = Ms(mi);
    Np_min = cfg.Np_min_factor*M;
    Np_sweep = Np_min .* logspace(0, 2, 200);
    eta_sweep = se.maspm_eta(M, Np_sweep, beta);
    delta_sir = pg.delta_sir_of_eta(eta_sweep, M, BER, BERp);
    r_sweep = pl.range_relative_to_reference(eta_sweep, M, BER, 1, eta_ref, Mref, BER, 1, gamma);
    h2(mi) = plot(20*log10(r_sweep), 10*log10(delta_sir), '-', 'Color', colors(mi,:), 'LineWidth',1.6);
end
% LoRa SF tick annotations, same reference as Fig. 1/5/8 of the paper
for i = 1:numel(cfg.SF_list)
    Mt = 2^cfg.SF_list(i);
    eta_t = se.lora_eta(Mt);
    r_t = pl.range_relative_to_reference(eta_t, Mt, BER, 1, eta_ref, Mref, BER, 1, gamma);
    xline(20*log10(r_t), ':', sprintf('SF%d', cfg.SF_list(i)), 'LabelOrientation','horizontal', 'FontSize', 7);
end
xlabel('range relative to LoRa w/ SF=7  (dB)');
ylabel('\Delta_{SIR}  (dB)');
title('Fig. 8 (bottom) reproduction: \Delta_{SIR} vs range (power-law path loss)');
grid on;
legend(h2, arrayfun(@(m) sprintf('%d-ASPM',m), Ms, 'UniformOutput', false), 'Location','northwest');

if ~exist('results','dir'); mkdir('results'); end
saveas(fig, fullfile('results','fig8_sir_margins.png'));
fprintf('Saved results/fig8_sir_margins.png\n');
