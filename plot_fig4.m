%PLOT_FIG4 Reproduce Fig. 4: Uncoded BER vs SNR for LoRa and M-ASPM.
%
% PAPER REFERENCE
%   Section II-B, Fig. 4 ("Uncoded BER vs SNR performances of LoRa
%   (dashed lines), single-sideband 16-ASPM (solid lines), and
%   single-sideband 64-ASPM (dotted lines) for noncoherent detection
%   in AWGN channel").
%
% EQUATIONS USED
%   Eq. (17): Pb(Gamma/eta; M) -- the noncoherent BER formula, via
%             ber_calculation.m.
%   Eq. (14): eta_tilde(M_tilde) = log2(M_tilde)/M_tilde for LoRa.
%   Eq. (18): equal-range/equal-BER matching between LoRa SF and
%             M-ASPM spectral efficiency, used to choose the family of
%             M-ASPM curves shown (one per LoRa SF, at the spectral
%             efficiency that makes M-ASPM hit the SAME target BER at
%             the SAME range/SNR as that LoRa SF -- exactly the
%             benchmark construction described in Sec. III-A).
%
% PARAMETERS USED
%   BER target = 1e-4 (cfg.BER), SF = 6..12 (cfg.SF_list), M = 16, 64
%   (cfg.M_list), SNR sweep cfg.SNR_dB_range.
%
% SIMULATION EXPLANATION
%   For each LoRa spreading factor SF, we compute the SNR Gamma_SF at
%   which LoRa achieves the target BER (via Eq. 17/18). We then find
%   the M-ASPM spectral efficiency eta_M (for M=16 and M=64) that
%   achieves the SAME target BER at that SAME Gamma_SF (again Eq. 18):
%   since Pb depends only on Gamma/eta, this is simply
%   eta_M = Gamma_SF / required_EbN0(M,BER). Plotting Pb(Gamma;M,eta_M)
%   over the full SNR sweep for each of these eta_M values reproduces
%   the family of M-ASPM curves shown in Fig. 4, each anchored to a
%   different LoRa SF's operating point, exactly as the paper's own
%   B/fb legend values are constructed (Sec. III-A, Eq. 26 discussion).
%
% NOTE ON EXACT LEGEND VALUES
%   The paper's legend lists specific B/fb (=1/eta) values (e.g.
%   14.1, 27.1, ... for 16-ASPM). These are reproduced here to within
%   the precision of Eq. (17)-(18) root-finding; minor numeric
%   differences from the paper's published values are expected because
%   the exact rounding/selection procedure used to pick those legend
%   numbers is not fully specified in the paper text.

clear; clc;
cfg = config();
se = spectral_efficiency();

SNR_dB = cfg.SNR_dB_range;
SNR_lin = 10.^(SNR_dB/10);

fig = figure('Position',[100 100 900 650], 'Color','w'); hold on; box on;
colors_lora = lines(numel(cfg.SF_list));

% --- LoRa curves (dashed) ---
h_lora = gobjects(1,numel(cfg.SF_list));
for i = 1:numel(cfg.SF_list)
    SF = cfg.SF_list(i);
    Mt = 2^SF;
    eta_t = se.lora_eta(Mt);
    Pb = ber_calculation(SNR_lin./eta_t, Mt);
    h_lora(i) = semilogy(SNR_dB, Pb, '--', 'Color', colors_lora(i,:), 'LineWidth', 1.3);
end

% --- M-ASPM curves, anchored to each LoRa SF operating point ---
line_styles = {'-','-.'};
markers = {'','',};
h_maspm = {};
labels_maspm = {};
for mi = 1:numel(cfg.M_list)
    M = cfg.M_list(mi);
    rho_M = se.required_EbN0(M, cfg.BER);
    for i = 1:numel(cfg.SF_list)
        SF = cfg.SF_list(i);
        Mt = 2^SF;
        eta_t = se.lora_eta(Mt);
        rho_t = se.required_EbN0(Mt, cfg.BER);
        Gamma_SF = eta_t * rho_t;             % SNR at which LoRa SF hits target BER
        eta_M = Gamma_SF / rho_M;              % M-ASPM eta matched to same (Gamma,BER)
        Pb = ber_calculation(SNR_lin./eta_M, M);
        hh = semilogy(SNR_dB, Pb, line_styles{mi}, 'Color', colors_lora(i,:), 'LineWidth', 1.1);
        h_maspm{end+1} = hh; %#ok<SAGROW>
        labels_maspm{end+1} = sprintf('%d-ASPM  B/f_b\\approx%.3g  (SF%d anchor)', M, 1/eta_M, SF); %#ok<SAGROW>
    end
end

set(gca, 'YScale', 'log');
ylim([1e-5 1]);
xlim([min(SNR_dB) max(SNR_dB)]);
xlabel('SNR = \Gamma (dB)');
ylabel('BER');
title('Fig. 4 reproduction: Uncoded BER vs SNR -- LoRa (dashed) vs 16-/64-ASPM (solid/dash-dot)');
grid on;

leg_labels = [arrayfun(@(s) sprintf('LoRa SF%d', s), cfg.SF_list, 'UniformOutput', false), labels_maspm];
leg_handles = [h_lora, h_maspm{:}];
legend(leg_handles, leg_labels, 'Location', 'eastoutside', 'FontSize', 7);

fprintf('--- Fig. 4 data summary ---\n');
for mi = 1:numel(cfg.M_list)
    M = cfg.M_list(mi);
    rho_M = se.required_EbN0(M, cfg.BER);
    fprintf('%d-ASPM required EbN0 for BER=%.0e: %.3f (%.2f dB)\n', M, cfg.BER, rho_M, 10*log10(rho_M));
end

if ~exist('results', 'dir'); mkdir('results'); end
saveas(fig, fullfile('results','fig4_ber_vs_snr.png'));
fprintf('Saved results/fig4_ber_vs_snr.png\n');
