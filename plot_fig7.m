%PLOT_FIG7 Reproduce Fig. 7: Product eta*Delta_SIR vs BER.
%
% PAPER REFERENCE
%   Section II-C, Fig. 7 ("Product eta*Delta_SIR in relation to chosen
%   constraints on BER and BER'").
%
% EQUATIONS USED
%   Eq. (19): Pb(Gamma/eta;M)=BER=Pb(Gamma'/eta;M)*(BER/BER').
%   Eq. (20): Delta_SIR = eta^-1 * const, const = 1/rho2 - 1/rho1,
%             where rho1=required_EbN0(M,BER), rho2=required_EbN0(M,BER').
%             Hence eta*Delta_SIR = const (a function of BER, BER', M only).
%
% PARAMETERS USED
%   M = 16, 64; two families of curves as in the paper:
%     (a) solid lines: BER' fixed (1e-3 for one curve, 1e-4 for the
%         paper's second solid-line pair is actually the SAME BER' as
%         the dashed pair's ratio -- see note below), BER swept on the
%         x-axis;
%     (b) dashed lines: fixed ratio BER'/BER = 10, BER swept on the
%         x-axis (so BER' = 10*BER moves together with BER).
%   BER sweep: 1e-6 to 1e-3 (matching the paper's Fig. 7 x-axis).
%
% SIMULATION EXPLANATION
%   For each value of BER on the x-axis (log-spaced sweep), we compute
%   the SIR-margin proportionality constant const(BER) =
%   processing_gain().sir_margin_const(M,BER,BER_prime) for BER_prime
%   chosen either as a FIXED value (1e-3) or as 10*BER (fixed ratio),
%   and plot eta*Delta_SIR = const in dB.

clear; clc;
cfg = config();
pg = processing_gain();

BER_sweep = logspace(-6, -3, 200);
Ms = [16, 64];
colors = lines(numel(Ms));

fig = figure('Position',[100 100 800 600], 'Color','w'); hold on; box on;

h = gobjects(1, 2*numel(Ms));
labels = cell(1, 2*numel(Ms));
idx = 0;
for mi = 1:numel(Ms)
    M = Ms(mi);

    % Solid: fixed BER' = 1e-3
    const_fixed = arrayfun(@(b) pg.sir_margin_const(M, b, 1e-3), BER_sweep);
    idx = idx+1;
    h(idx) = plot(BER_sweep, 10*log10(const_fixed), '-', 'Color', colors(mi,:), 'LineWidth', 1.6);
    labels{idx} = sprintf('%d-ASPM  BER''=10^{-3}', M);

    % Dashed: fixed ratio BER'/BER = 10
    const_ratio = arrayfun(@(b) pg.sir_margin_const(M, b, 10*b), BER_sweep);
    idx = idx+1;
    h(idx) = plot(BER_sweep, 10*log10(const_ratio), '--', 'Color', colors(mi,:), 'LineWidth', 1.6);
    labels{idx} = sprintf('%d-ASPM  BER''/BER=10', M);
end

set(gca,'XScale','log');
xlabel('BER');
ylabel('\eta\Delta_{SIR}  (dB)');
title('Fig. 7 reproduction: Product \eta\Delta_{SIR} vs BER');
grid on;
legend(h, labels, 'Location','best');
xlim([1e-6 1e-3]);

if ~exist('results','dir'); mkdir('results'); end
saveas(fig, fullfile('results','fig7_eta_deltasir_vs_ber.png'));
fprintf('Saved results/fig7_eta_deltasir_vs_ber.png\n');
