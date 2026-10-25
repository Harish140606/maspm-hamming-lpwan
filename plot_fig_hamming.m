%PLOT_FIG_HAMMING  NOVELTY ADDITION (not in the paper): BER improvement of an
%   outer Hamming FEC layer on noncoherent M-ASPM.
%
% MOTIVATION
%   Nikitin & Davidchack (IEEE Access 2023) analyse only UNCODED M-ASPM
%   (Eq. 17, Figs. 3-4). LPWAN links almost always run a light channel code;
%   this figure quantifies what the simplest possible perfect code -- a
%   binary Hamming(n,k) block code with hard-decision syndrome decoding
%   (hamming_code.m) -- buys on top of the paper's own BER equation.
%
% WHAT IS PLOTTED
%   BER vs Eb/N0 PER INFORMATION BIT for 16-ASPM and 64-ASPM: uncoded
%   (Eq. 17, ber_calculation.m) vs + Hamming(7,4) / (15,11) / (31,26) with
%   hard-decision syndrome decoding (coded_ber_hamming.m). The code rate
%   R = k/n is charged against the code (x-axis is information-bit Eb/N0),
%   so the curves show the NET effect. Markers are independent Monte-Carlo
%   points: the real codec (hamming_code.m) + interleaver
%   (block_interleaver.m) run over the Eq. (17) channel model
%   (mary_noncoherent_channel.m).
%
% KEY FINDING (see README section 9)
%   Uncoded noncoherent M-ASPM already sits very close to the noncoherent
%   orthogonal-signalling bound, so its BER waterfall is extremely steep.
%   A hard-decision Hamming code therefore trades rate for a large
%   BER-FLOOR reduction (1-2 orders of magnitude at a fixed link budget)
%   rather than a large SNR-domain coding gain: at a fixed target BER the
%   net gain is only ~ -0.5 dB for Hamming(7,4) (rate loss dominates) rising
%   to ~ +0.3..0.7 dB for the higher-rate Hamming(31,26). The
%   fixed-SNR BER reduction is what benefits LPWAN packet-delivery /
%   fade-margin, which is the practical contribution here.

clear; clc;
cfg = config();

EbN0_dB  = 0:0.5:14;
EbN0_lin = 10.^(EbN0_dB/10);

Mset = cfg.M_list;
mset = cfg.hamming.m_list;
cols = lines(numel(Mset));
sty  = {'-','--',':','-.'};   % uncoded, then one per entry of mset

fig = figure('Position',[100 100 950 680], 'Color','w'); hold on; box on;

leg_h = gobjects(0);
leg_t = {};
for mi = 1:numel(Mset)
    M = Mset(mi);

    Pb_unc = ber_calculation(EbN0_lin, M);
    h = semilogy(EbN0_dB, Pb_unc, sty{1}, 'Color', cols(mi,:), 'LineWidth', 1.9);
    leg_h(end+1) = h; %#ok<SAGROW>
    leg_t{end+1} = sprintf('%d-ASPM  uncoded (Eq. 17)', M); %#ok<SAGROW>

    for si = 1:numel(mset)
        mm = mset(si);
        n = 2^mm - 1; k = n - mm;
        Pb_cod = coded_ber_hamming(EbN0_lin, M, mm);
        ls = sty{mod(si, numel(sty)) + 1};
        h = semilogy(EbN0_dB, Pb_cod, ls, 'Color', cols(mi,:), 'LineWidth', 1.5);
        leg_h(end+1) = h; %#ok<SAGROW>
        leg_t{end+1} = sprintf('%d-ASPM  + Hamming(%d,%d)', M, n, k); %#ok<SAGROW>
    end
end

% --- Monte-Carlo validation (real codec over the Eq. 17 channel) ---
rng(cfg.rng_seed);
mc_dB = 2:2:12;
hc = hamming_code(4);   % validate the Hamming(15,11) curve (the useful net-gain case)
for mi = 1:numel(Mset)
    M   = Mset(mi);
    bps = round(log2(M));
    ber_u = nan(size(mc_dB));
    ber_c = nan(size(mc_dB));
    for j = 1:numel(mc_dB)
        EbN0 = 10^(mc_dB(j)/10);

        nb = ceil(4e5 / bps) * bps;
        tb = randi([0 1], 1, nb);
        ber_u(j) = mean(tb ~= mary_noncoherent_channel(tb, EbN0, M));

        B   = bps * 5000;                     % Hamming blocks -> coded length is a multiple of bps
        ib  = randi([0 1], 1, hc.k * B);
        cb  = hc.encode(ib);
        itl = block_interleaver(B, hc.n);     % spread symbol-error bursts across codewords
        rcb = itl.inv(mary_noncoherent_channel(itl.fwd(cb), EbN0 * hc.rate, M));
        ber_c(j) = mean(hc.decode(rcb) ~= ib);
    end
    semilogy(mc_dB, ber_u, 'o', 'Color', cols(mi,:), ...
        'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    semilogy(mc_dB, ber_c, 's', 'Color', cols(mi,:), ...
        'MarkerFaceColor', cols(mi,:), 'MarkerSize', 6, 'HandleVisibility', 'off');
end

set(gca, 'YScale', 'log');
grid on;
xlim([min(EbN0_dB) max(EbN0_dB)]);
ylim([1e-6 1]);
xlabel('E_b/N_0 per information bit (dB)');
ylabel('BER');
title({'NOVELTY: outer Hamming FEC on noncoherent M-ASPM (rate loss charged)', ...
       'lines = analytic (Eq. 17 + bounded-distance decoding);  markers = Monte-Carlo (circles uncoded, squares Hamming(15,11))'}, ...
      'Interpreter', 'none');
legend(leg_h, leg_t, 'Location', 'southwest', 'FontSize', 8);

% --- coding-gain summary at the paper's operating BERs ---
fprintf('--- Novelty: Hamming coding gain (net of rate loss) ---\n');
for target = [cfg.BER_prime, cfg.BER]
    for mi = 1:numel(Mset)
        M = Mset(mi);
        e_unc = ebn0_at_ber(EbN0_dB, ber_calculation(EbN0_lin, M), target);
        for si = 1:numel(mset)
            mm = mset(si); n = 2^mm - 1; k = n - mm;
            e_cod = ebn0_at_ber(EbN0_dB, coded_ber_hamming(EbN0_lin, M, mm), target);
            fprintf('  BER=%.0e  %d-ASPM  Hamming(%2d,%2d): gain = %+5.2f dB  (uncoded %.2f dB -> coded %.2f dB)\n', ...
                target, M, n, k, e_unc - e_cod, e_unc, e_cod);
        end
    end
end

if ~exist('results', 'dir'); mkdir('results'); end
saveas(fig, fullfile('results', 'fig_hamming_coded_ber.png'));
fprintf('Saved results/fig_hamming_coded_ber.png\n');
