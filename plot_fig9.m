%PLOT_FIG9 Reproduce Fig. 9: inter-PSF interference for 5 and 9 16-ASPM nodes.
%
% PAPER REFERENCE
%   Section III-B, Fig. 9 ("Illustration of impact of inter-PSF
%   interference for five and nine 16-ASPM transmitters placed at
%   particular distances from receiver. Pulse duty cycles are chosen
%   to ensure equality of received powers...SIR=-6dB for five nodes
%   and SIR=-9dB for nine nodes...IpI for 9 nodes is twice that for 5 nodes").
%
% EQUATIONS USED
%   Eq. (9)-(13): per-node designed pulse train and PSF (M=16, one PSF
%             channel per node, alternating up/down chirp, Sec. III-B).
%   Eq. (25): D_i/D_j = (d_i/d_j)^gamma for equal received power (used
%             here with D1 and the stated maximum duty cycle 93.8% to
%             place each node's relative distance d_i, gamma=2).
%
% PARAMETERS USED (cfg.fig9_*)
%   M=16 (cfg.M_fig9); gamma=2 (cfg.fig9_gamma, free-space, as stated);
%   D1 = 20% (5 nodes) / 10% (9 nodes) (cfg.fig9_D1); Dmax = 93.8% at
%   the outermost node in both cases (cfg.fig9_Dmax_at_dmax); IpI
%   Np = 280 (5 nodes) / 560 (9 nodes) (cfg.fig9_Np -- IMPLEMENTATION
%   ASSUMPTION: clean round numbers preserving the paper's stated 2:1
%   IpI ratio, see config.m for the full rationale).
%
% SIMULATION EXPLANATION
%   For N in {5,9}: node i's relative distance d_i is derived from its
%   duty cycle D_i (geometrically spaced so that D_1 and D_N match the
%   paper's stated values, with D_i/D_1=(d_i/d_1)^gamma, Eq. 25); each
%   node's PSF length L_i = round(D_i*Np), alternating up-/down-chirp
%   direction node-to-node (Sec. III-B guidance to avoid excess
%   peakedness of the aggregate interference, Fig. 11). All nodes'
%   waveforms (independent random payloads, same Np, no external noise
%   so SINR=SIR as in the paper's Fig. 9) are summed via
%   multi_user_channel.m, and EACH node's own matched filter is applied
%   to the composite signal, reproducing the per-node interference
%   traces on the right of Fig. 9.

clear; clc;
cfg = config();
psf = pulse_shaping_filter();
mfr = matched_filter_receiver();
mc = multi_user_channel();

M = cfg.M_fig9;
beta = cfg.beta_fig9;
Ns = cfg.Ns_fig9;
gamma = cfg.fig9_gamma;
n_spacing = 4;
deltaN = 0;
n_symbols = 4;

fig = figure('Position',[60 40 1300 950], 'Color','w');

for case_idx = 1:2
    N = cfg.fig9_n_nodes(case_idx);
    Np = cfg.fig9_Np(case_idx);
    D1 = cfg.fig9_D1(case_idx);
    Dmax = cfg.fig9_Dmax_at_dmax;

    % Node relative distances via D_i/D_1=(d_i/d_1)^gamma (Eq. 25),
    % geometrically spaced from d_1=1 (normalized) to d_N solving D_N=Dmax.
    dN_over_d1 = (Dmax/D1)^(1/gamma);
    d_rel = logspace(0, log10(dN_over_d1), N);
    D_i = D1 .* (d_rel).^gamma;
    L_i = max(8, round(D_i .* Np));

    signals = cell(1,N);
    zeta_list = cell(1,N);
    for i = 1:N
        dir_i = (-1)^(i+1); % alternate up/down chirp, Sec. III-B
        bits_i = generate_bits(n_symbols*log2(M), cfg.rng_seed + case_idx*100 + i);
        [z_i, info_i] = transmitter().run_chain(bits_i, M, Np, L_i(i), beta, Ns, dir_i, n_spacing, deltaN);
        signals{i} = z_i;
        zeta_list{i} = info_i.zeta_hat;
    end

    amplitudes = ones(1,N);   % equal PEAK amplitude; AVERAGE power set via L_i/Np (duty cycle)
    delays = zeros(1,N);
    r_composite = mc.combine(signals, amplitudes, delays);

    % Per-node matched-filter magnitude-squared output against the
    % composite (interference-affected) received signal.
    subplot(2,2,(case_idx-1)*2+1); hold on; box on;
    for i = 1:N
        offset = (i-1)*0.15;
        plot(real(signals{i}) + 3*offset, 'LineWidth', 0.8);
    end
    title(sprintf('Fig. 9(%s) reproduction: %d Tx signals (stacked), IpI=%d', ...
        char('I'+case_idx-1), N, Np));
    xlabel('sample k'); ylabel('Re\{signal\} (stacked, offset for clarity)');
    grid on;

    subplot(2,2,(case_idx-1)*2+2); hold on; box on;
    peak_mags = zeros(1,N);
    for i = 1:N
        y2_i = mfr.detect(r_composite, zeta_list{i});
        offset = (i-1)*1.3;
        plot(y2_i./max(y2_i(:)) + offset, 'LineWidth', 0.9);
        peak_mags(i) = max(y2_i);
    end
    title(sprintf('Fig. 9(%s) reproduction: per-node |(I+iQ)*(g_i-ih_i)|^2 (stacked)', char('I'+case_idx-1)));
    xlabel('sample k'); ylabel('normalized y_{nc}^2 (stacked, offset)');
    grid on;

    fprintf('N=%d nodes: d_rel=[%s]\n', N, sprintf('%.3f ', d_rel));
    fprintf('N=%d nodes: D_i=[%s] (target D1=%.2f, Dmax=%.3f)\n', N, sprintf('%.3f ', D_i), D1, Dmax);
    fprintf('N=%d nodes: L_i=[%s]\n', N, sprintf('%d ', L_i));
end

if ~exist('results','dir'); mkdir('results'); end
saveas(fig, fullfile('results','fig9_interference_illustration.png'));
fprintf('Saved results/fig9_interference_illustration.png\n');
