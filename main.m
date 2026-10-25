%MAIN Top-level driver for the M-ASPM MATLAB reproduction project.
%
% PURPOSE
%   Runs the complete project end to end:
%     1. Loads configuration (config.m).
%     2. Runs a bit-true transmit/receive demonstration (encoder ->
%        PSF pulse shaping -> AWGN channel -> matched-filter receiver
%        -> decoder) and reports the resulting BER against the
%        closed-form theoretical value (Eq. 17), for both M=16 and
%        M=64.
%     3. Prints a summary of the key analytical quantities reproduced
%        from the paper (required Eb/N0, spectral efficiencies, SIR
%        margins, N0dB) for M=16 and M=64.
%     4. Generates every figure reproduction (Figs. 2, 4, 5, 6, 7, 8,
%        9, 22), saving each as a PNG under results/.
%     5. NOVELTY ADDITION (not in the paper): adds an outer binary
%        Hamming(n,k) FEC layer around the M-ASPM bit stream and
%        quantifies its coding gain -- bit-stream Monte-Carlo over the
%        Eq. (17) channel (Part 1b), closed-form coding gain (Part 2b),
%        and a dedicated figure (plot_fig_hamming.m). See README
%        section 9. Toggle with cfg.hamming.enabled.
%
% USAGE
%   From the project root, in MATLAB:  >> main
%
% REQUIRED TOOLBOXES
%   - Signal Processing Toolbox   (besseli scaling, hilbert, xcorr)
%   - Communications Toolbox      (rcosdesign)
%   - Symbolic Math Toolbox       (only used defensively in earlier
%     development of ber_calculation.m; the shipped version uses a
%     pure double-precision Bessel-integral evaluation and does NOT
%     require the Symbolic Math Toolbox at run time)
%
% See README.md for the full project description, equation-by-equation
% documentation, implementation-assumption list, and numerical
% validation results.

clear; close all; clc;

fprintf('=========================================================\n');
fprintf(' M-ASPM MATLAB Reproduction Project\n');
fprintf(' Nikitin & Davidchack, IEEE Access, vol. 11, pp. 25642-25661, 2023\n');
fprintf('=========================================================\n\n');

cfg = config();

if ~exist('results', 'dir'); mkdir('results'); end

% ---------------------------------------------------------------
% Part 1: bit-true transmit/receive demonstration and BER check
% ---------------------------------------------------------------
fprintf('--- Part 1: bit-true chain demonstration (Eq. 9-13, Fig. 2-III) ---\n');
tx = transmitter();
mfr = matched_filter_receiver();
ch = awgn_channel();
pd = pulse_detector();

for M = cfg.M_list
    Np = 20*M;             % comfortably above Np_min=4M
    L = round(0.6*Np);      % substantial duty cycle for a clean demo
    beta = cfg.beta_default;
    Ns = 2/(1+beta);
    n_spacing = cfg.pulse_position_spacing;

    n_bits = 500*log2(M);
    bits = generate_bits(n_bits, cfg.rng_seed);
    [z, info] = tx.run_chain(bits, M, Np, L, beta, Ns, 1, n_spacing, 0);

    EbN0dB_demo = 8;
    EbN0_demo = 10^(EbN0dB_demo/10);
    rng(cfg.rng_seed);
    r = ch.add_noise(z, EbN0_demo, log2(M), Np, L);
    y2 = mfr.detect(r, info.zeta_hat);
    delta_k_hat = pd.find_peaks(y2, numel(info.mj), Np, 0, L-1);
    [bits_hat, ~] = aspm_decoder(delta_k_hat, M, n_spacing);

    ber_sim = sum(bits_hat ~= bits)/n_bits;
    ber_theory = ber_calculation(EbN0_demo, M);
    fprintf('M=%3d: constant-envelope=%d, EbN0=%d dB -> BER_sim=%.5f (theory Eq.17=%.5f)\n', ...
        M, info.constant_envelope_ok, EbN0dB_demo, ber_sim, ber_theory);
end
fprintf(['NOTE: BER_sim is somewhat higher than the idealized Eq. (17) theory\n' ...
         'because the bit-true PSF chirp uses an approximate NLFM design\n' ...
         '(pulse_shaping_filter.m, documented Implementation Assumption); see\n' ...
         'README.md "Numerical Validation" and "Known Limitations".\n\n']);

% ---------------------------------------------------------------
% Part 1b: NOVELTY ADDITION -- outer Hamming FEC (not in the paper)
% Runs the real interleaved Hamming codec (hamming_code.m +
% block_interleaver.m) end to end over the paper's own noncoherent
% M-ary channel model (Eq. 17, Monte-Carlo via
% mary_noncoherent_channel.m), so its error-correction behaviour is
% measured on genuinely independent channel errors rather than
% through the approximate PSF chirp. This part compares uncoded vs
% coded AT EQUAL CHANNEL Eb/N0 (i.e. equal energy per transmitted
% bit) to isolate the decoder's correction capability; the NET effect
% once the rate loss R = k/n is charged is reported analytically in
% Part 2b and plotted by plot_fig_hamming.m. See README.md section 9.
% ---------------------------------------------------------------
if isfield(cfg, 'hamming') && cfg.hamming.enabled
    fprintf('--- Part 1b: NOVELTY -- interleaved Hamming FEC (bit-stream Monte-Carlo, equal channel Eb/N0) ---\n');
    hc = hamming_code(cfg.hamming.m);
    fprintf('Outer code: Hamming(%d,%d), rate R = %.4f\n', hc.n, hc.k, hc.rate);
    for M = cfg.M_list
        bps = log2(M);
        for EbN0dB = [4 6]
            EbN0 = 10^(EbN0dB/10);

            nb = ceil(3e5/bps)*bps;
            tb = generate_bits(nb, cfg.rng_seed + EbN0dB);
            ber_u = mean(tb ~= mary_noncoherent_channel(tb, EbN0, M, cfg.rng_seed + EbN0dB + 1));

            B   = bps*4000;                      % Hamming blocks -> coded length is a multiple of bps
            ib  = generate_bits(hc.k*B, cfg.rng_seed + EbN0dB + 2);
            cb  = hc.encode(ib);
            itl = block_interleaver(B, hc.n);     % spread symbol-error bursts across codewords
            rcb = itl.inv(mary_noncoherent_channel(itl.fwd(cb), EbN0, M, cfg.rng_seed + EbN0dB + 3));
            ber_c = mean(hc.decode(rcb) ~= ib);

            if ber_c > 0
                ratio_str = sprintf('%.0fx fewer errors', ber_u/ber_c);
            else
                ratio_str = 'no coded errors observed';
            end
            fprintf('  M=%3d  channel Eb/N0=%2d dB:  BER uncoded=%.3e   Hamming(%d,%d)=%.3e   (%s)\n', ...
                    M, EbN0dB, ber_u, hc.n, hc.k, ber_c, ratio_str);
        end
    end
    fprintf('\n');
end

% ---------------------------------------------------------------
% Part 2: key analytical quantities (Eqs. 14-21, 26)
% ---------------------------------------------------------------
fprintf('--- Part 2: analytical summary (Eqs. 14-21, 26) ---\n');
se = spectral_efficiency();
pg = processing_gain();
for M = cfg.M_list
    rho = se.required_EbN0(M, cfg.BER);
    eta_max = se.maspm_eta_max(M, cfg.beta_default);
    N0dB = pg.N0dB(M, cfg.beta_default, cfg.BER, cfg.BER_prime);
    fprintf('M=%3d: required Eb/N0=%.3f (%.2f dB) for BER=%.0e | eta_max=%.5f | N0dB=%.2f\n', ...
        M, rho, 10*log10(rho), cfg.BER, eta_max, N0dB);
end
fprintf('\n');

% ---------------------------------------------------------------
% Part 2b: NOVELTY ADDITION -- closed-form Hamming coding gain
% (Eq. 17 + bounded-distance block decoding, coded_ber_hamming.m).
% ---------------------------------------------------------------
if isfield(cfg, 'hamming') && cfg.hamming.enabled
    fprintf('--- Part 2b: NOVELTY -- closed-form Hamming coding gain (Eq. 17 + FEC) ---\n');
    EbN0_grid_dB = 0:0.05:16;
    EbN0_grid    = 10.^(EbN0_grid_dB/10);
    for M = cfg.M_list
        e_unc = ebn0_at_ber(EbN0_grid_dB, ber_calculation(EbN0_grid, M), cfg.BER);
        for mm = cfg.hamming.m_list
            n = 2^mm - 1; k = n - mm;
            e_cod = ebn0_at_ber(EbN0_grid_dB, coded_ber_hamming(EbN0_grid, M, mm), cfg.BER);
            fprintf('  M=%3d: Hamming(%2d,%2d) net coding gain at BER=%.0e = %+.2f dB\n', ...
                M, n, k, cfg.BER, e_unc - e_cod);
        end
    end
    fprintf('\n');
end

% ---------------------------------------------------------------
% Part 3: generate every figure reproduction
% ---------------------------------------------------------------
fprintf('--- Part 3: generating figures (this may take a minute) ---\n');
% NOTE: each plot_fig*.m script begins with its own 'clear' statement
% (so it can also be run standalone). Because MATLAB's run() executes
% a script in the CALLER's workspace, looping over a cell array of
% script names here would have that 'clear' wipe the loop variable
% itself. We therefore call each script on its own line rather than
% through a loop, so no state needs to survive a script's internal
% 'clear'.
% NOTE: figures are intentionally left open (not closed) after each
% script so that all 8 figure windows remain visible once main.m
% finishes; each is still saved to results/ as it is created.
fprintf('Running plot_fig2.m ...\n');  run('plot_fig2.m');
fprintf('Running plot_fig4.m ...\n');  run('plot_fig4.m');
fprintf('Running plot_fig5.m ...\n');  run('plot_fig5.m');
fprintf('Running plot_fig6.m ...\n');  run('plot_fig6.m');
fprintf('Running plot_fig7.m ...\n');  run('plot_fig7.m');
fprintf('Running plot_fig8.m ...\n');  run('plot_fig8.m');
fprintf('Running plot_fig9.m ...\n');  run('plot_fig9.m');
fprintf('Running plot_fig22.m ...\n'); run('plot_fig22.m');
% NOTE: the plot_fig*.m scripts above each begin with their own 'clear',
% so 'cfg' from Part 1/2 no longer exists here -- reload it before the
% (novelty) Hamming figure's enable check.
cfg = config();
if isfield(cfg, 'hamming') && cfg.hamming.enabled
    fprintf('Running plot_fig_hamming.m (NOVELTY) ...\n'); run('plot_fig_hamming.m');
end

fprintf('\nAll figures generated under results/. Done.\n');
