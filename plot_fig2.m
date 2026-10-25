%PLOT_FIG2 Reproduce Fig. 2: single-sideband M-ary ASPM link illustration.
%
% PAPER REFERENCE
%   Section II, Fig. 2 ("Illustration of single-sideband M-ary ASPM
%   link with constant-envelope pulses and noncoherent detection").
%
% EQUATIONS USED
%   Eq. (9): designed pulse train x_hat[k].
%   Eq. (10): pulse arrival times kj (with the M=8 example of the
%             paper's own Fig. 2, 3-bit symbols shown as binary labels).
%   Eq. (11): constant-envelope PSF chirp zeta_hat[k] (real/imag parts g,h).
%   Eq. (12): pulse-shaped transmit signal z[k]=x_g[k]+i*x_h[k].
%   Eq. (13): constant-envelope check on z[k].
%   Fig. 2(III) formula: noncoherent detection statistic y_nc^2[k],
%             implemented via matched_filter_receiver.m.
%
% PARAMETERS USED
%   M = cfg.M_fig2 (=8, matching the paper's own Fig. 2 example),
%   Np, L chosen for clear visual separation of a short demo payload.
%
% SIMULATION EXPLANATION
%   A short random payload (8 symbols) is run through the full
%   transmitter -> AWGN channel -> matched-filter receiver chain, and
%   each stage's waveform is plotted in a stacked-panel layout mirroring
%   Fig. 2's three blocks: (I) designed pulse train and PSF, (II) the
%   transmitted constant-envelope pulse-shaped waveform, (III) the
%   noisy received signal and the noncoherent detection statistic
%   y_nc^2[k], with the correctly detected pulse positions marked.

clear; clc;
cfg = config();

M = cfg.M_fig2;          % 8, as in the paper's own Fig. 2 example
Np = 48;
L = 16;
beta = cfg.beta_fig9;
Ns = 2/(1+beta);
n_spacing = 4; % Eq. (27) spacing "n"
deltaN = 0;
n_symbols = 7;

bits = generate_bits(n_symbols*log2(M), cfg.rng_seed);

tx = transmitter();
[z, info] = tx.run_chain(bits, M, Np, L, beta, Ns, 1, n_spacing, deltaN);

mfr = matched_filter_receiver();
ch = awgn_channel();
rng(cfg.rng_seed);
EbN0_demo = 10^(10/10); % 10 dB, for a clean illustrative received waveform
r = ch.add_noise(z, EbN0_demo, log2(M), Np, L);
y2 = mfr.detect(r, info.zeta_hat);

pd = pulse_detector();
filter_delay = L-1;
delta_k_hat = pd.find_peaks(y2, n_symbols, Np, deltaN, filter_delay);
[bits_hat, mj_hat] = aspm_decoder(delta_k_hat, M, n_spacing);
fprintf('Fig.2 demo: symbol errors = %d/%d, bit errors = %d/%d\n', ...
    sum(mj_hat ~= info.mj), n_symbols, sum(bits_hat ~= bits), numel(bits));

fig = figure('Position',[80 60 1100 900], 'Color','w');

% --- Block I: designed pulse train x_hat[k] and PSF ---
subplot(4,1,1); hold on; box on;
kplot = 0:(numel(info.x_hat)-1);
stem(kplot, info.x_hat, 'filled', 'MarkerSize', 3);
for j = 1:n_symbols
    text(info.kj(j), info.signs(j)*1.25, sprintf('%s', dec2bin(info.mj(j)-1, log2(M))), ...
        'FontSize', 8, 'HorizontalAlignment','center');
end
ylim([-1.6 1.6]);
xlabel('sample k'); ylabel('x\_hat[k]');
title(sprintf('Fig. 2(I) reproduction: designed pulse train (Eq. 9-10), M=%d, symbol labels shown', M));
grid on;

% --- Block I (cont.): PSF real/imag parts ---
subplot(4,1,2); hold on; box on;
plot(0:L-1, info.g_hat, 'b-', 'LineWidth',1.3);
plot(0:L-1, info.h_hat, 'r-', 'LineWidth',1.3);
xlabel('sample k'); ylabel('amplitude');
title('Fig. 2(I) reproduction: PSF chirp \zeta\_hat[k] = g\_hat[k] + i h\_hat[k]  (Eq. 11)');
legend('g\_hat[k] (real)', 'h\_hat[k] = Hilbert\{g\_hat\}[k] (imag)', 'Location','best');
grid on;

% --- Block II: transmitted constant-envelope waveform ---
subplot(4,1,3); hold on; box on;
kz = 0:(numel(z)-1);
plot(kz, real(z), 'b-'); plot(kz, imag(z), 'r-'); plot(kz, abs(z), 'k-','LineWidth',1.3);
xlabel('sample k'); ylabel('amplitude');
title(sprintf('Fig. 2(II) reproduction: transmitted signal z[k]=x_g[k]+i x_h[k] (Eq. 12-13); constant envelope: %d', info.constant_envelope_ok));
legend('x_g[k] (real)','x_h[k] (imag)','|z[k]| (envelope)','Location','best');
grid on;

% --- Block III: received signal + noncoherent detection statistic ---
subplot(4,1,4); hold on; box on;
kr = 0:(numel(r)-1);
plot(kr, real(r), 'Color',[0.6 0.6 0.6]);
ky2 = 0:(numel(y2)-1);
plot(ky2, y2./max(y2), 'm-', 'LineWidth', 1.3);
det_pos = (0:n_symbols-1).*Np + deltaN + filter_delay + delta_k_hat;
plot(det_pos, ones(size(det_pos))*1.02, 'kv', 'MarkerFaceColor','k', 'MarkerSize',5);
xlabel('sample k'); ylabel('amplitude (Re part) / normalized y_{nc}^2');
title('Fig. 2(III) reproduction: Rx signal (gray) and noncoherent detection statistic y_{nc}^2[k] (magenta), detected pulses (\nabla)');
legend('Re\{r[k]\} (Rx signal+noise)', 'y_{nc}^2[k] (normalized)', 'detected pulse positions', 'Location','best');
grid on;

if ~exist('results','dir'); mkdir('results'); end
saveas(fig, fullfile('results','fig2_waveform_illustration.png'));
fprintf('Saved results/fig2_waveform_illustration.png\n');
