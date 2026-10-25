function rx_bits = mary_noncoherent_channel(tx_bits, EbN0_lin, M, seed)
%MARY_NONCOHERENT_CHANNEL Monte-Carlo of the noncoherent M-ary orthogonal
%   detection channel that underlies Eq. (17) of the paper (square-law
%   combining, pick-max decision), operating directly on a bit stream.
%
% NOVELTY-SUPPORT UTILITY.
%   ber_calculation.m gives the closed-form Eq. (17) BER; this function is
%   its sample-level Monte-Carlo counterpart. It lets the outer Hamming
%   codec (hamming_code.m) be validated end-to-end against GENUINELY
%   INDEPENDENT channel bit errors, without invoking the approximate PSF
%   chirp of the bit-true waveform chain (whose residual non-orthogonality
%   would otherwise contaminate the coding-gain measurement -- see README
%   "Known Limitations"). Used by plot_fig_hamming.m and main.m Part 1b.
%
% CHANNEL MODEL (identical statistics to Eq. 17 / ber_calculation.m):
%   Per symbol, per-symbol SNR gamma_s = log2(M) * (Eb/N0). The branch for
%   the transmitted symbol yields  |sqrt(gamma_s) + w|^2  with w ~ CN(0,1);
%   the other M-1 branches yield i.i.d. Exp(1) (central chi-square, 2 dof).
%   A symbol error occurs when any noise-only branch exceeds the signal
%   branch. On a symbol error the received symbol is drawn uniformly from
%   the M-1 wrong symbols (no Gray mapping assumed -- consistent with the
%   orthogonal-signalling Pb = M/(2(M-1))*Ps relation behind Eq. 17).
%
% INPUTS
%   tx_bits  - 1xN row vector of 0/1 (N a multiple of log2(M)). MSB-first
%              grouping, matching aspm_encoder.m's convention.
%   EbN0_lin - Eb/N0 (LINEAR) per bit of tx_bits.
%   M        - alphabet size (power of two).
%   seed     - optional RNG seed for reproducibility.
%
% OUTPUT
%   rx_bits  - 1xN received bit vector (same length as tx_bits).
%
% Requires no MATLAB toolboxes.

if nargin >= 4 && ~isempty(seed)
    rng(seed);
end

bps = round(log2(M));
if abs(bps - log2(M)) > 1e-9
    error('mary_noncoherent_channel:M', 'M must be a power of two.');
end
tx_bits = tx_bits(:).';
if mod(numel(tx_bits), bps) ~= 0
    error('mary_noncoherent_channel:len', ...
        'Number of bits (%d) must be a multiple of log2(M) = %d.', numel(tx_bits), bps);
end

nsym = numel(tx_bits) / bps;
w = 2.^(bps-1:-1:0);                          % MSB-first weights
sym = w * reshape(tx_bits, bps, nsym);       % 1 x nsym, values 0..M-1

gamma_s = bps * EbN0_lin;

% Signal branch: noncentral chi-square (2 dof), noncentrality 2*gamma_s.
sig = abs(sqrt(gamma_s) + (randn(1, nsym) + 1i*randn(1, nsym)) / sqrt(2)).^2;

% Competing M-1 branches: i.i.d. Exp(1); only the max matters.
max_noise = max(-log(rand(M-1, nsym)), [], 1);

sym_err = sig < max_noise;

rx_sym = sym;
e = find(sym_err);
if ~isempty(e)
    rx_sym(e) = mod(sym(e) + randi([1 M-1], 1, numel(e)), M);
end

% Unpack symbols back to MSB-first bits (division, no bitget dependency).
bmat = zeros(bps, nsym);
tmp = rx_sym;
for r = bps:-1:1
    bmat(r, :) = mod(tmp, 2);
    tmp = floor(tmp / 2);
end
rx_bits = reshape(bmat, 1, []);

end
