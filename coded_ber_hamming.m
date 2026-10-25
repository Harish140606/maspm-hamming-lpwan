function Pb_out = coded_ber_hamming(EbN0_lin, M, m)
%CODED_BER_HAMMING Post-decoding BER of an outer Hamming(n,k) code over the
%   noncoherent M-ary M-ASPM/AWGN channel of Eq. (17).
%
% NOVELTY ADDITION (not part of Nikitin & Davidchack, IEEE Access 2023).
%   Extends the paper's uncoded BER machinery (ber_calculation.m, Eq. 17)
%   with a hard-decision (bounded-distance) block-decoding analysis for the
%   Hamming codes implemented in hamming_code.m.
%
% INPUTS
%   EbN0_lin - Eb/N0 (LINEAR) per INFORMATION bit. Scalar or array.
%   M        - M-ASPM alphabet size (power of two).
%   m        - Hamming parameter: n = 2^m-1, k = n-m (default m = 3 -> (7,4)).
%
% OUTPUT
%   Pb_out   - approximate information-bit error probability after decoding,
%              same size as EbN0_lin.
%
% MODEL
%   1. Code rate R = k/n. For a fixed energy Eb per INFORMATION bit, each
%      transmitted (coded) bit carries R*Eb, so the raw channel bit-error
%      probability is
%          p = ber_calculation(EbN0_lin * R, M)              (Eq. 17)
%      i.e. the rate loss is charged honestly against the code.
%   2. A Hamming code is a perfect single-error-correcting code (t = 1):
%      a block is decoded correctly iff it contains 0 or 1 channel bit
%      errors; with >= 2 errors the decoder miscorrects.
%      NOTE: an M-ary symbol error corrupts up to log2(M) adjacent coded
%      bits at once. This model (and the classical block-code analysis it
%      uses) assumes those errors are spread across codewords by
%      interleaving (block_interleaver.m in the simulations); without
%      interleaving a short Hamming code sees bursts beyond its t = 1
%      reach and gains nothing.
%   3. The output BER is approximated by the standard t = 1 expression
%          Pb_out ~= (1/n) * sum_{i=2}^{n} i * C(n,i) * p^i * (1-p)^(n-i)
%      (an uncorrected block carrying i channel errors is taken to produce
%      ~ i information-bit errors, normalised by the block length n). This
%      is the classical textbook approximation -- e.g. B. Sklar, "Digital
%      Communications," 2nd ed., Eq. (6.46); J. G. Proakis, "Digital
%      Communications," ch. 8. It slightly overestimates Pb_out (safe
%      direction) and is exact in the leading p^2 term that dominates at
%      the BER = 1e-3..1e-4 operating points the paper cares about.
%
% This is the coded analogue of ber_calculation.m and is used by
% plot_fig_hamming.m and by main.m (Parts 1b, 2b).

if nargin < 3 || isempty(m), m = 3; end
n = 2^m - 1;
k = n - m;
R = k / n;

EbN0_lin = EbN0_lin(:);
p = ber_calculation(EbN0_lin * R, M);      % raw channel BER at coded-bit energy
p = p(:);

i = (2:n)';
logC = gammaln(n + 1) - gammaln(i + 1) - gammaln(n - i + 1);   % log C(n,i)

Pb_out = zeros(size(p));
for idx = 1:numel(p)
    pp = p(idx);
    if pp <= 0
        Pb_out(idx) = 0;
        continue;
    end
    if pp >= 1
        Pb_out(idx) = 1;
        continue;
    end
    logterm = logC + i .* log(pp) + (n - i) .* log1p(-pp);
    Pb_out(idx) = sum(i .* exp(logterm)) / n;
end

Pb_out = min(max(Pb_out, 0), 1);

end
