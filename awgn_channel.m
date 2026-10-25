function ch = awgn_channel()
%AWGN_CHANNEL Complex-baseband additive white Gaussian noise channel.
%
% PURPOSE
%   Returns a struct 'ch' with a function that adds complex AWGN to a
%   baseband signal at a specified Eb/N0 (or, equivalently, a target
%   uncoded BER operating point via ber_calculation.m/spectral
%   efficiency), matching the AWGN-channel assumption used throughout
%   the paper's BER analysis (Sec. II-B; Figs. 3, 4).
%
% INPUTS
%   (none) -- see ch.add_noise(...) below.
%
% OUTPUTS
%   ch - struct with field 'add_noise', a function handle:
%        r = ch.add_noise(z, EbN0_lin, bits_per_symbol, Np, L)
%
% THEORY / PAPER REFERENCE
%   Section II-B: Gamma = (Eb/N0)*(fb/B), and the receiver's AWGN
%   noise power spectral density is N0 (one-sided), consistent with
%   the standard complex-baseband AWGN channel model (Proakis,
%   "Digital Communications"). Given the transmit waveform's per-pulse
%   energy Ep = L/Np * (average power over one IpI) and Eb =
%   Ep/log2(M), noise variance per complex sample is set so that the
%   resulting Eb/N0 matches the requested value once processed by the
%   matched filter (matched_filter_receiver.m) whose peak output SNR
%   is exactly Eb/N0 * (matched-filter processing gain), consistent
%   with the standard matched-filter/AWGN relationship.
%
% EQUATION REFERENCES
%   Feeds the Eq. (17)-based BER evaluation via matched_filter_receiver.m.

ch = struct();
ch.add_noise = @(z, EbN0_lin, bits_per_symbol, Np, L) ...
    local_add_noise(z, EbN0_lin, bits_per_symbol, Np, L);

end

% -------------------------------------------------------------------
function r = local_add_noise(z, EbN0_lin, bits_per_symbol, Np, L)
    % Average TRANSMIT signal power over one IpI, for constant-envelope
    % pulses of unit peak amplitude and duty cycle D=L/Np (Sec. III,
    % "the average transmit power of an M-ASPM signal with constant-
    % envelope pulses is proportional to L/Np"; the PSF itself is
    % already normalized so |zeta_hat|=1/sqrt(L), i.e. z[k] has peak
    % pulse amplitude ~1 after Eq. (12)'s sqrt(L) scaling).
    D = L / Np;

    % Energy per bit Eb = (average power * IpI duration) / bits-per-pulse
    %                   = D * Np / bits_per_symbol   (unit-power pulses)
    Eb = D * Np / bits_per_symbol;

    % Required one-sided noise PSD N0 for the target Eb/N0.
    N0 = Eb / EbN0_lin;

    % Complex AWGN with per-sample variance N0 (standard complex-
    % baseband convention: E[|n[k]|^2] = N0).
    n = sqrt(N0/2) * (randn(size(z)) + 1i*randn(size(z)));

    r = z + n;
end
