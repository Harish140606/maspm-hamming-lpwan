function tx = transmitter()
%TRANSMITTER Full M-ASPM transmit chain (Fig. 2, blocks I-II).
%
% PURPOSE
%   Returns a struct 'tx' with a single high-level function that runs
%   the complete transmitter signal path described in Fig. 2 of the
%   paper: bits -> M-ASPM encoding (Eq. 9-10) -> PSF design (Eq. 11)
%   -> pulse shaping (Eq. 12) -> constant-envelope complex baseband
%   waveform z[k] = x_g[k] + i*x_h[k], ready for (conceptual) D/A
%   conversion and single-sideband upconversion to a carrier.
%
% INPUTS
%   (none) -- see tx.run_chain(...) below.
%
% OUTPUTS
%   tx - struct with field 'run_chain', a function handle:
%        [z, info] = tx.run_chain(bits, M, Np, L, beta, Ns, dir, n_spacing, deltaN)
%
% THEORY / PAPER REFERENCE
%   Section II, Fig. 2. We work throughout in the COMPLEX BASEBAND
%   EQUIVALENT representation of the single-sideband passband signal
%   Re(z(t))*sin(2*pi*fc*t) + Im(z(t))*cos(2*pi*fc*t) (Sec. II, text
%   following Eq. 12). This is the STANDARD baseband-equivalent
%   channel model used throughout digital communications simulation
%   (e.g. J. G. Proakis, "Digital Communications") and lets us skip
%   the (physically trivial, for simulation purposes) RF carrier
%   mixing stage, matching Fig. 2's D/A block conceptually while
%   keeping the simulation carrier-frequency-independent.
%   IMPLEMENTATION ASSUMPTION: explicit carrier upconversion/D-A/A-D
%   conversion are not modeled; this is standard practice whenever a
%   channel (AWGN, multi-user interference) is itself modeled at
%   complex baseband, since upconversion/downconversion by a common
%   carrier does not affect the SNR/SIR-determined performance metrics
%   this project reproduces (BER, spectral efficiency, capacity).
%
% EQUATION REFERENCES
%   Eq. (9), (10), (11), (12), (13).

tx = struct();
tx.run_chain = @(bits, M, Np, L, beta, Ns, dir, n_spacing, deltaN) ...
    local_run_chain(bits, M, Np, L, beta, Ns, dir, n_spacing, deltaN);

end

% -------------------------------------------------------------------
function [z, info] = local_run_chain(bits, M, Np, L, beta, Ns, dir, n_spacing, deltaN)

    if nargin < 9 || isempty(deltaN),   deltaN = 0; end
    if nargin < 8 || isempty(n_spacing), n_spacing = 4; end
    if nargin < 7 || isempty(dir),       dir = 1; end

    % --- Eq. (9)-(10): bits -> designed pulse train x_hat[k] ---
    [kj, mj, signs, x_hat, N_total] = aspm_encoder(bits, M, Np, deltaN, n_spacing); %#ok<ASGLU>

    % --- Eq. (11): design the constant-envelope PSF chirp ---
    psf = pulse_shaping_filter();
    [zeta_hat, g_hat, h_hat] = psf.design_chirp(L, beta, Ns, dir); %#ok<ASGLU>

    % --- Eq. (12): apply the PSF to the designed pulse train ---
    z = psf.apply_psf(x_hat, zeta_hat);

    % --- Eq. (13): constant-envelope verification (non-overlapping
    % pulses only, i.e. away from the very first/last pulse's filter
    % ramp-up/down at the signal edges) ---
    active_region = false(size(z));
    guard = L;
    if numel(z) > 2*guard
        active_region((guard+1):(numel(z)-guard)) = true;
    end
    [env, is_const] = psf.is_constant_envelope(z, active_region, 0.08); %#ok<ASGLU>

    info = struct();
    info.kj = kj;
    info.mj = mj;
    info.signs = signs;
    info.x_hat = x_hat;
    info.zeta_hat = zeta_hat;
    info.g_hat = g_hat;
    info.h_hat = h_hat;
    info.N_total = N_total;
    info.Np = Np;
    info.L = L;
    info.M = M;
    info.constant_envelope_ok = is_const;
    info.envelope = env;
end
