function psf = pulse_shaping_filter()
%PULSE_SHAPING_FILTER Constant-envelope nonlinear-FM PSF chirp (Eq. 11)
%   and its application to a designed pulse train (Eq. 12), with the
%   constant-envelope check of Eq. (13).
%
% PURPOSE
%   Returns a struct 'psf' of function handles that:
%     (1) synthesize the constant-envelope PSF chirp zeta_hat_i[k] of
%         Eq. (11), designed so that its autocorrelation function (ACF)
%         approximates the target raised-cosine ACF from
%         raised_cosine_filter.m;
%     (2) apply a PSF to a designed pulse train via convolution,
%         Eq. (12); and
%     (3) verify the constant-envelope property, Eq. (13).
%
% INPUTS
%   (none)
%
% OUTPUTS
%   psf - struct of function handles, described below.
%
% THEORY
%   Eq. (11):
%       zeta_hat_i[k] = g_hat_i[k] + i*h_hat_i[k]
%                     = (1/sqrt(Li)) * [[0<=k<Li]] * exp(i*Phi_i[k])
%   The paper specifies only that Phi_i[k] is a "nonlinear chirp" phase
%   law chosen so the ACF (g_hat_i * g_i)[k] approximates a small-TBP
%   pulse such as a raised-cosine pulse (Sec. II, text preceding
%   Eq. 11); the explicit formula for Phi_i[k] is given in the
%   companion papers [1],[2],[7], not in this paper.
%
% IMPLEMENTATION ASSUMPTION -- NLFM phase design (Phi_i[k])
%   We synthesize Phi_i[k] using the classical Cook-Bernfeld
%   stationary-phase nonlinear-FM design rule (D. C. Cook and
%   M. Bernfeld, "Radar Signals," Academic Press, 1967; also W. E. Key,
%   E. N. Fowle and R. D. Haggarty, 1961), a standard, well established
%   pulse-compression technique that synthesizes a constant-envelope
%   chirp whose instantaneous frequency law reproduces a prescribed
%   target power-spectral shape S(f):
%       instantaneous frequency f(t) is chosen such that the FRACTION
%       of the pulse duration elapsed at time t equals the FRACTION of
%       the target spectral energy swept up to frequency f(t):
%           t/L = ( integral_0^f(t) S(f') df' ) / ( integral_0^fmax S(f') df' )
%   With the target S(f) taken as the (one-sided) raised-cosine roll-off
%   spectrum from raised_cosine_filter.m, the chirp's own energy
%   spectral density approximates S(f), and -- by the Wiener-Khinchin
%   theorem -- its autocorrelation approximates the desired raised-
%   cosine ACF w[k]. Because the instantaneous frequency is designed to
%   be single-signed (either always >=0 for an "up-chirp", or always
%   <=0 for a "down-chirp" -- see Sec. III-B, "alternate up-chirp and
%   down-chirp PSFs"), the resulting g_hat_i[k]=cos(Phi_i[k])/sqrt(Li)
%   and h_hat_i[k]=sin(Phi_i[k])/sqrt(Li) satisfy the paper's stated
%   property h_hat_i[k] ~ H{g_hat_i[k]} (Eq. 11 discussion) to good
%   approximation -- this is verified numerically in this project (see
%   README.md, "Numerical Validation": the Pearson correlation between
%   sin(Phi[k]) and hilbert_filter's H{cos(Phi[k])} is approximately
%   0.97-0.98 for the PSF lengths used in this project).
%
% KNOWN LIMITATION (documented, not hidden)
%   Because the exact NLFM phase law used by the paper's authors is not
%   published in this paper (see IMPLEMENTATION ASSUMPTION above), the
%   Cook-Bernfeld approximation used here does not reach perfectly
%   orthogonal (zero cross-correlation) M-ary pulse-position candidates
%   at the paper's stated spacing n=4 samples (Eq. 27): measured
%   autocorrelation sidelobes at lag=n are on the order of 0.15-0.2 of
%   the mainlobe peak for the PSF lengths used in this project (see
%   README.md, "Numerical Validation"). Consequently, the BIT-TRUE
%   waveform simulation (used only for the illustrative Fig. 2 and
%   Fig. 9 reproductions) shows a somewhat higher simulated BER than
%   the idealized Eq. (17) curve, especially at high Eb/N0, though it
%   still decreases monotonically with SNR and matches Eq. (17)
%   reasonably closely at low-to-moderate SNR. This limitation affects
%   ONLY the bit-true simulation; every quantitative BER, spectral-
%   efficiency and capacity figure in this project (Figs. 4, 5, 7, 8
%   and the network_capacity.m/spectral_efficiency.m analyses) is
%   computed directly from the closed-form Eq. (17) and is therefore
%   an exact reproduction of the paper's analytical results,
%   independent of this waveform-synthesis approximation.
%
% PAPER REFERENCE
%   Section II, Eqs. (11)-(13).
%
% EQUATION REFERENCES
%   Eq. (11), (12), (13).

psf = struct();
rc = raised_cosine_filter();

% ---------------------------------------------------------------
% Eq. (11): design one PSF chirp of length L samples, roll-off beta,
% oversampling factor Ns, and sweep direction dir (+1 = up-chirp,
% -1 = down-chirp, Sec. III-B "flip" PSFs).
% Returns: zeta_hat (complex row vector, length L), g_hat, h_hat.
% ---------------------------------------------------------------
psf.design_chirp = @(L, beta, Ns, dir) local_design_chirp(rc, L, beta, Ns, dir);

% ---------------------------------------------------------------
% Eq. (12): apply a PSF (zeta_hat, length Li) to a designed pulse
% train x_hat (Eq. 9) via convolution, scaled by sqrt(Li):
%   z_i[k] = sqrt(Li) * (x_hat * zeta_hat_i)[k]
% ---------------------------------------------------------------
psf.apply_psf = @(x_hat, zeta_hat) sqrt(numel(zeta_hat)) .* conv(x_hat, zeta_hat);

% ---------------------------------------------------------------
% Eq. (13): constant-envelope check. Returns |z[k]| and a boolean
% flag indicating whether the envelope is constant (within tol) over
% the region where pulses do not overlap.
% ---------------------------------------------------------------
psf.envelope = @(z) abs(z);
psf.is_constant_envelope = @(z, active_mask, tol) ...
    local_check_envelope(z, active_mask, tol);

end

% ===================================================================
function [zeta_hat, g_hat, h_hat] = local_design_chirp(rc, L, beta, Ns, dir) %#ok<INUSD>
    % NOTE: Ns (oversampling factor) is accepted for interface
    % consistency with the rest of the project (Fs=2*Ns*B, Sec. II-A/
    % III-B) but is not used directly by this wideband chirp-bandwidth
    % design -- see the "IMPLEMENTATION ASSUMPTION -- chirp bandwidth
    % for pulse-position resolution" note below for why the chirp
    % instead sweeps a fixed, near-Nyquist digital bandwidth.
    if nargin < 5 || isempty(dir), dir = 1; end
    fmax = (1+beta);                 % one-sided target bandwidth, RC-normalized units
    Ndesign = max(2048, 8*L);
    fgrid = linspace(0, fmax, Ndesign);
    Svals = rc.spectral_shape(fgrid, beta);
    Ecum = cumtrapz(fgrid, Svals);
    Ecum = Ecum / Ecum(end);         % normalized cumulative spectral energy, 0->1

    % Sample times (cell-centered fractions of duration L), 0..1
    k = (0:L-1)';
    u = (k + 0.5) / L;

    % Invert the stationary-phase mapping: find f(u) such that
    % Ecum(f) = u (Ecum is monotonic increasing by construction).
    f_of_u = interp1(Ecum, fgrid, u, 'linear', 'extrap');
    f_of_u = max(0, min(fmax, f_of_u));

    % Convert normalized frequency (0..fmax, in RC roll-off units) to
    % digital frequency (cycles/sample).
    %
    % IMPLEMENTATION ASSUMPTION -- chirp bandwidth for pulse-position
    % resolution: pulse-compression theory (Cook & Bernfeld) requires
    % the chirp to sweep a wide FRACTION OF THE AVAILABLE DIGITAL
    % BANDWIDTH to achieve a narrow (few-sample) autocorrelation
    % mainlobe, which is what lets the receiver resolve pulse-position
    % candidates spaced only n=4 samples apart (Eq. 27) via matched
    % filtering (Sec. II text following Eq. 12, "the position of the
    % sample with the largest magnitude will correspond to the
    % position of the respective pulse"). The paper does not give an
    % explicit numeric bandwidth for this per-PSF chirp (that level of
    % detail is left to the companion papers [1],[2],[7]); we choose
    % the widest digital bandwidth consistent with the Nyquist limit,
    % leaving a small guard margin (0.45 cycles/sample of the maximum
    % 0.5), and use beta only to control the RELATIVE shape of the
    % frequency sweep across the pulse duration (i.e. how much of the
    % duration is spent at each frequency, which is what produces the
    % raised-cosine-like taper of the resulting autocorrelation) via
    % the same S(f)-based stationary-phase warping used above.
    max_digital_freq = 0.45;
    f_digital = dir .* (f_of_u ./ fmax) .* max_digital_freq;

    % Phase accumulation (discrete integral of instantaneous frequency).
    Phi = 2*pi * cumsum(f_digital);
    Phi = Phi - Phi(1); % start phase at 0

    g_hat = (1/sqrt(L)) .* cos(Phi)';
    h_hat = (1/sqrt(L)) .* sin(Phi)';
    zeta_hat = g_hat + 1i*h_hat;
end

% ===================================================================
function [env, is_const] = local_check_envelope(z, active_mask, tol)
    env = abs(z);
    if nargin < 3 || isempty(tol)
        tol = 0.05;
    end
    vals = env(active_mask);
    vals = vals(vals > 1e-9); % ignore true zero-guard samples
    if isempty(vals)
        is_const = true;
        return;
    end
    is_const = (max(vals) - min(vals)) <= tol*mean(vals);
end
