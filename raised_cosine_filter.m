function rc = raised_cosine_filter()
%RAISED_COSINE_FILTER Target raised-cosine ACF/spectrum for the M-ASPM PSF.
%
% PURPOSE
%   Returns a struct 'rc' of function handles that build the TARGET
%   autocorrelation function (ACF) referenced in Eq. (11) -- "a good
%   choice for the ACF would be...a raised-cosine (RC) pulse [13]"
%   (Sec. II, discussion preceding Eq. 11) -- and its associated
%   frequency-domain (roll-off) magnitude shape, which is what the
%   nonlinear-FM chirp PSF in pulse_shaping_filter.m is designed to
%   reproduce.
%
% INPUTS
%   (none)
%
% OUTPUTS
%   rc - struct of function handles, described below.
%
% THEORY
%   The paper's own roll-off factor beta and oversampling relation
%   Fs = 2*Ns*B, Ns = 2/(1+beta) (Sec. II-A, text following Eq. 13;
%   Sec. III-B) are exactly the parameters of the classical Nyquist
%   raised-cosine pulse (J. G. Proakis and D. G. Manolakis, "Digital
%   Signal Processing," 4th ed. -- the paper's own reference [13]).
%   By the Wiener-Khinchin theorem, if a chirp's autocorrelation
%   function is to equal a time-domain raised-cosine pulse w[k], then
%   the chirp's own energy spectral density |Z(f)|^2 must equal the
%   Fourier transform of w[k], i.e. the classical raised-cosine
%   FREQUENCY roll-off shape S(f) (a smooth, compact-support spectral
%   window with flat passband |f|<=(1-beta)/2 and a raised-cosine
%   taper out to |f|=(1+beta)/2). This S(f) is exactly the target
%   spectral magnitude used by pulse_shaping_filter.m's nonlinear-FM
%   (NLFM) chirp synthesis (Cook-Bernfeld stationary-phase method).
%
% PAPER REFERENCE
%   Section II, discussion preceding Eq. (11); Section II-A
%   (Eqs. 15-16, sample-rate/roll-off relations); Ref. [13]
%   (Proakis & Manolakis).
%
% EQUATION REFERENCES
%   Eq. (11) (the ACF w[k] that zeta_i[k] is designed to approximate).

rc = struct();

% Time-domain target ACF w[k]: the standard Nyquist raised-cosine
% pulse, via MATLAB's Communications Toolbox rcosdesign (an accepted
% standard implementation of the classical RC pulse referenced by the
% paper's own citation [13]).
%   beta - roll-off factor, 0<=beta<=1 (Eq. 11 context)
%   span - pulse truncation, in symbols, each of 'sps' samples
%   sps  - samples per symbol (oversampling)
rc.target_acf = @(beta, span, sps) rcosdesign(beta, span, sps, 'normal');

% Frequency-domain raised-cosine roll-off magnitude shape S(f), on a
% normalized frequency axis f in [-1,1] representing [-B,B] (baseband,
% single-sided bandwidth B at f=1-beta..1 tapering to 0 at f=1+... );
% here we normalize so the flat passband is |f|<=(1-beta), tapering to
% 0 at |f|=(1+beta), matching the standard normalized RC spectral
% shape (Proakis, Eq. 9.2-16-style definition) with unit "symbol rate".
rc.spectral_shape = @(f, beta) local_rc_spectrum(f, beta);

end

% -------------------------------------------------------------
function S = local_rc_spectrum(f, beta)
    absf = abs(f);
    S = zeros(size(f));
    if beta == 0
        S(absf <= 1) = 1;
        return;
    end
    flat = absf <= (1-beta);
    taper = (absf > (1-beta)) & (absf <= (1+beta));
    S(flat) = 1;
    S(taper) = 0.5*(1 + cos( (pi/(2*beta)) .* (absf(taper) - (1-beta)) ));
end
