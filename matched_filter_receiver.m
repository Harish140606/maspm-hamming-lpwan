function mfr = matched_filter_receiver()
%MATCHED_FILTER_RECEIVER Noncoherent M-ASPM matched-filter detection (Fig. 2-III).
%
% PURPOSE
%   Returns a struct 'mfr' with a function implementing the receiver
%   chain of Fig. 2(III): quadrature demodulation (already implicit,
%   since we work at complex baseband -- see transmitter.m), matched
%   filtering, and noncoherent combining into the detection statistic
%   y_nc[k]^2.
%
% INPUTS
%   (none) -- see mfr.detect(...) below.
%
% OUTPUTS
%   mfr - struct with field 'detect', a function handle:
%         y2 = mfr.detect(r, zeta_hat)
%
% THEORY
%   Fig. 2(III) explicitly gives the noncoherent detection statistic:
%       y_nc^2[k] = (I*g+Q*h)^2[k] + (Q*g-I*h)^2[k] = |(I+iQ)*(g-i h)|^2[k]
%   where r[k]=I[k]+i*Q[k] is the (baseband) received signal and
%   g[k]-i*h[k] = zeta_i[k] is the PSF's MATCHED FILTER, defined in
%   Sec. II as zeta_i[k] = zeta_hat_i*[-k] (time-reversed complex
%   conjugate of the transmit PSF zeta_hat_i[k]). Substituting this
%   definition shows the two expressions for y_nc^2 are IDENTICAL:
%       (I+iQ)*(g-ih) = r * zeta_i = r * conj(fliplr(zeta_hat_i))
%   so y_nc^2[k] = | (r * conj(fliplr(zeta_hat_i)))[k] |^2, i.e. simply
%   the squared magnitude of r[k] passed through the PSF's STANDARD
%   MATCHED FILTER. This function implements that equivalent, simpler
%   form directly (verified algebraically above; no change to the
%   figure's formula, only a simplification of its evaluation).
%
% PAPER REFERENCE
%   Section II, Fig. 2(III); Section II text defining zeta_i[k] as the
%   matched filter of zeta_hat_i[k].
%
% EQUATION REFERENCES
%   Fig. 2(III) noncoherent-detection formula; Section II definition
%   zeta_i[k] = zeta_hat_i*[-k].

mfr = struct();
mfr.matched_filter = @(zeta_hat) conj(fliplr(zeta_hat));
mfr.detect = @(r, zeta_hat) local_detect(r, zeta_hat);

end

% -------------------------------------------------------------------
function y2 = local_detect(r, zeta_hat)
    zeta = conj(fliplr(zeta_hat));   % matched filter, zeta_i[k]=zeta_hat_i*[-k]
    y = conv(r, zeta);
    y2 = abs(y).^2;
end
