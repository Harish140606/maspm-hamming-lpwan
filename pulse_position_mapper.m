function ppm = pulse_position_mapper()
%PULSE_POSITION_MAPPER Invertible symbol-to-offset map Delta_k[m] (Eq. 10).
%
% PURPOSE
%   Returns a struct 'ppm' of function handles implementing the
%   integer-valued, invertible pulse-position offset function
%   Delta_k[m] used in Eq. (10), together with its inverse (used by
%   the receiver, aspm_decoder.m, to recover the transmitted symbol
%   from a detected pulse position).
%
% INPUTS
%   (none)
%
% OUTPUTS
%   ppm - struct of function handles, described below.
%
% THEORY / PAPER REFERENCE
%   Eq. (10): kj = j*Np + Delta_N + Delta_k[mj], where Delta_k[m] is
%   "an integer-valued invertible function such that 0 <= Delta_k[m] <
%   Np and Delta_k[m] != Delta_k[l] for m != l" (Sec. II, text
%   following Eq. 10). The paper leaves the specific mapping free
%   subject to these constraints, except for the fully worked 16-ASPM
%   example of Eq. (27):
%       x_hat[k] = sum_j [[ k = j*Np + (8*a_j+4*b_j+2*c_j+d_j)*n ]],  n = 4
%   which is exactly Delta_k[m] = (m-1)*n with n = 4 (m-1 running over
%   the 16 possible 4-bit values 0..15).
%
% IMPLEMENTATION ASSUMPTION
%   We adopt this same uniform-spacing rule, Delta_k[m] = (m-1)*n, for
%   general M, with n = cfg.pulse_position_spacing (default 4, matching
%   Eq. 27). This satisfies the invertibility/non-collision requirement
%   for any n >= 1, and is consistent with the minimum allowed maximum
%   pulse-position offset delta_N >= 4*(M-1) given in Eq. (55) (which,
%   for the smallest allowed delta_N, is exactly (M-1)*4).
%
% EQUATION REFERENCES
%   Eq. (10), (27), (55).

ppm = struct();

% Forward map: symbol index m (1..M) -> offset Delta_k[m] (0..(M-1)*n).
% m is 1-indexed here (m=1 corresponds to the all-zero bit pattern).
ppm.forward = @(m, n) (m-1) .* n;

% Inverse map: offset -> symbol index m (1..M).
ppm.inverse = @(delta_k, n) round(delta_k ./ n) + 1;

% Convenience: minimum required Np for order M and spacing n so that
% Delta_k[M] = (M-1)*n < Np (Eq. 10 constraint 0<=Delta_k[m]<Np).
ppm.min_Np = @(M, n) (M-1)*n + 1;

end
