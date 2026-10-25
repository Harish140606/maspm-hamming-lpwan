function se = spectral_efficiency()
%SPECTRAL_EFFICIENCY LoRa and M-ASPM spectral-efficiency relations.
%
% PURPOSE
%   Returns a struct 'se' of function handles implementing Eqs. (14),
%   (15), (16) and (18) of the paper: LoRa's spectral efficiency, the
%   M-ASPM spectral efficiency as a function of (M, Np, beta), the
%   maximum M-ASPM spectral efficiency, and the Eb/N0 required to hit a
%   target BER for a given M (the numeric inversion needed for Eq. 18
%   and, downstream, Eq. 26).
%
% INPUTS
%   (none)
%
% OUTPUTS
%   se - struct of function handles, described field by field below.
%
% THEORY / PAPER REFERENCE
%   Section II-A ("M-ASPM Spectral Efficiency"), Eqs. (14)-(16);
%   Section II-B, Eq. (18); Section III-A, Eq. (26).
%
% EQUATION REFERENCES
%   Eq. (14): eta_tilde(M_tilde) = log2(M_tilde) / M_tilde        (LoRa)
%   Eq. (15): eta(M,Np) = 4*log2(M) / ((1+beta)*Np)                (M-ASPM)
%   Eq. (16): eta_max(M) = eta_tilde(M) / (1+beta), at Np = 4*M
%   Eq. (18): Pb(Gamma;M,Np) = Pb_tilde(Gamma;M_tilde) = BER
%             -> solved here as: required Eb/N0 for M at target BER.

se = struct();

% Eq. (14): LoRa spectral efficiency for spreading factor equivalent M_tilde = 2^SF.
se.lora_eta = @(M_tilde) log2(M_tilde) ./ M_tilde;

% Eq. (15): M-ASPM spectral efficiency for M-ary order M, IpI Np, RC roll-off beta.
se.maspm_eta = @(M, Np, beta) 4*log2(M) ./ ((1+beta) .* Np);

% Inverse of Eq. (15): IpI Np required to realize a target spectral efficiency eta.
se.maspm_Np_from_eta = @(M, eta, beta) 4*log2(M) ./ ((1+beta) .* eta);

% Eq. (16): maximum M-ASPM spectral efficiency, achieved at Np = 4M.
se.maspm_eta_max = @(M, beta) se.lora_eta(M) ./ (1+beta);

% Required Eb/N0 (linear) to achieve a target BER for M-ary noncoherent
% detection, i.e. numeric inversion of Eq. (17) used by Eq. (18).
% ber_calculation(EbN0,M) is monotonically decreasing in EbN0, so a
% bracketed bisection search is well posed.
se.required_EbN0 = @(M, target_BER) local_required_EbN0(M, target_BER);

end

% ------------------------------------------------------------------
function EbN0 = local_required_EbN0(M, target_BER)
% LOCAL_REQUIRED_EBN0 Solve ber_calculation(EbN0,M) = target_BER for
%   EbN0 (linear), via bisection (Eq. 17 inverted, used in Eq. 18).
    h = helper_functions();
    f = @(x) ber_calculation(x, M);
    lo = 1e-6;   % Eb/N0 -> 0 gives Pb -> (M-1)/(2(M-1)) ~ 0.5 (upper bound on Pb)
    hi = 1e4;    % Eb/N0 large gives Pb -> 0
    % Expand hi if needed so that f(hi) < target_BER (guarantees bracketing
    % since f is monotonically decreasing).
    iter_guard = 0;
    while f(hi) > target_BER && iter_guard < 60
        hi = hi * 2;
        iter_guard = iter_guard + 1;
    end
    EbN0 = h.bisection_solve(f, target_BER, lo, hi, 1e-10);
end
