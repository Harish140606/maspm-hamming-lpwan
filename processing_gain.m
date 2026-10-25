function pg = processing_gain()
%PROCESSING_GAIN SIR-margin, processing-gain and IpI/capacity relations.
%
% PURPOSE
%   Returns a struct 'pg' of function handles implementing the paper's
%   SIR-margin derivation (Eqs. 19-21) and the basic capacity-vs-range
%   relations from the introductory toy example (Eqs. 1-5).
%
% INPUTS
%   (none)
%
% OUTPUTS
%   pg - struct of function handles, described below.
%
% THEORY / PAPER REFERENCE
%   Section I, Eqs. (1)-(5) (capacity C(d), number of PSF channels
%   m(d), Clim); Section II-C, Eqs. (19)-(21) (SIR margin derivation
%   and its proportionality to IpI/range).
%
% EQUATION REFERENCES
%   Eq. (1):  C(d) = m(d)*C1(d)
%   Eq. (2):  SIR^-1 = (m-1)*alpha <= Delta_SIR
%   Eq. (3):  m(d) = 1 + Delta_SIR(d)/alpha
%   Eq. (4):  C(d) = ceil(1+Delta_SIR(d)/alpha)*C1(d) = C1(d)+Clim
%   Eq. (5):  C(d->inf) = Clim = const
%   Eq. (19): Pb(Gamma/eta;M) = BER = Pb(Gamma'/eta;M)*(BER/BER')
%   Eq. (20): Delta_SIR = 1/Gamma' - 1/Gamma = eta^-1 * const
%   Eq. (21): Delta_SIR ~ eta^-1 ~ Np ~ d^gamma

pg = struct();
se = spectral_efficiency();

% ---------------------------------------------------------------
% Eq. (19)-(20): SIR-margin proportionality constant
%   const = 1/rho2 - 1/rho1
% where rho1 = Eb/N0 required to reach BER at order M (the "clean"
% AWGN operating point), and rho2 = Eb/N0 required to reach the
% degraded BER' at order M (the operating point once inter-PSF
% interference has consumed some SIR margin). Both are obtained by
% inverting Eq. (17) via spectral_efficiency().required_EbN0.
% ---------------------------------------------------------------
pg.sir_margin_const = @(M, BER, BER_prime) local_sir_const(se, M, BER, BER_prime);

% ---------------------------------------------------------------
% Eq. (20)-(21): SIR margin as a function of spectral efficiency eta.
%   Delta_SIR(eta) = const / eta
% ---------------------------------------------------------------
pg.delta_sir_of_eta = @(eta, M, BER, BER_prime) ...
    local_sir_const(se, M, BER, BER_prime) ./ eta;

% ---------------------------------------------------------------
% Eq. (20)-(21) combined with Eq. (15): SIR margin as a function of
% IpI Np directly:  Delta_SIR(Np) = const*(1+beta)*Np / (4*log2(M))
% ---------------------------------------------------------------
pg.delta_sir_of_Np = @(Np, M, beta, BER, BER_prime) ...
    local_sir_const(se, M, BER, BER_prime) .* (1+beta) .* Np ./ (4*log2(M));

% ---------------------------------------------------------------
% N0dB: the IpI Np for which Delta_SIR(Np) = 1 (0 dB SIR margin),
% i.e. N0dB = Delta_SIR^-1(1) (Sec. IV, text following Eq. 34).
%   From Delta_SIR(Np) = const*(1+beta)*Np/(4 log2 M) = 1
%   =>  N0dB = 4*log2(M) / ((1+beta)*const)
% ---------------------------------------------------------------
pg.N0dB = @(M, beta, BER, BER_prime) ...
    4*log2(M) ./ ((1+beta) .* local_sir_const(se, M, BER, BER_prime));

% ---------------------------------------------------------------
% Eq. (3): number of usable PSF channels at range d (continuous,
% non-ceiling'd form so that Eq. (4)'s Clim limit can be exposed).
% ---------------------------------------------------------------
pg.m_of_delta_sir = @(delta_sir, alpha) 1 + delta_sir ./ alpha;

% ---------------------------------------------------------------
% Eq. (4): total capacity given C1(d) (nodes per single PSF channel)
% and the SIR margin at that range/IpI.
% ---------------------------------------------------------------
pg.capacity_C = @(C1, delta_sir, alpha) ceil(1 + delta_sir./alpha) .* C1;

% ---------------------------------------------------------------
% Eq. (5)/(46): Clim, the asymptotic capacity limit common to all
% ranges/IpIs -- Clim = Delta_SIR(Np)*C1(Np)/alpha = Np/(alpha*N0dB)*C1(Np).
% ---------------------------------------------------------------
pg.Clim = @(C1, delta_sir, alpha) delta_sir .* C1 ./ alpha;

end

% --------------------------------------------------------------
function c = local_sir_const(se, M, BER, BER_prime)
% LOCAL_SIR_CONST Eq. (19)-(20) proportionality constant
%   const = 1/rho2 - 1/rho1, rho1 = required_EbN0(M,BER),
%   rho2 = required_EbN0(M,BER_prime), BER_prime > BER.
    rho1 = se.required_EbN0(M, BER);
    rho2 = se.required_EbN0(M, BER_prime);
    c = 1./rho2 - 1./rho1;
end
