function pl = path_loss_model()
%PATH_LOSS_MODEL Power-law path loss and range<->SNR/range-matching tools.
%
% PURPOSE
%   Returns a struct 'pl' of function handles implementing the paper's
%   power-law path-loss model (Eq. 23) and the associated physical
%   range formula for M-ASPM (Eq. 26), including the LoRa-SF7 benchmark
%   range delta0(gamma) (Sec. III-A) used as the common x-axis reference
%   ("range relative to LoRa w/ SF = 7") in Figs. 1, 5, 8, 20-22.
%
% INPUTS
%   (none)
%
% OUTPUTS
%   pl - struct of function handles, described below.
%
% THEORY / PAPER REFERENCE
%   Section III, Eq. (23) (power-law path loss); Section III-A,
%   Eqs. (25)-(26) (physical range referenced to LoRa SF=7 benchmark).
%
% EQUATION REFERENCES
%   Eq. (23): P(d) = P(d0) * (d0/d)^gamma
%   Eq. (25): (d1/d2)^gamma = (D1*N1)/(D2*N2)
%   Eq. (26): d = (D*eta0/eta)^(1/gamma) * delta0(gamma)

pl = struct();
se = spectral_efficiency();

% Eq. (23): received power at distance d, given reference power P(d0)
% at reference distance d0, and path-loss exponent gamma.
pl.received_power = @(P_d0, d0, d, gamma) P_d0 .* (d0./d).^gamma;

% Eq. (23) inverted: since AWGN SNR Gamma is proportional to received
% power, Gamma(d)/Gamma(dref) = (dref/d)^gamma for a fixed transmit
% power and noise floor. This is the relation used throughout the
% paper (Sec. II-B, "for the power-law path loss Gamma ~ d^-gamma")
% to convert between range ratios and SNR (Gamma) ratios.
pl.snr_ratio_from_range_ratio = @(d, dref, gamma) (dref./d).^gamma;
pl.range_ratio_from_snr_ratio = @(Gamma_ratio, gamma) Gamma_ratio.^(-1./gamma);

% -----------------------------------------------------------------
% Eq. (18)/(26): "equal received-SNR-at-target-BER" range matching.
%
% For any modulation with spectral efficiency eta (M-ASPM at a given
% Np, or LoRa at a given SF, eta_tilde = log2(Mx)/Mx), operating with
% pulse duty cycle D (D=1 for full constant-power, spread-spectrum
% IpI-controlled range as in Sec. II-B; D<1 for duty-cycle range
% control as in Sec. III), the SNR required to reach a target BER is
%     Gamma_required = eta * required_EbN0(Mx, BER)   [since rho=Gamma/eta]
% and for power-law path loss with a reference transmitter of the same
% peak power, Gamma ~ D * d^-gamma (duty cycle D scales the *average*
% transmit power for constant-envelope pulses, Sec. III, first
% paragraph: "average transmit power is proportional to D"). The range
% at which this modulation achieves the target BER, relative to any
% reference range d_ref (itself defined by its own eta_ref, M_ref,
% D_ref), follows from equal average-received-power-normalized-to-Gamma
% at each transmitter's own operating point:
%     (d/d_ref)^gamma = (D*Gamma_ref) / (D_ref*Gamma)
% This directly generalizes Eq. (26) (which is the special case
% D_ref = 1, M_ref = 128 i.e. LoRa SF=7, giving d_ref = delta0(gamma)).
% -----------------------------------------------------------------
pl.required_gamma_snr = @(eta, Mx, BER, D) ...
    eta .* se.required_EbN0(Mx, BER) ./ D;

pl.range_relative_to_reference = @(eta, Mx, BER, D, eta_ref, Mref, BERref, Dref, gamma) ...
    local_range_ratio(pl, eta, Mx, BER, D, eta_ref, Mref, BERref, Dref, gamma);

% -----------------------------------------------------------------
% delta0(gamma): LoRa SF=7 (M_tilde = 128) benchmark range, expressed
% as a range RATIO relative to itself, i.e. delta0(gamma)/delta0(gamma)
% = 1 by definition (Sec. III-A). Provided for clarity/symmetry with
% the paper's notation; range_relative_to_reference() should be called
% with eta_ref = lora_eta(128), Mref = 128, D_ref = 1 to place any
% other operating point on the delta0(gamma)-normalized axis used in
% Figs. 1, 5, 8, 20-22.
% -----------------------------------------------------------------
pl.delta0_ratio = 1;

end

% -------------------------------------------------------------
function ratio = local_range_ratio(pl, eta, Mx, BER, D, eta_ref, Mref, BERref, Dref, gamma)
    Gamma_x   = pl.required_gamma_snr(eta, Mx, BER, D);
    Gamma_ref = pl.required_gamma_snr(eta_ref, Mref, BERref, Dref);
    % (d/d_ref)^gamma = Gamma_ref / Gamma_x  (higher required SNR at
    % the SAME transmit power/noise floor means a SHORTER usable range)
    ratio = (Gamma_ref ./ Gamma_x) .^ (1./gamma);
end
