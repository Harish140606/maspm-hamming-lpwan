function cfg = config()
%CONFIG Central configuration/parameter file for the M-ASPM project.
%
% PURPOSE
%   Returns a single struct 'cfg' holding every physical, coding and
%   simulation parameter used across the project, so that every script
%   and function draws its numbers from one authoritative place.
%
% INPUTS
%   (none)
%
% OUTPUTS
%   cfg - struct with fields described inline below.
%
% THEORY / PAPER REFERENCE
%   A. V. Nikitin and R. L. Davidchack, "M-Ary Aggregate Spread Pulse
%   Modulation With Pulse-Shaping Power Control for Highly Scalable
%   LPWANs," IEEE Access, vol. 11, pp. 25642-25661, 2023.
%
%   Parameters below are taken directly from the paper wherever the
%   paper gives an explicit numeric value. Every parameter that the
%   paper does NOT pin down numerically is marked "IMPLEMENTATION
%   ASSUMPTION" with the reasoning documented next to it (see README.md,
%   section "Implementation Assumptions" for the consolidated list).

cfg = struct();

% -------------------------------------------------------------------
% Bit-error-rate operating points (used throughout Sections II-C, III,
% IV; e.g. Figs. 4, 7, 8, 9, 20, 22). The paper's default working pair
% is BER = 1e-4 (uncoded AWGN target) degraded to BER' = 1e-3 under
% inter-PSF interference.
% -------------------------------------------------------------------
cfg.BER      = 1e-4;   % target uncoded AWGN BER (paper default, Sec. II-B)
cfg.BER_prime = 1e-3;   % degraded BER under co-channel/inter-PSF interference (Sec. II-C)

% -------------------------------------------------------------------
% Co-PSF / co-SF channel-utilization constraint alpha, Eq. (2).
% The paper states that for pure ALOHA random access the throughput-
% maximizing value is alpha = 1/2 (Sec. I, discussion after Eq. (2)),
% and this is the value used for all numeric examples (Figs. 1, 12-22).
% -------------------------------------------------------------------
cfg.alpha = 1/2;

% -------------------------------------------------------------------
% Path-loss exponent gamma (power-law path loss, Eq. (23)).
% The paper uses gamma = 2 (free space) as the default illustrative
% value (Sec. III-B, Fig. 9) and also sweeps gamma in {2,3,4} for
% Figs. 12-19. gamma = 1.6-3.1 range is the empirical Bonn measurement
% used in Fig. 21 ([17]).
% -------------------------------------------------------------------
cfg.gamma       = 2;          % default path-loss exponent (free space)
cfg.gamma_sweep = [2, 3, 4];  % used by network_capacity comparisons

% -------------------------------------------------------------------
% M-ary orders considered. M = 16 and M = 64 are the two M-ASPM cases
% analyzed in detail throughout the paper (Figs. 1,3-5,7-9,17-22).
% LoRa is compared using its equivalent M_tilde = 2^SF for SF = 6..12
% (Sec. I, "M-ASPM Spectral efficiency"; Fig. 1, 3, 4, 5).
% -------------------------------------------------------------------
cfg.M_list      = [16, 64];
cfg.SF_list     = 6:12;                 % LoRa spreading factors used
cfg.M_tilde_SF  = 2.^cfg.SF_list;       % equivalent LoRa M for each SF
cfg.SF_benchmark   = 7;                 % benchmark range delta0(gamma), Sec. III-A
cfg.M_tilde_bench  = 2^cfg.SF_benchmark; % = 128

% -------------------------------------------------------------------
% Raised-cosine roll-off factor beta for the PSF autocorrelation
% function (ACF), 0 <= beta <= 1 (Eq. (11), Sec. II-A).
% The paper uses beta = 1/3 in the general spectral-efficiency
% discussion (Eq. (15)-(16)) and beta = 1/4 specifically for the
% Fig. 9/10/11 waveform-level simulation (Sec. III-B). Both are kept
% here; individual scripts select whichever the corresponding figure
% in the paper specifies.
% -------------------------------------------------------------------
cfg.beta_default = 1/3;   % used in Eqs. (15),(16) numeric examples
cfg.beta_fig9     = 1/4;   % used in the Fig. 9/10/11 waveform example

% -------------------------------------------------------------------
% Oversampling factor Ns = 2/(1+beta), 1 <= Ns < 2 (Sec. II-A, text
% following Eq. (13); Sec. III-B, Eq. for eta = 2*Ns*log2(M)/Np).
% -------------------------------------------------------------------
cfg.Ns_default = 2/(1+cfg.beta_default);
cfg.Ns_fig9     = 2/(1+cfg.beta_fig9);

% -------------------------------------------------------------------
% Minimum IpI for noncoherent M-ASPM: Np_min = 4*M (Sec. II-A, text
% after Eq. (16)).
% -------------------------------------------------------------------
cfg.Np_min_factor = 4;    % Np_min = Np_min_factor * M

% -------------------------------------------------------------------
% Pulse-position encoding parameters (Eq. (10), Eq. (27)).
% IMPLEMENTATION ASSUMPTION: the paper leaves the invertible mapping
% Delta_k[m] free subject only to 0 <= Delta_k[m] < Np and
% Delta_k[m] != Delta_k[l] for m != l (Sec. II, text after Eq. (10)).
% For the explicit 16-ASPM example (Eq. (27)) the paper uses a uniform
% spacing n = 4 samples between adjacent symbol positions, i.e.
% Delta_k[m] = (m-1)*n with n = 4. We adopt this same uniform-spacing
% rule for general M, which is consistent with the minimum allowed
% pulse-position offset delta_N >= 4*(M-1) given in Eq. (55).
% -------------------------------------------------------------------
cfg.pulse_position_spacing = 4;   % "n" in Eq. (27)

% -------------------------------------------------------------------
% PSF chirp design (Eq. (11)). The paper specifies only that ζ_i[k] is
% a constant-envelope "nonlinear chirp" whose autocorrelation function
% (ACF) approximates a small-time-bandwidth-product pulse such as a
% raised-cosine pulse (Sec. II-A, and Sec. II, "a good choice for the
% ACF would be ... a raised-cosine (RC) pulse [13]"); the exact phase
% law Phi_i[k] is given in the companion papers [1],[2],[7], not in
% this paper.
% IMPLEMENTATION ASSUMPTION: we synthesize Phi_i[k] using the classical
% stationary-phase / Cook-Bernfeld nonlinear-FM design rule (Cook &
% Bernfeld, "Radar Signals," 1967; also Key, Fowle & Haggarty, 1961),
% which is the standard method for building a constant-envelope chirp
% whose power spectrum approximates a prescribed magnitude shape
% (here, the square root of the raised-cosine spectrum). This is a
% widely used, standard pulse-compression technique and is documented
% here as an external reference used only to fill an implementation
% gap left open by the paper.
% -------------------------------------------------------------------
cfg.chirp_design_oversample = 8;  % samples per "chirp design" unit (numerical NLFM synthesis resolution)

% -------------------------------------------------------------------
% AWGN / SNR sweep ranges for BER curves (Figs. 3, 4).
% -------------------------------------------------------------------
cfg.EbN0_dB_range  = -5:0.25:15;   % for Fig. 3-style Eb/N0 sweeps
cfg.SNR_dB_range   = -30:0.25:0;   % for Fig. 4-style Gamma (SNR) sweeps

% -------------------------------------------------------------------
% Bit-true waveform simulation defaults (used by transmitter.m,
% matched_filter_receiver.m, plot_fig2.m, plot_fig9.m).
% -------------------------------------------------------------------
cfg.n_payload_symbols = 64;   % number of M-ary symbols per demo payload

% Fig. 2 illustration uses M = 8 explicitly (Sec. II, "In the example
% of Fig. 2, M = 8").
cfg.M_fig2 = 8;

% Fig. 9 uses noncoherent 16-ASPM (M = 16) explicitly (Sec. III-B).
cfg.M_fig9 = 16;

% Fig. 9(I): 5 nodes, SIR = -6 dB; Fig. 9(II): 9 nodes, SIR = -9 dB,
% with IpI for 9 nodes twice that for 5 nodes (Sec. III-B, Fig. 9
% caption). Minimum-range duty cycles D1 = 20% (5 nodes) and D1 = 10%
% (9 nodes); maximum duty cycle 93.8% for both (Sec. III-B text).
cfg.fig9_n_nodes        = [5, 9];
cfg.fig9_D1             = [0.20, 0.10];
cfg.fig9_Dmax_at_dmax    = 0.938;
cfg.fig9_gamma           = 2;      % free-space path loss (explicitly stated)
% IMPLEMENTATION ASSUMPTION: exact Np values (paper states eta = 1/140
% and eta = 1/280, implying particular Np through a sample-rate
% relation not fully pinned down numerically in this section) are
% replaced with clean round IpI values preserving the stated 2:1 ratio
% and the qualitative SIR/duty-cycle relationships; see README.md.
cfg.fig9_Np              = [280, 560];

% -------------------------------------------------------------------
% Fig. 20/22 hexagonal-cell areal-coverage scenario (Sec. V).
% d0 = LoRa SF6 range, dmax = LoRa SF12 range, alpha = 1/2 (stated).
% IMPLEMENTATION ASSUMPTION: the paper does not state the path-loss
% exponent gamma used for the specific numeric example in Fig. 20/22
% (unlike Fig. 21, which explicitly uses Bonn-measured gamma = 1.8 and
% 3.1). We use the paper's own default gamma = 2 (free space), the
% value used elsewhere for single, unambiguous numeric illustrations
% (e.g. Fig. 9). The resulting channel counts will therefore be
% reproductions of the same qualitative relationship described in the
% paper (roughly halving the number of partitions/duty cycles for a
% small, single-digit-percent capacity loss) rather than an exact
% match to the paper's reported "70 to 66" channel counts.
% -------------------------------------------------------------------
cfg.fig22_gamma  = 2;
cfg.fig22_SF_d0  = 6;
cfg.fig22_SF_dmax = 12;
cfg.fig22_alpha  = 1/2;

% -------------------------------------------------------------------
% NOVELTY ADDITION (beyond the paper): outer Hamming FEC.
% The paper analyses only UNCODED M-ASPM (Eq. 17, Figs. 3-4). This
% project adds a classical binary Hamming(2^m-1, 2^m-1-m) block code
% around the pulse-position bit stream (hamming_code.m,
% coded_ber_hamming.m, plot_fig_hamming.m). Hamming is the
% lowest-complexity perfect single-error-correcting code -- a natural
% fit for the paper's low-power LPWAN target -- and needs no MATLAB
% toolbox. See README.md section 9.
%   m = 3  -> Hamming(7,4),   rate 0.571
%   m = 4  -> Hamming(15,11), rate 0.733
% -------------------------------------------------------------------
cfg.hamming = struct();
cfg.hamming.enabled = true;      % gate for the Part 1b/2b novelty sections in main.m
cfg.hamming.m       = 4;         % default code used by the bit-stream demo (Hamming(15,11))
cfg.hamming.m_list  = [3, 4, 5]; % codes compared in plot_fig_hamming.m / Part 2b
                                 % Hamming(7,4), (15,11), (31,26)

% -------------------------------------------------------------------
% Random number seed for reproducibility of bit-true simulations.
% -------------------------------------------------------------------
cfg.rng_seed = 2023;

end
