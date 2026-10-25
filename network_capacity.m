function nc = network_capacity()
%NETWORK_CAPACITY Areal-coverage capacity analysis for M-ASPM networks.
%
% PURPOSE
%   Returns a struct 'nc' of function handles implementing:
%     (a) the two-range "toy example" of areal coverage (Eqs. 6-8),
%     (b) the general equal-power vs. equal-IpI PSF-channel-partition
%         machinery for a uniform areal node density within a
%         hexagonal cell (Eqs. 28-41), including the "flip PSF pair"
%         variant (Eq. 41) used to reproduce Fig. 22.
%
% INPUTS
%   (none)
%
% OUTPUTS
%   nc - struct of function handles, described field by field below.
%
% THEORY / PAPER REFERENCE
%   Section I-A (toy example, Eqs. 6-8); Section III-C (node density
%   functions, Eqs. 28-32); Section IV (equal-power vs. equal-IpI
%   partitioning, Eqs. 33-53); Section IV-A-2 and Fig. 22 (flip-PSF
%   pairs, Eq. 41).
%
% IMPLEMENTATION ASSUMPTION -- hexagonal areal geometry
%   The paper describes (Sec. IV) "uniform areal coverage beyond range
%   d0 within a hexagonal cell with circumradius dmax" but does not
%   give the explicit closed-form radial node-density function phi(r)
%   for this geometry (only the general definition, Eq. 31, and the
%   two-point toy-example special case, Eq. 32). We derive phi(r) for
%   a uniform areal density (Phi(phi,r) = Phi0 = const for d0<=r<=
%   hexagon boundary) from standard hexagonal-cell geometry (e.g.
%   T. S. Rappaport, "Wireless Communications: Principles and
%   Practice," 2nd ed. -- already cited by the paper itself as [23]):
%   for a regular hexagon with circumradius dmax, the polar boundary is
%       R(phi) = apothem / cos( wrap_to_pm_pi_6(phi) ),  apothem = dmax*cos(pi/6)
%   Integrating Phi0 over the available angular range at each radius r
%   gives (Eq. 31):
%       phi(r) = Phi0 * r * 2*pi                          for d0<=r<=apothem
%       phi(r) = Phi0 * r * ( 2*pi - 12*acos(apothem/r) )  for apothem<r<=dmax
%       phi(r) = 0                                          otherwise
%   The overall normalization constant Phi0 cancels in every quantity
%   this module computes (m-quantile partition ranges, Eq. 38, and the
%   bracket-averages, Eq. 30, are both invariant to a uniform rescaling
%   of phi), so Phi0 = 1 is used throughout without loss of generality.
%
% EQUATION REFERENCES
%   Eq. (6)-(8), (28)-(41).

nc = struct();

% ===================================================================
% PART A -- Toy example of areal coverage (Section I-A, Eqs. 6-8)
% ===================================================================

% Eq. (6): total capacity for two ranges d1 (fraction kappa of nodes)
% and d2 (fraction 1-kappa), general (unconstrained) m1, m2.
nc.toy_capacity_general = @(m1, C1_d1, m2, C1_d2) m1.*C1_d1 + m2.*C1_d2;

% Eq. (7): equal-TRANSMIT-POWER toy example. Given d1, d2, kappa,
% gamma, alpha, and Delta_SIR(d2) (via processing_gain), returns
% [Ctot, m1, m2] using C1(d) ~ 1/d^gamma (arbitrary proportionality
% constant K=1; only RATIOS are meaningful, consistent with the
% paper's relative-scale figures).
nc.toy_equal_power = @(d1, d2, kappa, gamma, alpha, delta_sir_d2) ...
    local_toy_equal_power(d1, d2, kappa, gamma, alpha, delta_sir_d2);

% Eq. (8): equal-IpI toy example (same IpI at d1 as at d2, transmit
% power reduced by (d1/d2)^gamma at the shorter range instead).
nc.toy_equal_ipi = @(d1, d2, kappa, gamma, alpha, delta_sir_d2) ...
    local_toy_equal_ipi(d1, d2, kappa, gamma, alpha, delta_sir_d2);

% ===================================================================
% PART B -- Hexagonal-cell uniform areal coverage (Section IV)
% ===================================================================

% Eq. (31) radial node density for a uniform areal density within a
% hexagonal cell of circumradius dmax, restricted to d0 <= r <= dmax
% (Phi0 = 1, see "IMPLEMENTATION ASSUMPTION" above).
nc.hex_radial_density = @(r, d0, dmax) local_hex_phi(r, d0, dmax);

% Eq. (38): partition ranges d_0=d0 < d_1 < ... < d_m=dmax such that
% each partition carries an equal fraction 1/m of the total node mass
% (m-quantiles of phi(r), Sec. IV, citing [22]).
nc.equal_ipi_partition_ranges = @(d0, dmax, m) ...
    local_equal_mass_partition(d0, dmax, m);

% Eq. (30): bracket-average <f(r)>_i of a function f(r) over the range
% interval [d_lo, d_hi], weighted by phi(r).
nc.bracket_average = @(f, d_lo, d_hi) local_bracket_average(f, d_lo, d_hi);

% Eq. (39): inverse-SIR for each PSF channel i = 1..m under the
% equal-IpI partitioning (non-flip PSFs), as a function of the
% partition ranges and path-loss exponent gamma.
nc.equal_ipi_sir_inv = @(d_edges, gamma) local_equal_ipi_sir_inv(d_edges, gamma, 1);

% Eq. (41): inverse-SIR for "flip" PSF-pair partitioning (factor of 2
% relative to Eq. 39, reflecting two co-located, equal-duty-cycle
% channels per range interval).
nc.equal_ipi_sir_inv_flip = @(d_edges, gamma) local_equal_ipi_sir_inv(d_edges, gamma, 2);

% Maximum number of partitions m (or PSF-pairs, for the flip case)
% that satisfies the SIR constraint max_i(SIR_inv_i) <= threshold,
% where threshold = Nm/(alpha*N0dB) (Eq. 39) or the same threshold
% used in Eq. (41) for pairs. Performed by an increasing search over m
% (SIR_inv typically increases with m, i.e. more channels -> more
% mutual interference -> the constraint eventually fails).
nc.solve_max_m = @(d0, dmax, gamma, threshold, use_flip) ...
    local_solve_max_m(d0, dmax, gamma, threshold, use_flip);

% Eq. (46): Clim = C1(Nm) * Nm / (alpha*N0dB), the common asymptotic
% capacity limit (same functional form as processing_gain().Clim but
% kept here too for interface completeness with the paper's Sec. IV-C
% notation).
nc.Clim = @(C1_Nm, Nm, alpha, N0dB) C1_Nm .* Nm ./ (alpha .* N0dB);

% Eq. (44): relative number of equal-power nodes per PSF channel i,
% given the harmonic-mean-normalized d_i^gamma values.
nc.relative_nodes_equal_power = @(d_edges, gamma) local_relative_nodes_equal_power(d_edges, gamma);

% Eq. (45): relative number of equal-IpI nodes per PSF channel (always
% uniform, 1/m, independent of geometry or gamma).
nc.relative_nodes_equal_ipi = @(m) ones(1,m) ./ m;

end

% ===================================================================
% Local implementations
% ===================================================================

function [Ctot, m1, m2] = local_toy_equal_power(d1, d2, kappa, gamma, alpha, delta_sir_d2)
    % Eq. (7)
    m2 = (1-kappa) * (1 + delta_sir_d2/alpha);
    m1 = (d1/d2)^gamma * (kappa/(1-kappa)) * m2;
    C1_d1 = 1/d1^gamma;  % C1(d) ~ 1/d^gamma, arbitrary unit constant
    C1_d2 = 1/d2^gamma;
    Ctot = m1*C1_d1 + m2*C1_d2;
end

function [Ctot, m1, m2] = local_toy_equal_ipi(d1, d2, kappa, gamma, alpha, delta_sir_d2)
    % Eq. (8)
    m2 = (1-kappa) * (1 + delta_sir_d2/alpha);
    m1 = (kappa/(1-kappa)) * m2;
    C1_d2 = 1/d2^gamma;   % same IpI at d1 as d2 -> same C1
    Ctot = (m1+m2) * C1_d2;
end

function phi = local_hex_phi(r, d0, dmax)
    apothem = dmax * cos(pi/6);
    phi = zeros(size(r));
    in_range = (r >= d0) & (r <= dmax);
    inner = in_range & (r <= apothem);
    outer = in_range & (r > apothem);
    phi(inner) = r(inner) * 2*pi;
    ratio = apothem ./ r(outer);
    ratio = min(ratio, 1); % numerical guard
    phi(outer) = r(outer) .* ( 2*pi - 12*acos(ratio) );
end

function d_edges = local_equal_mass_partition(d0, dmax, m)
    % Numerically build the CDF of phi(r) over [d0,dmax] on a fine grid
    % and invert it at the m-quantiles k/m, k=0..m (Eq. 38).
    Ngrid = 20000;
    r = linspace(d0, dmax, Ngrid);
    phi = local_hex_phi(r, d0, dmax);
    cdf = [0, cumsum( (phi(1:end-1)+phi(2:end))/2 .* diff(r) )];
    cdf = cdf / cdf(end);
    d_edges = zeros(1, m+1);
    d_edges(1) = d0;
    d_edges(end) = dmax;
    for i = 1:(m-1)
        target = i/m;
        d_edges(i+1) = interp1(cdf, r, target, 'linear');
    end
end

function avg = local_bracket_average(f, d_lo, d_hi)
    % Generic Eq. (30)-style bracket average of f(r) over [d_lo,d_hi]
    % under a locally-uniform weighting. NOTE: the SIR-margin
    % computations in this module (Eqs. 39/41) do NOT call this
    % generic helper -- they use local_partition_moment(), which
    % correctly weights by the hexagon's GLOBAL radial density
    % phi(r; d0,dmax) restricted to the sub-range. This generic
    % handle is exposed for interface completeness / other simple,
    % geometry-independent uses of Eq. (30).
    Ngrid = 2000;
    r = linspace(d_lo, d_hi, Ngrid);
    num = trapz(r, f(r));
    den = (d_hi - d_lo);
    avg = num/den;
end

function m_val = local_partition_moment(d_lo, d_hi, d0, dmax, gamma, d_target)
    % <(d_target/r)^gamma>_[d_lo,d_hi], weighted by the GLOBAL hexagon
    % density phi(r) restricted to the sub-range [d_lo,d_hi] (Eq. 30).
    Ngrid = 2000;
    r = linspace(d_lo, d_hi, Ngrid);
    phi = local_hex_phi(r, d0, dmax);
    f = (d_target ./ r).^gamma;
    num = trapz(r, f.*phi);
    den = trapz(r, phi);
    m_val = num/den;
end

function sir_inv = local_equal_ipi_sir_inv(d_edges, gamma, flip_factor)
    % Eq. (39) [flip_factor=1] / Eq. (41) [flip_factor=2].
    m = numel(d_edges) - 1;
    d0 = d_edges(1);
    dmax = d_edges(end);
    moments = zeros(1,m); % <(d_j/r)^gamma>_j for j=1..m
    for j = 1:m
        moments(j) = local_partition_moment(d_edges(j), d_edges(j+1), d0, dmax, gamma, d_edges(j+1));
    end
    S = sum(moments);
    sir_inv = zeros(1,m);
    for i = 1:m
        sir_inv(i) = flip_factor * S / moments(i) - 1;
    end
end

function m_max = local_solve_max_m(d0, dmax, gamma, threshold, use_flip)
    flip_factor = 1;
    if use_flip, flip_factor = 2; end
    m_max = 1;
    for m_try = 1:400
        d_edges = local_equal_mass_partition(d0, dmax, m_try);
        sir_inv = local_equal_ipi_sir_inv(d_edges, gamma, flip_factor);
        if max(sir_inv) <= threshold
            m_max = m_try;
        else
            break;
        end
    end
end

function rel = local_relative_nodes_equal_power(d_edges, gamma)
    % Eq. (44): relative # of equal-power nodes per PSF channel,
    % using each partition's outer edge d_i as the representative range.
    m = numel(d_edges)-1;
    di_gamma = d_edges(2:end).^gamma;
    h = helper_functions();
    dbar_gamma = h.harmmean_local(di_gamma);
    rel = (1/m) * dbar_gamma ./ di_gamma;
end
