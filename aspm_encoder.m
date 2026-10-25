function [kj, mj, signs, x_hat, N_total] = aspm_encoder(bits, M, Np, deltaN, n_spacing)
%ASPM_ENCODER M-ASPM pulse-position encoder (Eqs. 9-10).
%
% PURPOSE
%   Maps an information bit sequence to M-ary symbols and, from those,
%   to a "designed" sparse pulse train x_hat[k] whose nonzero,
%   alternating-sign samples mark the pulse arrival times (Eq. 9-10).
%
% INPUTS
%   bits      - 1xN_bits row vector of 0/1 bits (from generate_bits.m).
%               N_bits must be a multiple of log2(M).
%   M         - M-ary alphabet size (number of states per pulse).
%   Np        - average interpulse interval (IpI), in samples.
%   deltaN    - integer offset Delta_N in Eq. (10) (IMPLEMENTATION
%               ASSUMPTION: the paper leaves Delta_N as an arbitrary
%               integer with no stated typical value; we use
%               deltaN = 0, the simplest choice, with no loss of
%               generality since it is a pure time-origin shift).
%   n_spacing - pulse-position spacing "n" for pulse_position_mapper.m
%               (Eq. 27 uses n = 4).
%
% OUTPUTS
%   kj      - 1xN_sym vector of pulse arrival sample indices (Eq. 10).
%   mj      - 1xN_sym vector of symbol indices (1..M) encoded by each pulse.
%   signs   - 1xN_sym vector of (-1)^j alternating signs (Eq. 9), j=0-based.
%   x_hat   - sparse row vector (length N_total) with x_hat[kj]=signs (Eq. 9).
%   N_total - length of x_hat (= last kj + Np, enough room for one more
%             average IpI so pulse-shaping convolution has no edge issues).
%
% THEORY
%   Eq. (9):  x_hat[k] = sum_j [[k=kj]] * (-1)^j
%   Eq. (10): kj = j*Np + Delta_N + Delta_k[mj]
%   Symbol grouping: each pulse encodes log2(M) bits, mj in {1,...,M}
%   (Sec. II, "For mj in {1,2,...,M} this pulse train encodes log2(M)
%   bits per pulse").
%
% PAPER REFERENCE
%   Section II, Eqs. (9)-(10), (27).
%
% EQUATION REFERENCES
%   Eq. (9), (10), (27).

bits_per_symbol = log2(M);
if abs(bits_per_symbol - round(bits_per_symbol)) > 1e-9
    error('aspm_encoder:M', 'M must be a power of two.');
end
bits_per_symbol = round(bits_per_symbol);

if mod(numel(bits), bits_per_symbol) ~= 0
    error('aspm_encoder:bits', 'Number of bits must be a multiple of log2(M).');
end

N_sym = numel(bits) / bits_per_symbol;
bits_matrix = reshape(bits, bits_per_symbol, N_sym); % column j = bits for symbol j

% Convert each column of bits (MSB first) to an integer 0..M-1, then to
% a 1-indexed symbol mj = value+1, matching Eq. (27)'s convention
% (8*a+4*b+2*c+d for a 4-bit symbol, i.e. MSB-first binary value).
weights = 2.^(bits_per_symbol-1:-1:0)';
symbol_values = weights' * bits_matrix;   % 1xN_sym, values 0..M-1
mj = symbol_values + 1;                   % 1xN_sym, values 1..M

ppm = pulse_position_mapper();
delta_k = ppm.forward(mj, n_spacing);     % Eq. (27)-style offset

j_index = 0:(N_sym-1);
kj = j_index .* Np + deltaN + delta_k;    % Eq. (10)
signs = (-1).^j_index;                    % Eq. (9)

N_total = kj(end) + Np;                   % headroom for PSF convolution
x_hat = aggregate_pulse_generator(kj, signs, N_total);

end
